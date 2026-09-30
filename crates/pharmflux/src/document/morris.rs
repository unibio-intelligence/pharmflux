use super::*;
use pharmflux_core::{identity::canonical_hash, morris::*};
impl CompiledDocument {
    pub fn execute_morris(
        self: &Arc<Self>,
        request: &MorrisExecutionRequest,
    ) -> Result<MorrisExecutionResult, Error> {
        let (mut result, design) = match &request.problem {
            MorrisProblem::Explicit(r) => (self.morris(r)?, None),
            MorrisProblem::Generated(r) => {
                let generated = generate_morris_design(r)?;
                (
                    self.morris(&generated.request)?,
                    Some(MorrisDesignProvenance {
                        generator: generated.generator,
                        seed: generated.seed,
                        design_request_hash: generated.design_request_hash,
                    }),
                )
            }
        };
        result.analysis_request_hash = canonical_hash(
            &serde_json::json!({"request":request,"analysis_request_hash":result.analysis_request_hash}),
        )?;
        Ok(MorrisExecutionResult {
            schema: MorrisResultSchema::V01,
            result,
            design,
        })
    }

    pub fn morris(self: &Arc<Self>, request: &MorrisRequest) -> Result<MorrisResult, Error> {
        let d = request.ranges.len();
        let count = request.scan.points.len();
        if d == 0
            || d > 11
            || request.trajectories.len() < 2
            || request.trajectories.len() > 12
            || count > 25
            || request.num_levels < 4
            || request.num_levels > 100
            || request.num_levels % 2 != 0
        {
            return Err(invalid("Morris requires 1..11 parameters, 2..12 trajectories, at most 25 points and an even 4..100-level grid"));
        }
        let output = self
            .output_names
            .iter()
            .position(|n| n == &request.output)
            .ok_or_else(|| invalid("unknown Morris output"))?;
        let mut names = BTreeSet::new();
        let mut bounds = Vec::new();
        for r in &request.ranges {
            if !names.insert(r.parameter.clone())
                || !self.parameters.iter().any(|p| p.name == r.parameter)
            {
                return Err(invalid(
                    "Morris range names must be unique model parameters",
                ));
            }
            let unit = Unit::parse(&r.lower.unit)?;
            let lo = quantity(&r.lower)?.in_unit(unit)?;
            let hi = quantity(&r.upper)?.in_unit(unit)?;
            if !lo.is_finite() || !hi.is_finite() || !(hi - lo).is_finite() || hi <= lo {
                return Err(invalid("Morris ranges require finite increasing bounds"));
            }
            bounds.push((lo, hi - lo, unit));
        }
        let mut points = BTreeMap::new();
        let levels = (request.num_levels - 1) as f64;
        let delta = request.num_levels as f64 / (2. * levels);
        for (index, p) in request.scan.points.iter().enumerate() {
            if p.parameters.keys().cloned().collect::<BTreeSet<_>>() != names {
                return Err(invalid(
                    "each Morris point must specify exactly the selected parameters",
                ));
            }
            let mut x = Vec::new();
            for (r, (lo, span, unit)) in request.ranges.iter().zip(&bounds) {
                let value = (quantity(&p.parameters[&r.parameter])?.in_unit(*unit)? - lo) / span;
                if !value.is_finite()
                    || !(-1e-12..=1. + 1e-12).contains(&value)
                    || (value * levels - (value * levels).round()).abs() > 1e-8
                {
                    return Err(invalid("Morris point is outside its normalized grid"));
                }
                x.push((value * levels).round() / levels);
            }
            if points.insert(p.id.clone(), (index, x)).is_some() {
                return Err(invalid("duplicate Morris point ID"));
            }
        }
        let mut used = BTreeSet::new();
        let mut steps = Vec::new();
        for path in &request.trajectories {
            if path.len() != d + 1 {
                return Err(invalid("Morris trajectory must contain d+1 points"));
            }
            for id in path {
                if !points.contains_key(id) || !used.insert(id) {
                    return Err(invalid(
                        "Morris trajectories require distinct known point IDs",
                    ));
                }
            }
            let mut changed = BTreeSet::new();
            for pair in path.windows(2) {
                let (a, x) = &points[&pair[0]];
                let (b, y) = &points[&pair[1]];
                let changes: Vec<_> = x
                    .iter()
                    .zip(y)
                    .enumerate()
                    .filter_map(|(i, (x, y))| {
                        if (y - x).abs() > 1e-12 {
                            Some((i, y - x))
                        } else {
                            None
                        }
                    })
                    .collect();
                if changes.len() != 1 {
                    return Err(invalid("Morris step must change exactly one parameter"));
                }
                let (p, step) = changes[0];
                if (step.abs() - delta).abs() > 1e-8 || !changed.insert(p) {
                    return Err(invalid(
                        "Morris step size or parameter permutation is invalid",
                    ));
                }
                steps.push((*a, *b, p, step));
            }
        }
        if used.len() != count {
            return Err(invalid(
                "every scan point must belong to a Morris trajectory",
            ));
        }
        let scan = self.scan(&request.scan)?;
        let first = &scan.results[0].result;
        let mut values = Vec::new();
        for p in &scan.results {
            if p.result.times != first.times || p.result.sides != first.sides {
                return Err(invalid(
                    "Morris row requires identical time/side grids across points",
                ));
            }
            values.push(
                *p.result.outputs[output]
                    .values
                    .get(request.row)
                    .ok_or_else(|| invalid("Morris output row is out of range"))?,
            );
        }
        let mut effects = vec![Vec::new(); d];
        for (a, b, p, delta) in steps {
            let e = (values[b] - values[a]) / delta;
            if !e.is_finite() {
                return Err(Error::new(
                    ErrorCode::Domain,
                    "Morris elementary effect overflow",
                ));
            }
            effects[p].push(e);
        }
        let mut summaries = Vec::new();
        for (range, effects) in request.ranges.iter().zip(effects) {
            let mut mu = 0.;
            let mut m2 = 0.;
            let mut mu_star = 0.;
            for (i, e) in effects.iter().enumerate() {
                let n = (i + 1) as f64;
                let difference = e - mu;
                mu += difference / n;
                m2 += difference * (e - mu);
                mu_star += (e.abs() - mu_star) / n;
            }
            let sigma = libm::sqrt(m2 / (effects.len() - 1) as f64);
            if ![mu, mu_star, sigma].iter().all(|v| v.is_finite()) {
                return Err(Error::new(ErrorCode::Domain, "Morris summary overflow"));
            }
            summaries.push(ElementaryEffects {
                parameter: range.parameter.clone(),
                unit: self.output_units[output].clone(),
                effects,
                mu,
                mu_star,
                sigma,
            });
        }
        Ok(MorrisResult {
            analysis_request_hash: canonical_hash(
                &serde_json::json!({"algorithm":"morris_explicit_grid_v1","request":request}),
            )?,
            scan,
            effects: summaries,
        })
    }
}
