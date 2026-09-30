//! Pharmacology execution engine. Experimental M0 API.
pub mod capabilities;
pub mod compile;
pub mod document;
pub mod language;
pub mod linear_pk;
pub mod runtime;
pub mod steady_state;
use diffsol::{NalgebraLU, NalgebraMat, OdeBuilder, OdeSolverMethod};
use pharmflux_core::run::Solver;
pub use pharmflux_core::{Budgets, Error, ErrorCode, Event, Protocol, Sample};
use std::cell::{Cell, RefCell};

/// Exact vendored Diffsol source revision, including both local Rosenbrock methods.
pub(crate) const DIFFSOL_BACKEND_VERSION: &str =
    concat!("0.17.1+sha256:", env!("PHARMFLUX_DIFFSOL_SOURCE_SHA256"));

/// Constant-coefficient interval propagation used by the common event driver.
/// The driver retains ownership of event ordering, observations and budgets.
pub trait IntervalPropagator {
    fn propagate(&self, amounts: &[f64], rates: &[f64], duration: f64) -> Result<Vec<f64>, Error>;
}

/// Models supply equations and directional Jacobian products; the event driver
/// supplies rates in the model state/time basis.
pub trait Model {
    /// Opt in only when this propagator describes the complete current model.
    fn interval_propagator(&self) -> Option<&dyn IntervalPropagator> {
        None
    }

    fn check_delivery(&self, _event: &Event) -> Result<(), Error> {
        Ok(())
    }
    fn validate(&self) -> Result<(), Error>;
    /// Binding discontinuities are fixed in advance and never inferred from solver callbacks.
    fn discontinuities(&self) -> &[f64] {
        &[]
    }
    /// Times where the ODE right-hand side changes. A model may have additional
    /// smooth checkpoints in `discontinuities` that still split integration.
    fn rhs_jump_times(&self) -> &[f64] {
        self.discontinuities()
    }
    /// Called after the pre-event state is captured, before resets and doses.
    fn enter_boundary(&self, _time: f64) -> Result<(), Error> {
        Ok(())
    }
    fn initial(&self) -> Vec<f64>;
    /// Validate a physical state after accepted steps, dense output and discrete
    /// updates. This deliberately excludes Newton iterates and rejected steps.
    fn check_state(&self, _time: f64, _state: &[f64]) -> Result<(), Error> {
        Ok(())
    }

    fn dose_scale(&self, state: usize) -> Option<f64>;
    /// Apply a fixed-time bolus, including any augmented-state jump derivatives.
    fn apply_bolus(
        &self,
        _time: f64,
        target: usize,
        amount: f64,
        state: &mut [f64],
    ) -> Result<(), Error> {
        let scale = self
            .dose_scale(target)
            .ok_or_else(|| Error::new(ErrorCode::InvalidInput, "unsupported bolus target"))?;
        *state
            .get_mut(target)
            .ok_or_else(|| Error::new(ErrorCode::InvalidInput, "invalid bolus target"))? +=
            amount * scale;
        Ok(())
    }
    /// Constant state resets have zero parameter derivative for the reset state.
    fn apply_reset(
        &self,
        _time: f64,
        target: usize,
        value: f64,
        state: &mut [f64],
    ) -> Result<(), Error> {
        *state
            .get_mut(target)
            .ok_or_else(|| Error::new(ErrorCode::InvalidInput, "invalid reset target"))? = value;
        Ok(())
    }
    /// Receive the unscaled input rate so zero dose scales retain their derivatives.
    fn accumulate_infusion_rate(
        &self,
        _time: f64,
        target: usize,
        rate: f64,
        rates: &mut [f64],
    ) -> Result<(), Error> {
        let scale = self
            .dose_scale(target)
            .ok_or_else(|| Error::new(ErrorCode::InvalidInput, "unsupported infusion target"))?;
        *rates
            .get_mut(target)
            .ok_or_else(|| Error::new(ErrorCode::InvalidInput, "invalid infusion target"))? +=
            rate * scale;
        Ok(())
    }
    fn rhs(&self, time: f64, x: &[f64], rate: &[f64], out: &mut [f64]) -> Result<(), Error>;
    fn jac_mul(
        &self,
        time: f64,
        x: &[f64],
        rates: &[f64],
        v: &[f64],
        out: &mut [f64],
    ) -> Result<(), Error>;
}

fn sorted_unique(values: &mut Vec<f64>) {
    values.sort_by(f64::total_cmp);
    values.dedup();
}

fn incoming_f64(time: f64) -> f64 {
    if time == 0.0 {
        -f64::from_bits(1)
    } else if time > 0.0 {
        f64::from_bits(time.to_bits() - 1)
    } else {
        f64::from_bits(time.to_bits() + 1)
    }
}

trait ProtocolExt {
    fn validate<M: Model>(&self, model: &M) -> Result<(), Error>;
    fn boundaries(&self) -> Vec<f64>;
    fn has_event(&self, t: f64) -> bool;
}
impl ProtocolExt for Protocol {
    fn validate<M: Model>(&self, model: &M) -> Result<(), Error> {
        model.validate()?;
        if !self.start.is_finite()
            || !self.end.is_finite()
            || self.end <= self.start
            || !self.rtol.is_finite()
            || !self.atol.is_finite()
            || self.rtol <= 0.0
            || self.atol <= 0.0
        {
            return Err("end and tolerances must be finite and positive".into());
        }
        let initial = model.initial();
        model.check_state(self.start, &initial)?;
        if initial.is_empty() || initial.iter().any(|v| !v.is_finite()) {
            return Err("invalid initial state".into());
        }
        if self
            .samples
            .iter()
            .any(|t| !t.is_finite() || *t < self.start || *t > self.end)
        {
            return Err("sample outside simulation window".into());
        }
        let mut reset_orders = std::collections::BTreeSet::new();
        for e in &self.events {
            if !e.time().is_finite()
                || e.time() < self.start
                || e.time() > self.end
                || e.target() >= initial.len()
            {
                return Err("invalid dose".into());
            }
            if let Event::Reset {
                value, time, order, ..
            } = e
            {
                if !value.is_finite() {
                    return Err("reset value must be finite".into());
                }
                if !reset_orders.insert((if *time == 0.0 { 0 } else { time.to_bits() }, *order)) {
                    return Err("simultaneous resets require distinct order values".into());
                }
                continue;
            }
            model.check_delivery(e)?;
            let amount = e.amount().expect("dose variant");
            if !amount.is_finite() || amount < 0.0 {
                return Err("invalid dose amount".into());
            }
            let scale = model
                .dose_scale(e.target())
                .ok_or("unsupported dose target")?;
            if !scale.is_finite() || scale < 0.0 || !(amount * scale).is_finite() {
                return Err("invalid dose scale or overflow".into());
            }
            if let Event::Infusion {
                time,
                amount,
                duration,
                ..
            } = e
            {
                if !duration.is_finite()
                    || *duration <= 0.0
                    || time + duration > self.end
                    || time + duration <= *time
                    || !(amount / duration * scale).is_finite()
                {
                    return Err("invalid infusion duration or rate".into());
                }
            }
        }
        Ok(())
    }
    fn boundaries(&self) -> Vec<f64> {
        let mut times = vec![self.start, self.end];
        for e in &self.events {
            times.push(e.time());
            if let Event::Infusion { time, duration, .. } = e {
                times.push(time + duration);
            }
        }
        sorted_unique(&mut times);
        times
    }
    fn has_event(&self, t: f64) -> bool {
        self.events.iter().any(|e| {
            e.time() == t || matches!(e, Event::Infusion{time,duration,..} if time+duration == t)
        })
    }
}

/// Segment only at true dosing discontinuities, never at ordinary sample times.
/// Every boundary has explicit pre/post output; observations elsewhere are post.
pub fn simulate<M: Model>(model: &M, protocol: &Protocol) -> Result<Vec<Sample>, Error> {
    simulate_with_budgets(model, protocol, Budgets::default())
}

pub fn simulate_with_budgets<M: Model>(
    model: &M,
    protocol: &Protocol,
    budgets: Budgets,
) -> Result<Vec<Sample>, Error> {
    simulate_with_absolute_tolerances(model, protocol, budgets, None)
}

pub(crate) fn simulate_with_absolute_tolerances<M: Model>(
    model: &M,
    protocol: &Protocol,
    budgets: Budgets,
    absolute_tolerances: Option<&[f64]>,
) -> Result<Vec<Sample>, Error> {
    simulate_with_solver(
        model,
        protocol,
        budgets,
        absolute_tolerances,
        Solver::DiffsolBdf,
    )
}

pub(crate) fn simulate_with_solver<M: Model>(
    model: &M,
    protocol: &Protocol,
    budgets: Budgets,
    absolute_tolerances: Option<&[f64]>,
    solver_choice: Solver,
) -> Result<Vec<Sample>, Error> {
    let event_count = protocol
        .events
        .len()
        .checked_add(model.discontinuities().len())
        .ok_or_else(|| Error::new(ErrorCode::WorkBudget, "event count overflow"))?;
    if event_count > budgets.events {
        return Err(Error::new(
            ErrorCode::WorkBudget,
            "event count exceeds budget",
        ));
    }
    if model
        .discontinuities()
        .iter()
        .any(|t| !t.is_finite() || *t < protocol.start || *t > protocol.end)
        || model.discontinuities().windows(2).any(|w| w[0] >= w[1])
    {
        return Err("invalid binding discontinuities".into());
    }
    protocol.validate(model)?;
    let mut state = model.initial();
    let n = state.len();
    let scalar = [protocol.atol];
    let tolerances = absolute_tolerances.unwrap_or(&scalar);
    if (absolute_tolerances.is_some() && tolerances.len() != n)
        || tolerances.iter().any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err("absolute tolerances must be finite, positive and match every state".into());
    }
    // Each event adds at most two boundaries, with pre/post rows at each.
    let row_limit = event_count
        .checked_mul(4)
        .and_then(|v| v.checked_add(protocol.samples.len()))
        .and_then(|v| v.checked_add(2))
        .and_then(|v| v.checked_mul(n));
    if row_limit.is_none_or(|v| v > budgets.output_values) {
        return Err(Error::new(
            ErrorCode::OutputBudget,
            "requested output exceeds budget",
        ));
    }
    let callbacks = Cell::new(0u64);
    let first_error = RefCell::new(None::<Error>);
    let guarded =
        |time: f64, out: &mut [f64], evaluate: &dyn Fn(&mut [f64]) -> Result<(), Error>| {
            let count = callbacks.get().saturating_add(1);
            callbacks.set(count);
            if first_error.borrow().is_none() {
                let result = if count > budgets.solver_callbacks {
                    Err(Error::new(
                        ErrorCode::WorkBudget,
                        "solver callback budget exhausted",
                    ))
                } else {
                    evaluate(out).and_then(|()| {
                        if out.iter().any(|v| !v.is_finite()) {
                            Err(Error::new(ErrorCode::Domain, "non-finite model evaluation"))
                        } else {
                            Ok(())
                        }
                    })
                };
                if let Err(mut error) = result {
                    error.time = Some(time);
                    *first_error.borrow_mut() = Some(error);
                }
            }
            if first_error.borrow().is_some() {
                out.fill(f64::NAN);
            }
        };
    let solver_error = |error: diffsol::DiffsolError| {
        first_error
            .borrow()
            .clone()
            .unwrap_or_else(|| Error::new(ErrorCode::Solver, error.to_string()))
    };
    let mut boundaries = protocol.boundaries();
    boundaries.extend_from_slice(model.discontinuities());
    sorted_unique(&mut boundaries);
    let mut samples = protocol.samples.clone();
    sorted_unique(&mut samples);
    let mut output = Vec::new();
    let mut stopped_infusions = std::collections::BTreeSet::new();
    for (index, &t) in boundaries.iter().enumerate() {
        model.check_state(t, &state)?;
        if protocol.has_event(t) || model.discontinuities().contains(&t) {
            output.push(Sample {
                time: t,
                side: "pre".into(),
                values: state.clone(),
            });
        }
        model.enter_boundary(t)?;
        let mut resets: Vec<_> = protocol
            .events
            .iter()
            .filter_map(|e| match e {
                Event::Reset {
                    time,
                    target,
                    value,
                    order,
                    active_inputs,
                } if *time == t => Some((*order, *target, *value, *active_inputs)),
                _ => None,
            })
            .collect();
        resets.sort_by_key(|(order, ..)| *order);
        for (_, target, value, active_inputs) in resets {
            model.apply_reset(t, target, value, &mut state)?;
            model.check_state(t, &state)?;
            if active_inputs == pharmflux_core::ActiveInputs::StopTarget {
                for (id, event) in protocol.events.iter().enumerate() {
                    if matches!(event,Event::Infusion{time,target:input_target,duration,..} if *input_target==target&&*time<t&&t<*time+*duration)
                    {
                        stopped_infusions.insert(id);
                    }
                }
            }
        }
        for e in &protocol.events {
            if let Event::Bolus {
                time,
                target,
                amount,
            } = e
            {
                if *time == t {
                    model.apply_bolus(t, *target, *amount, &mut state)?;
                    model.check_state(t, &state)?;
                }
            }
        }
        if state.iter().any(|x| !x.is_finite()) {
            return Err("combined bolus overflow".into());
        }
        output.push(Sample {
            time: t,
            side: "post".into(),
            values: state.clone(),
        });
        if index + 1 == boundaries.len() {
            break;
        }
        let end = boundaries[index + 1];
        let mut rates = vec![0.0; n];
        for (id, e) in protocol.events.iter().enumerate() {
            if stopped_infusions.contains(&id) {
                continue;
            }
            if let Event::Infusion {
                time,
                target,
                amount,
                duration,
            } = e
            {
                if *time <= t && t < time + duration {
                    model.accumulate_infusion_rate(t, *target, amount / duration, &mut rates)?;
                }
            }
        }
        if rates.iter().chain(state.iter()).any(|x| !x.is_finite()) {
            return Err("combined event overflow".into());
        }
        let mut eval: Vec<f64> = samples
            .iter()
            .copied()
            .filter(|x| *x > t && *x < end)
            .collect();
        eval.push(end);
        if let Some(propagator) = model.interval_propagator() {
            // Every sample is evaluated from the same interval origin; adding
            // observations must not alter subsequent trajectory values.
            let initial = state.clone();
            for time in eval {
                let count = callbacks.get().saturating_add(1);
                callbacks.set(count);
                if count > budgets.solver_callbacks {
                    let mut error = Error::new(
                        ErrorCode::WorkBudget,
                        "analytic propagation budget exhausted",
                    );
                    error.time = Some(time);
                    return Err(error);
                }
                let values =
                    propagator
                        .propagate(&initial, &rates, time - t)
                        .map_err(|mut error| {
                            error.time = Some(time);
                            error
                        })?;
                if values.len() != n || values.iter().any(|v| !v.is_finite()) {
                    let mut error =
                        Error::new(ErrorCode::Domain, "invalid analytic propagation result");
                    error.time = Some(time);
                    return Err(error);
                }
                model.check_state(time, &values)?;
                if time == end {
                    state = values;
                } else {
                    output.push(Sample {
                        time,
                        side: "post".into(),
                        values,
                    });
                }
            }
            continue;
        }
        let initial = state.clone();
        // Integrate elapsed interval time so tiny accepted steps can advance
        // even when the scientific clock has a large positive/negative origin.
        // Model callbacks and all reported observations retain absolute time.
        let interval_duration = end - t;
        if !interval_duration.is_finite() || interval_duration <= 0.0 {
            return Err(Error::new(
                ErrorCode::InvalidInput,
                "invalid interval duration",
            ));
        }
        let is_rosenbrock = matches!(
            solver_choice,
            Solver::DiffsolRosenbrock23 | Solver::DiffsolRodas5p
        );
        let stop_is_discontinuity =
            protocol.has_event(end) || model.rhs_jump_times().contains(&end);
        let previous_rhs_time = Cell::new(None::<(f64, f64)>);
        let problem = OdeBuilder::<NalgebraMat<f64>>::new()
            .t0(0.0)
            .rtol(protocol.rtol)
            .atol(tolerances.iter().copied())
            .rhs_implicit(
                |x, _p, elapsed, y| {
                    let mut time = t + elapsed;
                    if is_rosenbrock && first_error.borrow().is_none() {
                        let collapsed_probe = previous_rhs_time
                            .get()
                            .is_some_and(|(previous_elapsed, previous_absolute)| {
                                let local_gap = (previous_elapsed - elapsed).abs();
                                let ulp_margin = 8.0
                                    * f64::EPSILON
                                    * (1.0 + previous_elapsed.abs().max(elapsed.abs()));
                                local_gap > ulp_margin && previous_absolute == time
                            });
                        let rounded_into_jump =
                            stop_is_discontinuity && elapsed < interval_duration && time >= end;
                        if rounded_into_jump {
                            let incoming = incoming_f64(end);
                            let resolution_limit =
                                f64::EPSILON.sqrt() * (1.0 + interval_duration);
                            if incoming > t && end - incoming <= resolution_limit {
                                time = incoming;
                            }
                        }
                        let unresolved_jump = rounded_into_jump && time >= end;
                        if collapsed_probe || unresolved_jump {
                            let mut error = Error::new(
                                ErrorCode::Unsupported,
                                "Rosenbrock time probes are not representable on the absolute model clock",
                            );
                            error.time = Some(time);
                            *first_error.borrow_mut() = Some(error);
                        }
                        previous_rhs_time.set(Some((elapsed, time)));
                    }
                    guarded(time, y, &|out| {
                        model.rhs(time, x, &rates, out)
                    })
                },
                |x, _p, elapsed, v, y| {
                    let time = t + elapsed;
                    guarded(time, y, &|out| {
                        model.jac_mul(time, x, &rates, v, out)
                    })
                },
            )
            .init(
                move |_p, _time, y| {
                    for i in 0..n {
                        y[i] = initial[i];
                    }
                },
                n,
            )
            .build()
            .map_err(&solver_error)?;
        if let Some(error) = first_error.borrow().clone() {
            return Err(error);
        }
        // All methods share the same event driver, dense observation handling,
        // callback budget, and validation of physical states. Only the solver
        // construction and nonlinear failure control differ.
        macro_rules! integrate_interval {
            ($constructor:expr, $configure:expr, $stop:expr) => {{
                let mut solver = $constructor.map_err(&solver_error)?;
                if let Some(error) = first_error.borrow().clone() {
                    return Err(error);
                }
                // Permit any representable step; progress and callback guards
                // below bound the solve when tolerances are very small.
                solver.config_mut().minimum_timestep = 0.0;
                $configure(&mut solver);
                $stop(&mut solver).map_err(|source| {
                    let mut error = solver_error(source);
                    if error.time.is_none() {
                        error.time = Some(t);
                    }
                    error
                })?;
                let mut column = 0;
                while column < eval.len() {
                    let previous_time = solver.state().t;
                    solver.step().map_err(|source| {
                        let mut error = solver_error(source);
                        if error.time.is_none() {
                            error.time = Some(t + previous_time);
                        }
                        error
                    })?;
                    if let Some(error) = first_error.borrow().clone() {
                        return Err(error);
                    }
                    let reached = solver.state().t;
                    model.check_state(t + reached, solver.state().y.batch_as_slice(0))?;
                    if !reached.is_finite() || reached <= previous_time {
                        let mut error = Error::new(
                            ErrorCode::Solver,
                            format!("solver made no forward progress: elapsed_previous={previous_time}, elapsed_reached={reached}, next_step={}, interval_origin={t}, interval_end={end}", solver.state().h),
                        );
                        error.time = Some(t + previous_time);
                        return Err(error);
                    }
                    while column < eval.len() && eval[column] - t <= reached {
                        let time = eval[column];
                        let values = solver
                            .interpolate(time - t)
                            .map_err(&solver_error)?
                            .batch_as_slice(0)
                            .to_vec();
                        model.check_state(time, &values)?;
                        column += 1;
                        if values.iter().any(|x| !x.is_finite()) {
                            return Err("nonfinite solve".into());
                        }
                        if time == end {
                            state = values;
                        } else {
                            output.push(Sample {
                                time,
                                side: "post".into(),
                                values,
                            });
                        }
                    }
                }
            }};
        }
        match solver_choice {
            Solver::DiffsolBdf => integrate_interval!(
                problem.bdf::<NalgebraLU<f64>>(),
                |solver: &mut diffsol::Bdf<_, _>| {
                    solver.config_mut().maximum_newton_fails =
                        budgets.solver_callbacks.min(2_000) as usize;
                },
                |solver: &mut diffsol::Bdf<_, _>| solver.set_stop_time(interval_duration)
            ),
            Solver::DiffsolEsdirk34 => integrate_interval!(
                problem.esdirk34::<NalgebraLU<f64>>(),
                |solver: &mut diffsol::Sdirk<_, _, _>| {
                    solver.config_mut().maximum_newton_fails =
                        budgets.solver_callbacks.min(2_000) as usize;
                },
                |solver: &mut diffsol::Sdirk<_, _, _>| solver.set_stop_time(interval_duration)
            ),
            Solver::DiffsolTrBdf2 => integrate_interval!(
                problem.tr_bdf2::<NalgebraLU<f64>>(),
                |solver: &mut diffsol::Sdirk<_, _, _>| {
                    solver.config_mut().maximum_newton_fails =
                        budgets.solver_callbacks.min(2_000) as usize;
                },
                |solver: &mut diffsol::Sdirk<_, _, _>| solver.set_stop_time(interval_duration)
            ),
            Solver::DiffsolTsit45 => integrate_interval!(
                problem.tsit45(),
                |_solver: &mut _| {},
                |solver: &mut diffsol::ExplicitRk<_, _>| solver.set_stop_time(interval_duration)
            ),
            Solver::DiffsolRosenbrock23 => integrate_interval!(
                problem.rosenbrock23::<NalgebraLU<f64>>(),
                |_solver: &mut _| {},
                |solver: &mut diffsol::Rosenbrock23<_, _>| {
                    if stop_is_discontinuity {
                        solver.set_discontinuity_stop_time(interval_duration)
                    } else {
                        solver.set_stop_time(interval_duration)
                    }
                }
            ),
            Solver::DiffsolRodas5p => integrate_interval!(
                problem.rodas5p::<NalgebraLU<f64>>(),
                |_solver: &mut _| {},
                |solver: &mut diffsol::Rodas5P<_, _>| {
                    if stop_is_discontinuity {
                        solver.set_discontinuity_stop_time(interval_duration)
                    } else {
                        solver.set_stop_time(interval_duration)
                    }
                }
            ),
            Solver::AnalyticLinearPk => {
                return Err(Error::new(
                    ErrorCode::Unsupported,
                    "analytic_linear_pk requires a linear_pk declaration",
                ));
            }
        }
    }
    Ok(output)
}

#[cfg(test)]
mod rosenbrock_integration_tests {
    use super::*;
    use std::cell::{Cell, RefCell};

    struct TimeForcing {
        origin: f64,
    }

    impl Model for TimeForcing {
        fn validate(&self) -> Result<(), Error> {
            Ok(())
        }
        fn initial(&self) -> Vec<f64> {
            vec![0.0]
        }
        fn dose_scale(&self, _: usize) -> Option<f64> {
            None
        }
        fn rhs(&self, time: f64, _: &[f64], _: &[f64], out: &mut [f64]) -> Result<(), Error> {
            out[0] = time - self.origin;
            Ok(())
        }
        fn jac_mul(
            &self,
            _: f64,
            _: &[f64],
            _: &[f64],
            _: &[f64],
            out: &mut [f64],
        ) -> Result<(), Error> {
            out[0] = 0.0;
            Ok(())
        }
    }

    struct TimeJump;

    struct NonlinearDecay {
        jacobian_times: RefCell<Vec<f64>>,
    }

    impl Model for NonlinearDecay {
        fn validate(&self) -> Result<(), Error> {
            Ok(())
        }
        fn initial(&self) -> Vec<f64> {
            vec![2.0]
        }
        fn dose_scale(&self, _: usize) -> Option<f64> {
            None
        }
        fn rhs(&self, _: f64, x: &[f64], _: &[f64], out: &mut [f64]) -> Result<(), Error> {
            out[0] = -x[0] * x[0];
            Ok(())
        }
        fn jac_mul(
            &self,
            time: f64,
            x: &[f64],
            _: &[f64],
            v: &[f64],
            out: &mut [f64],
        ) -> Result<(), Error> {
            self.jacobian_times.borrow_mut().push(time);
            out[0] = -2.0 * x[0] * v[0];
            Ok(())
        }
    }

    impl Model for TimeJump {
        fn validate(&self) -> Result<(), Error> {
            Ok(())
        }
        fn discontinuities(&self) -> &[f64] {
            &[1.0]
        }
        fn initial(&self) -> Vec<f64> {
            vec![0.0]
        }
        fn dose_scale(&self, _: usize) -> Option<f64> {
            Some(1.0)
        }
        fn rhs(&self, time: f64, _: &[f64], _: &[f64], out: &mut [f64]) -> Result<(), Error> {
            out[0] = if time < 1.0 { 0.0 } else { 1.0 };
            Ok(())
        }
        fn jac_mul(
            &self,
            _: f64,
            _: &[f64],
            _: &[f64],
            _: &[f64],
            out: &mut [f64],
        ) -> Result<(), Error> {
            out[0] = 0.0;
            Ok(())
        }
    }

    struct EndpointProbe {
        saw_exact_end: Cell<bool>,
        checkpoint: bool,
    }

    impl Model for EndpointProbe {
        fn validate(&self) -> Result<(), Error> {
            Ok(())
        }
        fn discontinuities(&self) -> &[f64] {
            if self.checkpoint {
                &[1.0]
            } else {
                &[]
            }
        }
        fn rhs_jump_times(&self) -> &[f64] {
            &[]
        }
        fn initial(&self) -> Vec<f64> {
            vec![0.0]
        }
        fn dose_scale(&self, _: usize) -> Option<f64> {
            None
        }
        fn rhs(&self, time: f64, _: &[f64], _: &[f64], out: &mut [f64]) -> Result<(), Error> {
            if time == 1.0 {
                self.saw_exact_end.set(true);
            }
            out[0] = 1.0;
            Ok(())
        }
        fn jac_mul(
            &self,
            _: f64,
            _: &[f64],
            _: &[f64],
            _: &[f64],
            out: &mut [f64],
        ) -> Result<(), Error> {
            out[0] = 0.0;
            Ok(())
        }
    }

    fn protocol(start: f64, end: f64, samples: Vec<f64>) -> Protocol {
        Protocol {
            start,
            end,
            samples,
            events: vec![],
            rtol: 1e-8,
            atol: 1e-10,
        }
    }

    #[test]
    fn smooth_absolute_time_forcing_and_declared_jump() {
        for solver in [Solver::DiffsolRosenbrock23, Solver::DiffsolRodas5p] {
            let smooth = simulate_with_solver(
                &TimeForcing { origin: 2.0 },
                &protocol(2.0, 6.0, vec![2.0, 6.0]),
                Budgets::default(),
                None,
                solver,
            )
            .unwrap();
            assert!((smooth.last().unwrap().values[0] - 8.0).abs() < 1e-6);

            let covariate_only_jump = simulate_with_solver(
                &TimeJump,
                &protocol(0.0, 2.0, vec![1.0, 2.0]),
                Budgets::default(),
                None,
                solver,
            )
            .unwrap();
            let before_change = covariate_only_jump
                .iter()
                .find(|row| row.time == 1.0 && row.side == "pre")
                .unwrap();
            let after_change = covariate_only_jump
                .iter()
                .find(|row| row.time == 1.0 && row.side == "post")
                .unwrap();
            assert!(before_change.values[0].abs() < 1e-8);
            assert!(after_change.values[0].abs() < 1e-8);
            assert!((covariate_only_jump.last().unwrap().values[0] - 1.0).abs() < 1e-6);

            let jumped = simulate_with_solver(
                &TimeJump,
                &Protocol {
                    events: vec![Event::Bolus {
                        time: 1.0,
                        target: 0,
                        amount: 2.0,
                    }],
                    ..protocol(0.0, 2.0, vec![0.5, 1.0, 1.5, 2.0])
                },
                Budgets::default(),
                None,
                solver,
            )
            .unwrap();
            let pre = jumped
                .iter()
                .find(|row| row.time == 1.0 && row.side == "pre")
                .unwrap();
            let post = jumped
                .iter()
                .find(|row| row.time == 1.0 && row.side == "post")
                .unwrap();
            assert!(pre.values[0].abs() < 1e-8, "{solver:?}: {pre:?}");
            assert!((post.values[0] - 2.0).abs() < 1e-8, "{solver:?}: {post:?}");
            let final_amount = jumped.last().unwrap().values[0];
            assert!(
                (final_amount - 3.0).abs() < 1e-6,
                "{solver:?}: {final_amount}"
            );
        }
    }

    #[test]
    fn ordinary_final_target_evaluates_rhs_at_exact_endpoint() {
        for solver in [Solver::DiffsolRosenbrock23, Solver::DiffsolRodas5p] {
            for checkpoint in [false, true] {
                let model = EndpointProbe {
                    saw_exact_end: Cell::new(false),
                    checkpoint,
                };
                let rows = simulate_with_solver(
                    &model,
                    &protocol(0.0, 1.0, vec![1.0]),
                    Budgets::default(),
                    None,
                    solver,
                )
                .unwrap();
                assert!(
                    model.saw_exact_end.get(),
                    "{solver:?} skipped the smooth endpoint with checkpoint={checkpoint}"
                );
                assert!((rows.last().unwrap().values[0] - 1.0).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn nonrepresentable_absolute_time_probes_fail_explicitly() {
        for origin in [-1e16, 1e16] {
            for solver in [Solver::DiffsolRosenbrock23, Solver::DiffsolRodas5p] {
                let result = simulate_with_solver(
                    &TimeForcing { origin },
                    &protocol(origin, origin + 4.0, vec![origin, origin + 4.0]),
                    Budgets::default(),
                    None,
                    solver,
                );
                let error = result.unwrap_err();
                assert_eq!(error.code, ErrorCode::Unsupported, "{solver:?}: {error}");
                assert!(error.message.contains("not representable"));
            }
        }
    }

    #[test]
    fn nonlinear_jacobian_is_refreshed_within_a_pharmflux_interval() {
        for solver in [Solver::DiffsolRosenbrock23, Solver::DiffsolRodas5p] {
            let model = NonlinearDecay {
                jacobian_times: RefCell::new(Vec::new()),
            };
            let rows = simulate_with_solver(
                &model,
                &protocol(0.0, 1.0, vec![1.0]),
                Budgets::default(),
                None,
                solver,
            )
            .unwrap();
            let final_value = rows.last().unwrap().values[0];
            let exact = 2.0 / 3.0;
            assert!(
                (final_value - exact).abs() < 2e-5,
                "{solver:?}: {final_value}"
            );
            assert!(
                model
                    .jacobian_times
                    .borrow()
                    .iter()
                    .any(|time| *time > 0.0 && *time < 1.0),
                "{solver:?} reused its initial nonlinear Jacobian"
            );
        }
    }
}
