use super::*;
use pharmflux_core::{identity::canonical_hash, scan::*, Budgets};
impl CompiledDocument {
    /// Evaluate 1..=25 explicit parameter points in request order.
    /// Total budgets are conservatively partitioned before execution.
    pub fn scan(self: &Arc<Self>, request: &ScanRequest) -> Result<ScanResult, Error> {
        let count = request.points.len();
        if count == 0 || count > 25 {
            return Err(invalid("scan requires 1..25 explicit points"));
        }
        let total = request.total_budgets;
        if total.solver_callbacks > u32::MAX as u64
            || total.events > u32::MAX as usize
            || total.output_values > u32::MAX as usize
        {
            return Err(invalid(
                "scan budgets must fit the shared unsigned 32-bit range",
            ));
        }
        let mut ids = BTreeSet::new();
        for point in &request.points {
            if point.id.is_empty() || point.id.len() > 128 || !ids.insert(&point.id) {
                return Err(invalid(
                    "scan point IDs must be unique and contain 1..128 bytes",
                ));
            }
            let mut parameters = request.run.parameters.clone();
            parameters.extend(point.parameters.clone());
            self.bind_with_initial_states(
                &parameters,
                &request.run.regimen.covariates,
                &request.run.initial_states,
            )
            .map_err(|mut e| {
                e.expression = Some(format!(
                    "points.{}.{}",
                    point.id,
                    e.expression.unwrap_or_default()
                ));
                e
            })?;
        }
        let budgets = Budgets {
            solver_callbacks: request
                .run
                .budgets
                .solver_callbacks
                .min(total.solver_callbacks / count as u64),
            events: request.run.budgets.events.min(total.events / count),
            output_values: request
                .run
                .budgets
                .output_values
                .min(total.output_values / count),
        };
        let mut results = Vec::with_capacity(count);
        for point in &request.points {
            let mut run = request.run.clone();
            run.parameters.extend(point.parameters.clone());
            run.budgets = budgets;
            let result = self.execute(&run).map_err(|mut e| {
                e.expression = Some(format!(
                    "points.{}.{}",
                    point.id,
                    e.expression.unwrap_or_default()
                ));
                e
            })?;
            results.push(ScanPointResult {
                id: point.id.clone(),
                parameters: run.parameters,
                result,
            });
        }
        Ok(ScanResult {
            schema: ScanResultSchema::V01,
            scan_request_hash: canonical_hash(
                &serde_json::json!({"algorithm":"explicit_parameter_scan_v1","request":request}),
            )?,
            results,
        })
    }
}
