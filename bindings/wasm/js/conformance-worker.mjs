// M0-only worker. Public workers will expose the validated model-IR API.
import init, { run_tmdd_case, run_failure_probe, CompiledSimulation, CompiledSensitivities } from '../pkg/web/pharmflux_wasm.js';
const ready = init();
// One compiled text model per worker; values and regimens are never cached.
let cachedSource = null;
let cachedModel = null;
let sensitivityKey = null;
let sensitivityModel = null;
self.onmessage = async ({ data }) => {
  if (data?.type !== 'run' || !Number.isSafeInteger(data.requestId) || data.requestId <= 0) return;
  const { requestId, payload } = data;
  let trajectory;
  try {
    if (!payload || typeof payload !== 'object') throw JSON.stringify({code:'invalid_input',message:'Run payload is required'});
    const callbacks = payload.callbacks === undefined ? 1000000 : payload.callbacks;
    if (!Number.isInteger(callbacks) || callbacks < 0 || callbacks > 4294967295) {
      throw JSON.stringify({code:'invalid_input',message:'Callback budget must be an unsigned 32-bit integer'});
    }
    await ready;
    if (payload.kind === 'sensitivities' || payload.kind === 'fit') {
      if (typeof payload.source !== 'string' || !payload.request) throw JSON.stringify({code:'invalid_input',message:'Sensitivity source and request are required'});
      if (payload.callbacks !== undefined || payload.protocol !== undefined || payload.regimen !== undefined || payload.overrides !== undefined) throw JSON.stringify({code:'invalid_input',message:'Sensitivities use only their request envelope'});
      let parameters = payload.request.with_respect_to;
      if (payload.kind === 'fit') {
        const problem = payload.request.problem;
        const fit = problem?.kind === 'individual' ? problem.request : problem?.kind === 'pooled' ? problem.request?.fit : null;
        parameters = fit?.objective?.simulation?.with_respect_to;
        if (!Array.isArray(parameters)) throw JSON.stringify({code:'invalid_input',message:'Fit parameter selection is required'});
      }
      const key = JSON.stringify([payload.source, parameters]);
      if (key !== sensitivityKey) {
        const next = new CompiledSensitivities(payload.source, JSON.stringify(parameters));
        sensitivityModel?.free(); sensitivityModel = next; sensitivityKey = key;
      }
      self.postMessage({type:'started',requestId});
      const result = JSON.parse((payload.kind === 'fit' ? sensitivityModel.fit(JSON.stringify(payload.request)) : sensitivityModel.execute(JSON.stringify(payload.request))));
      self.postMessage({type:'result',requestId,result});
      return;
    }
    let names, units;
    if (payload.kind === 'text' || payload.kind === 'run' || payload.kind === 'steady-state' || payload.kind === 'iterate-periodic' || payload.kind === 'scan' || payload.kind === 'morris') {
      if (typeof payload.source !== 'string') throw JSON.stringify({code:'invalid_input',message:'Text source is required'});
      if (payload.source !== cachedSource) {
        const next = CompiledSimulation.from_text(payload.source);
        cachedModel?.free(); cachedModel = next; cachedSource = payload.source;
      }
      names = cachedModel.output_names(); units = cachedModel.output_units();
      // Acknowledge after compilation, immediately before synchronous VM solving.
      self.postMessage({type:'started',requestId});
      if (payload.kind === 'scan' || payload.kind === 'morris') {
        if (payload.callbacks !== undefined || payload.protocol !== undefined || payload.regimen !== undefined || payload.overrides !== undefined) throw JSON.stringify({code:'invalid_input',message:'Scans use only their request envelope'});
        const result = JSON.parse((payload.kind === 'morris' ? cachedModel.morris(JSON.stringify(payload.request)) : cachedModel.scan(JSON.stringify(payload.request))));
        self.postMessage({type:'result',requestId,result});
        return;
      }
      if (payload.kind === 'steady-state') {
        if (payload.callbacks !== undefined || payload.protocol !== undefined || payload.regimen !== undefined || payload.overrides !== undefined) throw JSON.stringify({code:'invalid_input',message:'Steady state uses only its request envelope'});
        const result = JSON.parse(cachedModel.periodic_steady_state(JSON.stringify(payload.request)));
        self.postMessage({type:'result',requestId,result});
        return;
      }
      if (payload.kind === 'iterate-periodic') {
        if (payload.callbacks !== undefined || payload.protocol !== undefined || payload.regimen !== undefined || payload.overrides !== undefined) throw JSON.stringify({code:'invalid_input',message:'Steady state uses only its request envelope'});
        const result = JSON.parse(cachedModel.iterate_periodic_state(JSON.stringify(payload.request)));
        self.postMessage({type:'result',requestId,result});
        return;
      }
      if (payload.kind === 'run') {
        if (payload.callbacks !== undefined || payload.protocol !== undefined || payload.regimen !== undefined || payload.overrides !== undefined) throw JSON.stringify({code:'invalid_input',message:'Scientific runs use only their request envelope'});
        trajectory = cachedModel.execute(JSON.stringify(payload.request));
      } else {
      if ((payload.regimen === undefined) === (payload.protocol === undefined)) throw JSON.stringify({code:'invalid_input',message:'Provide exactly one protocol or regimen'});
      trajectory = payload.regimen === undefined
        ? cachedModel.simulate(JSON.stringify(payload.overrides ?? {}), JSON.stringify(payload.protocol), callbacks)
        : cachedModel.simulate_regimen(JSON.stringify(payload.overrides ?? {}), JSON.stringify(payload.regimen), callbacks);
      }
    } else if (payload.kind === 'probe') {
      self.postMessage({type:'started',requestId});
      trajectory = run_failure_probe(payload.domainError, callbacks);
    } else if (payload.kind === 'tmdd') {
      self.postMessage({type:'started',requestId});
      trajectory = run_tmdd_case(JSON.stringify(payload.case), callbacks);
    } else {
      throw JSON.stringify({ code: 'invalid_input', message: 'Unknown conformance operation' });
    }
    const metadata=trajectory.metadata_json();
    const result = {
      ...(metadata ? JSON.parse(metadata) : {}),
      ...(names ? {names, units} : {}),
      times: trajectory.times(), sides: trajectory.sides(),
      columns: Array.from({ length: trajectory.state_count() }, (_, i) => trajectory.column(i))
    };
    self.postMessage({ type: 'result', requestId, result },
      [result.times.buffer, result.sides.buffer, ...result.columns.map(column => column.buffer)]);
  } catch (error) {
    let diagnostic;
    try { diagnostic = JSON.parse(String(error)); }
    catch { diagnostic = { code: 'worker_failure', message: String(error) }; }
    self.postMessage({ type: 'error', requestId, error: diagnostic });
  } finally {
    trajectory?.free();
  }
};
