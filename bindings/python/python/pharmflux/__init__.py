"""Compile once and run versioned pharmacology requests with the native engine."""
import json
from . import _native

__version__ = _native.__version__
__all__ = ["CompiledSensitivities", "CompiledModel", "PharmfluxError", "PopulationFitDataset", "population_fit_dataset", "format_model", "__version__"]


class PharmfluxError(ValueError):
    """A structured compiler or scientific execution failure with no partial result."""
    def __init__(self, diagnostic):
        self.diagnostic = diagnostic
        self.code = diagnostic.get("code", "internal")
        self.time = diagnostic.get("time")
        self.expression = diagnostic.get("expression")
        self.line = diagnostic.get("line")
        self.column = diagnostic.get("column")
        super().__init__(diagnostic.get("message", self.code))


def _encode(value):
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        raise TypeError("expected model text, JSON text, or a dict")
    try:
        return json.dumps(value, allow_nan=False)
    except (ValueError, TypeError) as error:
        raise PharmfluxError({"code": "invalid_input", "message": str(error)}) from error


def _call(function, *args):
    try:
        return function(*args)
    except _native.NativeError as error:
        raise PharmfluxError(json.loads(str(error))) from None


class CompiledModel:
    """Immutable compiled model; independent request bindings on each run.

    Compilation and native execution release the GIL. A compiled model can be
    shared across Python threads; returned dicts and numeric lists are independent.
    """
    __slots__ = ("_model",)

    def __init__(self, source):
        self._model = _call(_native.CompiledModel, _encode(source))

    @property
    def model_content_hash(self):
        return self._model.model_content_hash

    @property
    def output_names(self):
        return self._model.output_names

    @property
    def output_units(self):
        return self._model.output_units

    def morris_json(self, request):
        return _call(self._model.morris_json, _encode(request))

    def morris(self, request):
        """Execute explicit or seeded Morris screening with bounded scans."""
        return json.loads(self.morris_json(request))

    def scan_json(self, request):
        return _call(self._model.scan_json, _encode(request))

    def scan(self, request):
        """Run up to 25 explicit parameter points under aggregate budgets."""
        return json.loads(self.scan_json(request))

    def run(self, request):
        """Return a pharmflux.result/v0.1 dict, including execution identity."""
        return json.loads(self.run_json(request))

    def periodic_steady_state(self, request):
        """Return periodic pre-event amounts, residuals and execution identity."""
        return json.loads(_call(self._model.periodic_steady_state_json, _encode(request)))

    def iterate_periodic_state(self, request):
        """Iterate a phase-periodic model; return verified state and cycle count."""
        return json.loads(_call(self._model.iterate_periodic_state_json, _encode(request)))

    def run_json(self, request):
        """Return the complete result as JSON for storage or interprocess transfer."""
        return _call(self._model.execute_json, _encode(request))

    def fit_scalar_json(self, request):
        """Fit one bounded parameter through native simulations without sensitivities."""
        return _call(self._model.fit_scalar_json, _encode(request))

    def fit_scalar(self, request):
        return json.loads(self.fit_scalar_json(request))


def format_model(source):
    """Return canonical text after validating current engine support."""
    return _call(_native.format_model, _encode(source))


class CompiledSensitivities:
    """Compile one ordered parameter selection; each run binds independent values.

    Native compilation/execution release the GIL. Instances can be shared across
    threads; the request's with_respect_to must match the compiled selection.
    """
    __slots__ = ("_model",)

    def __init__(self, source, with_respect_to):
        self._model = _call(_native.CompiledSensitivities, _encode(source), with_respect_to)

    def run_json(self, request):
        return _call(self._model.execute_json, _encode(request))

    def run(self, request):
        return json.loads(self.run_json(request))

    def fit_json(self, request):
        """Return a versioned fit result; inspect fit.status for convergence."""
        return _call(self._model.fit_json, _encode(request))

    def fit(self, request):
        """Fit an individual, pooled, or population request in native code."""
        return json.loads(self.fit_json(request))

    def fit_population_data(self, dataset):
        """Fit a prepared FOCEI/SAEM dataset and retain source-row provenance."""
        if not isinstance(dataset, PopulationFitDataset):
            raise TypeError("dataset must come from population_fit_dataset")
        result = self.fit(dataset.request)
        population = result.get("population")
        if not isinstance(population, dict):
            raise PharmfluxError({"code": "internal", "message": "population fit returned no population result"})
        rows = population.get("fitted_observations")
        mapped = None
        if rows is not None:
            by_row = {(item["subject_index"], item["observation_index"]): item
                      for item in dataset.mapping}
            if len(rows) != len(by_row):
                raise PharmfluxError({"code": "internal", "message": "fitted observation count differs from source mapping"})
            mapped = []
            seen = set()
            for row in rows:
                key = (row["subject_index"], row["observation_index"])
                source = by_row.get(key)
                if key in seen or source is None or row["row"] != source["result_row"]:
                    raise PharmfluxError({"code": "internal", "message": "fitted observation row differs from source mapping"})
                seen.add(key)
                mapped.append({**row, "subject_id": source["subject_id"],
                               "source_row": source["source_row"]})
        return {"result": result, "subject_ids": list(dataset.subject_ids),
                "mapping": [item.copy() for item in dataset.mapping],
                "fitted_observations": mapped}


from .dataset import PopulationFitDataset, population_fit_dataset
