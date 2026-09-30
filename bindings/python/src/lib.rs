//! Python ownership and error boundary; no scientific implementation lives here.
use pharmflux::document::CompiledDocument;
use pyo3::{exceptions::PyException, prelude::*};
use serde_json::{json, Value};
use std::sync::Arc;
pyo3::create_exception!(_native, NativeError, PyException);

fn scientific(error: pharmflux_core::Error) -> Value {
    json!(error)
}
fn diagnostic(error: pharmflux::language::Diagnostic) -> Value {
    json!({"code":error.code,"message":error.message,"line":error.line,"column":error.column,"hint":error.hint})
}
fn invalid(message: impl ToString) -> Value {
    json!({"code":"invalid_input","message":message.to_string()})
}
fn detached<T: Send>(py: Python<'_>, f: impl FnOnce() -> Result<T, Value> + Send) -> PyResult<T> {
    let result = py.detach(|| std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)));
    match result {
        Ok(Ok(value)) => Ok(value),
        Ok(Err(error)) => Err(NativeError::new_err(error.to_string())),
        Err(_) => Err(NativeError::new_err(
            json!({"code":"internal","message":"Rust execution panicked; no result returned"})
                .to_string(),
        )),
    }
}
#[pyclass(frozen, name = "CompiledModel", module = "pharmflux._native")]
struct PythonModel {
    inner: Arc<CompiledDocument>,
}
#[pymethods]
impl PythonModel {
    #[new]
    fn new(py: Python<'_>, source: String) -> PyResult<Self> {
        detached(py, move || {
            let inner = if source.trim_start().starts_with('{') {
                CompiledDocument::from_json(&source).map_err(scientific)?
            } else {
                pharmflux::language::parse(&source)
                    .map_err(diagnostic)?
                    .compile()
                    .map_err(diagnostic)?
            };
            Ok(Self { inner })
        })
    }
    fn morris_json(&self, py: Python<'_>, request: String) -> PyResult<String> {
        let inner = self.inner.clone();
        detached(py, move || {
            if request.len() > 1_000_000 {
                return Err(invalid("request exceeds one million bytes"));
            }
            let request = serde_json::from_str(&request).map_err(invalid)?;
            let result = inner.execute_morris(&request).map_err(scientific)?;
            serde_json::to_string(&result).map_err(invalid)
        })
    }
    fn scan_json(&self, py: Python<'_>, request: String) -> PyResult<String> {
        let inner = self.inner.clone();
        detached(py, move || {
            if request.len() > 1_000_000 {
                return Err(invalid("request exceeds one million bytes"));
            }
            let request = serde_json::from_str(&request).map_err(invalid)?;
            let result = inner.scan(&request).map_err(scientific)?;
            serde_json::to_string(&result).map_err(invalid)
        })
    }
    fn execute_json(&self, py: Python<'_>, request: String) -> PyResult<String> {
        let inner = self.inner.clone();
        detached(py, move || {
            if request.len() > 1_000_000 {
                return Err(invalid("request exceeds one million bytes"));
            }
            let request = serde_json::from_str(&request).map_err(invalid)?;
            let result = inner.execute(&request).map_err(scientific)?;
            serde_json::to_string(&result).map_err(invalid)
        })
    }
    fn fit_scalar_json(&self, py: Python<'_>, request: String) -> PyResult<String> {
        let inner = self.inner.clone();
        detached(py, move || {
            if request.len() > 1_000_000 {
                return Err(invalid("request exceeds one million bytes"));
            }
            let request = serde_json::from_str(&request).map_err(invalid)?;
            let result = inner.fit_scalar(&request).map_err(scientific)?;
            serde_json::to_string(&result).map_err(invalid)
        })
    }
    fn periodic_steady_state_json(&self, py: Python<'_>, request: String) -> PyResult<String> {
        let inner = self.inner.clone();
        detached(py, move || {
            if request.len() > 1_000_000 {
                return Err(invalid("request exceeds one million bytes"));
            }
            let request = serde_json::from_str(&request).map_err(invalid)?;
            let result = inner.periodic_steady_state(&request).map_err(scientific)?;
            serde_json::to_string(&result).map_err(invalid)
        })
    }
    fn iterate_periodic_state_json(&self, py: Python<'_>, request: String) -> PyResult<String> {
        let inner = self.inner.clone();
        detached(py, move || {
            if request.len() > 1_000_000 {
                return Err(invalid("request exceeds one million bytes"));
            }
            let request = serde_json::from_str(&request).map_err(invalid)?;
            let result = inner.iterate_periodic_state(&request).map_err(scientific)?;
            serde_json::to_string(&result).map_err(invalid)
        })
    }
    #[getter]
    fn model_content_hash(&self) -> String {
        self.inner.model_content_hash().into()
    }
    #[getter]
    fn output_names(&self) -> Vec<String> {
        self.inner.output_names().to_vec()
    }
    #[getter]
    fn output_units(&self) -> Vec<String> {
        self.inner.output_units().to_vec()
    }
}
#[pyclass(frozen, name = "CompiledSensitivities", module = "pharmflux._native")]
struct PythonSensitivities {
    inner: Arc<pharmflux::document::sensitivity::CompiledSensitivityDocument>,
}
#[pymethods]
impl PythonSensitivities {
    #[new]
    fn new(py: Python<'_>, source: String, parameters: Vec<String>) -> PyResult<Self> {
        detached(py, move || {
            if source.len() > pharmflux::document::MODEL_DOCUMENT_BYTE_LIMIT {
                return Err(invalid("model exceeds four million bytes"));
            }
            let document = if source.trim_start().starts_with('{') {
                serde_json::from_str(&source).map_err(invalid)?
            } else {
                pharmflux::language::parse(&source)
                    .map_err(diagnostic)?
                    .document
            };
            let inner = pharmflux::document::sensitivity::CompiledSensitivityDocument::compile(
                &document,
                &parameters,
            )
            .map_err(scientific)?;
            Ok(Self { inner })
        })
    }
    fn fit_json(&self, py: Python<'_>, request: String) -> PyResult<String> {
        let inner = self.inner.clone();
        detached(py, move || {
            if request.len() > 1_000_000 {
                return Err(invalid("request exceeds one million bytes"));
            }
            let request = serde_json::from_str(&request).map_err(invalid)?;
            let result = inner.execute_fit(&request).map_err(scientific)?;
            serde_json::to_string(&result).map_err(invalid)
        })
    }
    fn execute_json(&self, py: Python<'_>, request: String) -> PyResult<String> {
        let inner = self.inner.clone();
        detached(py, move || {
            if request.len() > 1_000_000 {
                return Err(invalid("request exceeds one million bytes"));
            }
            let request = serde_json::from_str(&request).map_err(invalid)?;
            let result = inner.run(&request).map_err(scientific)?;
            serde_json::to_string(&result).map_err(invalid)
        })
    }
}
#[pyfunction]
fn format_model(py: Python<'_>, source: String) -> PyResult<String> {
    detached(py, move || {
        if source.len() > pharmflux::document::MODEL_DOCUMENT_BYTE_LIMIT {
            return Err(invalid("model exceeds four million bytes"));
        }
        let document = if source.trim_start().starts_with('{') {
            serde_json::from_str(&source).map_err(invalid)?
        } else {
            pharmflux::language::parse(&source)
                .map_err(diagnostic)?
                .document
        };
        pharmflux::language::format(&document).map_err(diagnostic)
    })
}
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PythonModel>()?;
    m.add_class::<PythonSensitivities>()?;
    m.add_function(wrap_pyfunction!(format_model, m)?)?;
    m.add("NativeError", m.py().get_type::<NativeError>())?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
