use std::{
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "pharmflux-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn write(&self, name: &str, bytes: impl AsRef<[u8]>) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}
fn run(args: &[&std::ffi::OsStr]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_pharmflux"))
        .args(args)
        .output()
        .unwrap()
}
#[test]
fn model_size_budget_is_separate_from_requests_and_text() {
    use pharmflux::document::{CompiledDocument, MODEL_DOCUMENT_BYTE_LIMIT};
    let temp = Temp::new();
    let mut source =
        std::fs::read_to_string(root().join("conformance/models/synthetic-one-compartment.json"))
            .unwrap();
    source.extend(std::iter::repeat_n(
        ' ',
        MODEL_DOCUMENT_BYTE_LIMIT - source.len(),
    ));
    assert!(CompiledDocument::from_json(&source).is_ok());
    let path = temp.write("large.json", &source);
    let output = run(&["validate".as_ref(), path.as_os_str()]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    source.push(' ');
    assert!(CompiledDocument::from_json(&source).is_err());
    let too_large = temp.write("too-large.json", source);
    assert_eq!(
        run(&["validate".as_ref(), too_large.as_os_str()])
            .status
            .code(),
        Some(2)
    );
    let request = temp.write("oversize-request.json", vec![b' '; 1_000_001]);
    let output = run(&["run".as_ref(), path.as_os_str(), request.as_os_str()]);
    let error: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(error["error"]["code"], "input_limit");
}
#[test]
fn packaged_command_matches_library_scientific_identity_and_results() {
    let source = root().join("conformance/models/synthetic-pbpk-24.pfx");
    let request = root().join("conformance/requests/synthetic-pbpk-24.json");
    let output = run(&["run".as_ref(), source.as_os_str(), request.as_os_str()]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    let actual: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let model = pharmflux::language::parse(&std::fs::read_to_string(source).unwrap())
        .unwrap()
        .compile()
        .unwrap();
    let expected = model
        .execute(&serde_json::from_slice(&std::fs::read(request).unwrap()).unwrap())
        .unwrap();
    assert_eq!(actual, serde_json::to_value(expected).unwrap());
}
#[test]
fn format_accepts_json_roundtrips_and_never_modifies_source() {
    let path = root().join("conformance/models/synthetic-one-compartment.json");
    let before = std::fs::read(&path).unwrap();
    let output = run(&["format".as_ref(), path.as_os_str()]);
    assert!(output.status.success());
    let temp = Temp::new();
    let formatted = temp.write("formatted model.pfx", &output.stdout);
    let second = run(&["format".as_ref(), formatted.as_os_str()]);
    assert!(second.status.success());
    assert_eq!(output.stdout, second.stdout);
    let valid = run(&["validate".as_ref(), formatted.as_os_str()]);
    assert!(valid.status.success());
    let report: serde_json::Value = serde_json::from_slice(&valid.stdout).unwrap();
    assert_eq!(report["scope"], "model_compilation");
    assert_eq!(report["valid"], true);
    assert_eq!(before, std::fs::read(path).unwrap());
}
#[test]
fn scientific_and_input_failures_emit_only_structured_stderr() {
    let temp = Temp::new();
    let model = root().join("conformance/models/synthetic-invariant-loss.pfx");
    let mut request: serde_json::Value = serde_json::from_slice(
        &std::fs::read(root().join("conformance/requests/synthetic-invariant-loss.json")).unwrap(),
    )
    .unwrap();
    request["regimen"]["end"]["value"] = 2.into();
    let bad = temp.write("crossing.json", serde_json::to_vec(&request).unwrap());
    let output = run(&["run".as_ref(), model.as_os_str(), bad.as_os_str()]);
    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    let error: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(error["error"]["code"], "invariant");
    assert!(error["error"]["time"].as_f64().unwrap() > 1.);
    let malformed = temp.write("invalid.json", b"{not json}");
    let huge = temp.write("huge.pfx", vec![b'x'; 1_000_001]);
    for args in [
        vec!["wat".as_ref()],
        vec!["run".as_ref(), model.as_os_str(), malformed.as_os_str()],
        vec!["validate".as_ref(), huge.as_os_str()],
        vec!["validate".as_ref(), temp.0.as_os_str()],
    ] {
        let output = run(&args);
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        assert!(serde_json::from_slice::<serde_json::Value>(&output.stderr).is_ok());
    }
    let help = run(&["--help".as_ref()]);
    assert!(help.status.success());
    assert!(help.stderr.is_empty());
}

#[test]
fn steady_state_command_matches_typed_library_and_fails_without_output() {
    let source = root().join("conformance/models/synthetic-linear-pk.json");
    let request = root().join("conformance/requests/synthetic-linear-pk-steady.json");
    let output = run(&[
        "steady-state".as_ref(),
        source.as_os_str(),
        request.as_os_str(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let actual: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let model = pharmflux::document::CompiledDocument::from_json(
        &std::fs::read_to_string(&source).unwrap(),
    )
    .unwrap();
    let q = serde_json::from_slice(&std::fs::read(&request).unwrap()).unwrap();
    let expected = model.periodic_steady_state(&q).unwrap();
    assert_eq!(actual, serde_json::to_value(expected).unwrap());
    let mut invalid: serde_json::Value =
        serde_json::from_slice(&std::fs::read(request).unwrap()).unwrap();
    invalid["run"]["budgets"]["solver_callbacks"] = serde_json::json!(0);
    let tmp = Temp::new();
    let path = tmp.write("invalid.json", invalid.to_string());
    let failure = run(&[
        "steady-state".as_ref(),
        source.as_os_str(),
        path.as_os_str(),
    ]);
    assert_eq!(failure.status.code(), Some(3));
    assert!(failure.stdout.is_empty());
    let error: serde_json::Value = serde_json::from_slice(&failure.stderr).unwrap();
    assert_eq!(error["error"]["code"], "work_budget");
}

#[test]
fn nonlinear_cycle_command_matches_typed_result() {
    let source = root().join("conformance/models/synthetic-saturable-periodic.json");
    let request = root().join("conformance/requests/synthetic-saturable-periodic.json");
    let output = run(&[
        "iterate-periodic".as_ref(),
        source.as_os_str(),
        request.as_os_str(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let compiled =
        pharmflux::document::CompiledDocument::from_json(&std::fs::read_to_string(source).unwrap())
            .unwrap();
    let q = serde_json::from_slice(&std::fs::read(request).unwrap()).unwrap();
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&output.stdout).unwrap(),
        serde_json::to_value(compiled.iterate_periodic_state(&q).unwrap()).unwrap()
    );
}

#[test]
fn fit_command_recovers_individual_and_pooled_parameters_and_reports_limits() {
    use pharmflux_core::fit::{FitProblem, FitRequest, FitStatus, PooledGaussianFitRequest};
    let source = root().join("conformance/models/synthetic-one-compartment.json");
    let mut request: FitRequest = serde_json::from_slice(
        &std::fs::read(root().join("conformance/requests/synthetic-fit.json")).unwrap(),
    )
    .unwrap();
    let document = serde_json::from_slice(&std::fs::read(&source).unwrap()).unwrap();
    let model = pharmflux::document::sensitivity::CompiledSensitivityDocument::compile(
        &document,
        request.parameters(),
    )
    .unwrap();
    let temp = Temp::new();
    for pooled in [false, true] {
        if pooled {
            let FitProblem::Individual(fit) = request.problem.clone() else {
                panic!()
            };
            let mut second = fit.objective.clone();
            second.simulation.run.regimen.administrations[0]
                .amount
                .value *= 2.;
            for observation in &mut second.observations {
                observation.value.value *= 2.;
            }
            request.problem = FitProblem::Pooled(PooledGaussianFitRequest {
                fit,
                additional_subjects: vec![second],
            });
        }
        let path = temp.write("fit.json", serde_json::to_vec(&request).unwrap());
        let output = run(&["fit".as_ref(), source.as_os_str(), path.as_os_str()]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let actual: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        let expected = model.execute_fit(&request).unwrap();
        let fit = expected.fit.as_ref().expect("individual or pooled fit");
        assert_eq!(fit.status, FitStatus::Converged, "{expected:?}");
        assert!((fit.parameters["cl"].value - 0.6).abs() < 1e-5);
        assert!((fit.parameters["v"].value - 4.5).abs() < 1e-5);
        assert_eq!(actual, serde_json::to_value(expected).unwrap());
    }
    let FitProblem::Pooled(ref mut pooled) = request.problem else {
        panic!()
    };
    pooled.fit.max_evaluations = 1;
    let path = temp.write("limit.json", serde_json::to_vec(&request).unwrap());
    let output = run(&["fit".as_ref(), source.as_os_str(), path.as_os_str()]);
    assert!(output.status.success());
    let actual: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(actual["fit"]["status"], "evaluation_limit");
    let mut bad = serde_json::to_value(&request).unwrap();
    bad["schema"] = "pharmflux.fit/v99".into();
    let path = temp.write("bad.json", bad.to_string());
    let output = run(&["fit".as_ref(), source.as_os_str(), path.as_os_str()]);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
}

#[test]
fn scan_command_matches_typed_results_and_returns_no_partial_results() {
    let source = root().join("conformance/models/synthetic-one-compartment.json");
    let request = root().join("conformance/requests/synthetic-scan.json");
    let output = run(&["scan".as_ref(), source.as_os_str(), request.as_os_str()]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let compiled = pharmflux::document::CompiledDocument::from_json(
        &std::fs::read_to_string(&source).unwrap(),
    )
    .unwrap();
    let q = serde_json::from_slice(&std::fs::read(&request).unwrap()).unwrap();
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&output.stdout).unwrap(),
        serde_json::to_value(compiled.scan(&q).unwrap()).unwrap()
    );
    let mut bad: serde_json::Value =
        serde_json::from_slice(&std::fs::read(request).unwrap()).unwrap();
    bad["points"][1]["parameters"]["cl"]["unit"] = "h".into();
    let temp = Temp::new();
    let path = temp.write("invalid.json", bad.to_string());
    let output = run(&["scan".as_ref(), source.as_os_str(), path.as_os_str()]);
    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
}

#[test]
fn morris_command_preserves_generated_and_explicit_design_evidence() {
    use pharmflux_core::morris::*;
    let source = root().join("conformance/models/synthetic-one-compartment.json");
    let path = root().join("conformance/requests/synthetic-morris.json");
    let mut request: MorrisExecutionRequest =
        serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    let compiled = pharmflux::document::CompiledDocument::from_json(
        &std::fs::read_to_string(&source).unwrap(),
    )
    .unwrap();
    let tmp = Temp::new();
    let mut previous_effects = None;
    let mut previous_hash = None;
    for explicit in [false, true] {
        if explicit {
            let MorrisProblem::Generated(ref r) = request.problem else {
                panic!()
            };
            request.problem = MorrisProblem::Explicit(generate_morris_design(r).unwrap().request);
        }
        let path = tmp.write("morris.json", serde_json::to_vec(&request).unwrap());
        let output = run(&["morris".as_ref(), source.as_os_str(), path.as_os_str()]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let expected = compiled.execute_morris(&request).unwrap();
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&output.stdout).unwrap(),
            serde_json::to_value(&expected).unwrap()
        );
        assert_eq!(expected.design.is_some(), !explicit);
        let effects = serde_json::to_value(&expected.result.effects).unwrap();
        if let Some(previous) = previous_effects {
            assert_eq!(previous, effects);
        }
        if let Some(previous) = previous_hash {
            assert_ne!(previous, expected.result.analysis_request_hash);
        }
        previous_effects = Some(effects);
        previous_hash = Some(expected.result.analysis_request_hash);
    }
    let mut invalid = serde_json::to_value(&request).unwrap();
    invalid["schema"] = "pharmflux.morris/v99".into();
    let path = tmp.write("invalid.json", invalid.to_string());
    let output = run(&["morris".as_ref(), source.as_os_str(), path.as_os_str()]);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
}
