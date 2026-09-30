//! Offline file and process adapter. Scientific libraries remain free of file I/O.
use pharmflux::{document::CompiledDocument, language};
use pharmflux_core::{
    fit::FitRequest,
    model::ModelDocument,
    run::{CycleIterationRequest, RunRequest, SensitivityRequest, SteadyStateRequest},
};
use serde_json::{json, Value};
use std::{
    ffi::OsString,
    fs::File,
    io::{Read, Write},
    path::Path,
    process::ExitCode,
    sync::Arc,
};
const INPUT_LIMIT: u64 = 1_000_000;
const HELP: &str = "pharmflux — offline pharmacology modeling\n\nUsage:\n  pharmflux run MODEL REQUEST.json\n  pharmflux morris MODEL REQUEST.json\n  pharmflux scan MODEL REQUEST.json\n  pharmflux fit MODEL REQUEST.json\n  pharmflux sensitivities MODEL REQUEST.json\n  pharmflux steady-state MODEL REQUEST.json\n  pharmflux iterate-periodic MODEL REQUEST.json\n  pharmflux validate MODEL\n  pharmflux format MODEL\n  pharmflux --version\n\nMODEL is a .pfx text or model JSON file (detected by content).\nJSON models are limited to 4,000,000 bytes; text models and requests\nare limited to 1,000,000 bytes. Results go to stdout;\nstructured errors go to stderr. Files are never modified.\n";
struct Failure {
    exit: u8,
    value: Value,
}
impl Failure {
    fn new(exit: u8, code: &str, message: impl ToString) -> Self {
        Self {
            exit,
            value: json!({"schema":"pharmflux.cli-error/v0.1","error":{"code":code,"message":message.to_string()}}),
        }
    }
    fn scientific(error: pharmflux_core::Error) -> Self {
        Self {
            exit: 3,
            value: json!({"schema":"pharmflux.cli-error/v0.1","error":error}),
        }
    }
    fn language(error: language::Diagnostic) -> Self {
        Self {
            exit: 3,
            value: json!({"schema":"pharmflux.cli-error/v0.1","error":{"code":error.code,"message":error.message,"line":error.line,"column":error.column,"hint":error.hint}}),
        }
    }
}
fn read(path: &Path, limit: u64) -> Result<String, Failure> {
    // Reject known pipes/devices before opening, which could otherwise block.
    let metadata = std::fs::metadata(path)
        .map_err(|e| Failure::new(2, "input_io", format!("{}: {e}", path.display())))?;
    if !metadata.is_file() {
        return Err(Failure::new(2, "input_io", "input must be a regular file"));
    }

    let file = File::open(path)
        .map_err(|e| Failure::new(2, "input_io", format!("{}: {e}", path.display())))?;
    let metadata = file
        .metadata()
        .map_err(|e| Failure::new(2, "input_io", e))?;
    if !metadata.is_file() {
        return Err(Failure::new(2, "input_io", "input must be a regular file"));
    }
    if metadata.len() > limit {
        return Err(Failure::new(
            2,
            "input_limit",
            format!("input exceeds {limit} bytes"),
        ));
    }
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| Failure::new(2, "input_io", e))?;
    if bytes.len() as u64 > limit {
        return Err(Failure::new(
            2,
            "input_limit",
            format!("input exceeds {limit} bytes"),
        ));
    }
    String::from_utf8(bytes).map_err(|e| Failure::new(2, "input_encoding", e))
}
fn model(source: &str) -> Result<(ModelDocument, Arc<CompiledDocument>), Failure> {
    if !source.trim_start().starts_with('{') && source.len() as u64 > INPUT_LIMIT {
        return Err(Failure::new(
            2,
            "input_limit",
            "model text exceeds one million bytes",
        ));
    }
    if source.trim_start().starts_with('{') {
        let document: ModelDocument =
            serde_json::from_str(source).map_err(|e| Failure::new(2, "invalid_json", e))?;
        let compiled = CompiledDocument::compile(&document).map_err(Failure::scientific)?;
        Ok((document, compiled))
    } else {
        let parsed = language::parse(source).map_err(Failure::language)?;
        let compiled = parsed.compile().map_err(Failure::language)?;
        Ok((parsed.document, compiled))
    }
}
fn command(args: &[OsString]) -> Result<String, Failure> {
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        return Ok(HELP.into());
    }
    if args.len() == 1 && args[0] == "--version" {
        return Ok(format!("pharmflux {}\n", env!("CARGO_PKG_VERSION")));
    }
    let action = args.first().and_then(|a| a.to_str()).unwrap_or("");
    let expected = match action {
        "run" | "steady-state" | "iterate-periodic" | "sensitivities" | "fit" | "scan"
        | "morris" => 3,
        "validate" | "format" => 2,
        _ => return Err(Failure::new(2, "usage", HELP)),
    };
    if args.len() != expected {
        return Err(Failure::new(2, "usage", HELP));
    }
    let source = read(
        Path::new(&args[1]),
        pharmflux::document::MODEL_DOCUMENT_BYTE_LIMIT as u64,
    )?;
    let (document, compiled) = model(&source)?;
    let result = match action {
        "run" => {
            let request = read(Path::new(&args[2]), INPUT_LIMIT)?;
            let request: RunRequest =
                serde_json::from_str(&request).map_err(|e| Failure::new(2, "invalid_json", e))?;
            let result = compiled.execute(&request).map_err(Failure::scientific)?;
            serde_json::to_value(result).map_err(|e| Failure::new(2, "serialization", e))?
        }
        "morris" => {
            let request = read(Path::new(&args[2]), INPUT_LIMIT)?;
            let request =
                serde_json::from_str(&request).map_err(|e| Failure::new(2, "invalid_json", e))?;
            let result = compiled
                .execute_morris(&request)
                .map_err(Failure::scientific)?;
            serde_json::to_value(result).map_err(|e| Failure::new(2, "serialization", e))?
        }
        "scan" => {
            let request = read(Path::new(&args[2]), INPUT_LIMIT)?;
            let request =
                serde_json::from_str(&request).map_err(|e| Failure::new(2, "invalid_json", e))?;
            let result = compiled.scan(&request).map_err(Failure::scientific)?;
            serde_json::to_value(result).map_err(|e| Failure::new(2, "serialization", e))?
        }
        "fit" => {
            let request = read(Path::new(&args[2]), INPUT_LIMIT)?;
            let request: FitRequest =
                serde_json::from_str(&request).map_err(|e| Failure::new(2, "invalid_json", e))?;
            let model = pharmflux::document::sensitivity::CompiledSensitivityDocument::compile(
                &document,
                request.parameters(),
            )
            .map_err(Failure::scientific)?;
            let result = model.execute_fit(&request).map_err(Failure::scientific)?;
            serde_json::to_value(result).map_err(|e| Failure::new(2, "serialization", e))?
        }
        "sensitivities" => {
            let request = read(Path::new(&args[2]), INPUT_LIMIT)?;
            let request: SensitivityRequest =
                serde_json::from_str(&request).map_err(|e| Failure::new(2, "invalid_json", e))?;
            let model = pharmflux::document::sensitivity::CompiledSensitivityDocument::compile(
                &document,
                &request.with_respect_to,
            )
            .map_err(Failure::scientific)?;
            let result = model.run(&request).map_err(Failure::scientific)?;
            serde_json::to_value(result).map_err(|e| Failure::new(2, "serialization", e))?
        }
        "steady-state" => {
            let request = read(Path::new(&args[2]), INPUT_LIMIT)?;
            let request: SteadyStateRequest =
                serde_json::from_str(&request).map_err(|e| Failure::new(2, "invalid_json", e))?;
            let result = compiled
                .periodic_steady_state(&request)
                .map_err(Failure::scientific)?;
            serde_json::to_value(result).map_err(|e| Failure::new(2, "serialization", e))?
        }
        "iterate-periodic" => {
            let request = read(Path::new(&args[2]), INPUT_LIMIT)?;
            let request: CycleIterationRequest =
                serde_json::from_str(&request).map_err(|e| Failure::new(2, "invalid_json", e))?;
            let result = compiled
                .iterate_periodic_state(&request)
                .map_err(Failure::scientific)?;
            serde_json::to_value(result).map_err(|e| Failure::new(2, "serialization", e))?
        }
        "validate" => {
            json!({"schema":"pharmflux.validation/v0.1","valid":true,"scope":"model_compilation","model_id":document.model_id,"model_content_hash":compiled.model_content_hash(),"states":document.states.len(),"outputs":compiled.output_names(),"output_units":compiled.output_units()})
        }
        "format" => return language::format(&document).map_err(Failure::language),
        _ => unreachable!(),
    };
    serde_json::to_string_pretty(&result)
        .map(|s| s + "\n")
        .map_err(|e| Failure::new(2, "serialization", e))
}
fn main() -> ExitCode {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    let result = command(&args).and_then(|output| {
        std::io::stdout()
            .lock()
            .write_all(output.as_bytes())
            .map_err(|e| Failure::new(2, "output_io", e))
    });
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            // Best effort if stderr itself has been closed by the caller.
            let _ = writeln!(std::io::stderr().lock(), "{}", error.value);
            ExitCode::from(error.exit)
        }
    }
}
