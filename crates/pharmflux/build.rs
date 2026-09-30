use sha2::{Digest, Sha256};
use std::{
    env, fs,
    path::{Path, PathBuf},
    process::Command,
};
fn sources(path: &Path, files: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(path).expect("source directory") {
        let path = entry.expect("source entry").path();
        if path.is_dir() {
            sources(&path, files);
        } else if path.extension().is_some_and(|s| s == "rs") {
            files.push(path);
        }
    }
}
fn all_files(path: &Path, files: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(path).expect("vendored dependency directory") {
        let path = entry.expect("vendored dependency entry").path();
        if path.is_dir() {
            all_files(&path, files);
        } else if path.is_file() {
            files.push(path);
        }
    }
}
fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap()).join("../..");
    let mut files = Vec::new();
    for directory in [
        "crates/pharmflux-core/src",
        "crates/pharmflux/src",
        "bindings/wasm/src",
    ] {
        let dir = root.join(directory);
        println!("cargo:rerun-if-changed={}", dir.display());
        sources(&dir, &mut files);
    }
    for file in [
        "Cargo.toml",
        "Cargo.lock",
        ".cargo/config.toml",
        "crates/pharmflux-core/Cargo.toml",
        "crates/pharmflux/Cargo.toml",
        "bindings/wasm/Cargo.toml",
        "crates/pharmflux/build.rs",
        "vendor/DIFFSOL-PROVENANCE.md",
    ] {
        files.push(root.join(file));
    }
    files.sort();
    let mut hash = Sha256::new();
    for path in files {
        println!("cargo:rerun-if-changed={}", path.display());
        let relative = path
            .strip_prefix(&root)
            .unwrap()
            .to_str()
            .unwrap()
            .as_bytes();
        let data = fs::read(&path).expect("source bytes");
        hash.update((relative.len() as u64).to_be_bytes());
        hash.update(relative);
        hash.update((data.len() as u64).to_be_bytes());
        hash.update(data);
    }
    let vendored_root = root.join("vendor/diffsol");
    println!("cargo:rerun-if-changed={}", vendored_root.display());
    let mut vendored_files = Vec::new();
    all_files(&vendored_root, &mut vendored_files);
    vendored_files.sort();
    let mut vendored_hash = Sha256::new();
    for path in vendored_files {
        println!("cargo:rerun-if-changed={}", path.display());
        let relative = path
            .strip_prefix(&root)
            .unwrap()
            .to_str()
            .unwrap()
            .as_bytes();
        let data = fs::read(&path).expect("vendored dependency bytes");
        for digest in [&mut hash, &mut vendored_hash] {
            digest.update((relative.len() as u64).to_be_bytes());
            digest.update(relative);
            digest.update((data.len() as u64).to_be_bytes());
            digest.update(&data);
        }
    }
    println!(
        "cargo:rustc-env=PHARMFLUX_DIFFSOL_SOURCE_SHA256={:x}",
        vendored_hash.finalize()
    );
    println!(
        "cargo:rustc-env=PHARMFLUX_SOURCE_SHA256=sha256:{:x}",
        hash.finalize()
    );
    for key in ["TARGET", "PROFILE"] {
        println!("cargo:rustc-env=PHARMFLUX_{key}={}", env::var(key).unwrap());
    }
    println!("cargo:rerun-if-env-changed=CARGO_ENCODED_RUSTFLAGS");
    println!(
        "cargo:rustc-env=PHARMFLUX_RUSTFLAGS_SHA256=sha256:{:x}",
        Sha256::digest(
            env::var("CARGO_ENCODED_RUSTFLAGS")
                .unwrap_or_default()
                .as_bytes()
        )
    );
    let rustc = Command::new(env::var("RUSTC").unwrap())
        .arg("--version")
        .output()
        .expect("rustc version");
    assert!(rustc.status.success());
    println!(
        "cargo:rustc-env=PHARMFLUX_RUSTC_VERSION={}",
        String::from_utf8(rustc.stdout).unwrap().trim()
    );
}
