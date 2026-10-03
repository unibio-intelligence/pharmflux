#!/usr/bin/env python3
"""Prepare a standalone R source package from the public workspace.

Registry dependencies must already exist in CARGO_HOME for the offline build.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def prepare(root, destination):
    destination.mkdir(parents=True, exist_ok=False)
    package = destination / 'pharmflux'
    package.mkdir()
    def copy(source, target):
        if source.is_symlink():
            raise ValueError(f'Symlink is not a source file: {source}')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    for name in ['DESCRIPTION', 'NAMESPACE', 'LICENSE', 'src/entrypoint.c']:
        copy(root / 'bindings/r' / name, package / name)
    for directory, suffix in [('R', '.R'), ('man', '.Rd')]:
        for source in sorted((root / 'bindings/r' / directory).glob('*' + suffix)):
            copy(source, package / directory / source.name)
    engine = package / 'src/engine'
    for name in ['Cargo.toml', 'Cargo.lock', '.cargo/config.toml', 'LICENSE']:
        copy(root / name, engine / name)
    members = ['crates/pharmflux-core', 'crates/pharmflux', 'bindings/wasm',
               'bindings/cli', 'bindings/python', 'bindings/r/src/rust']
    for member in members:
        copy(root / member / 'Cargo.toml', engine / member / 'Cargo.toml')
        for source in sorted((root / member / 'src').rglob('*.rs')):
            copy(source, engine / source.relative_to(root))
        if (root / member / 'build.rs').exists():
            copy(root / member / 'build.rs', engine / member / 'build.rs')
    copy(root / 'crates/pharmflux-core/examples/export_schemas.rs',
         engine / 'crates/pharmflux-core/examples/export_schemas.rs')
    copy(root / 'vendor/DIFFSOL-PROVENANCE.md', engine / 'vendor/DIFFSOL-PROVENANCE.md')
    for source in sorted((root / 'vendor/diffsol').rglob('*')):
        if source.is_file():
            copy(source, engine / source.relative_to(root))
    (package / 'src/Makevars').write_text('''CARGO ?= cargo
CARGO_TARGET_DIR = $(CURDIR)/rust-target
PKG_LIBS = $(CARGO_TARGET_DIR)/release/libpharmflux_r.a
.PHONY: rustlib
all: $(SHLIB)
$(SHLIB): rustlib
rustlib:
\t$(CARGO) build --manifest-path engine/bindings/r/src/rust/Cargo.toml --package pharmflux-r --release --offline --locked --target-dir $(CARGO_TARGET_DIR)
''')
    hashes = {str(path.relative_to(package)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(package.rglob('*')) if path.is_file()}
    (destination / 'source-receipt.json').write_text(json.dumps({
        'classification': 'source package from the supplied public workspace',
        'files': hashes}, indent=2) + '\n')
    return package


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--r', default='R')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    package = prepare(root, args.destination.resolve())
    subprocess.run([args.r, 'CMD', 'build', '--no-build-vignettes', str(package)],
                   cwd=package.parent, check=True)
