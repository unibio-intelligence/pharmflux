# Generated structural schemas

These Draft 7 schemas describe the current serialized Rust records for models,
run requests, results, expressions and regimens. They are generated from the
same serde types used by execution, with the optional `pharmflux-core/schema`
feature and schemars 0.8.22. They are draft contracts for the implemented subset.

Schema validation checks document shape, required fields, enums, expression
argument shapes and rejection of additional properties where the Rust type
rejects them. It does not establish scientific validity or runnable support.
Compilation and binding still validate units, symbol references, domains,
limits and supported capabilities. Execution validates tolerance modes and
rejects non-null seeds. Structurally valid documents can fail these checks.

Generate or verify without a solver dependency:

```sh
cargo run -p pharmflux-core --features schema --example export_schemas -- spec/schema
cargo run -p pharmflux-core --features schema --example export_schemas -- spec/schema --check
```

The generator writes files only in the example executable. Library callers can
obtain schema values from `pharmflux_core::schema::documents()`. The feature is
absent from the default WASM build; schema tooling need not ship to browsers.
