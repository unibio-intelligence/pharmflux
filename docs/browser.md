# Browser integration and privacy

The supplied browser client runs model compilation and simulation in a Web
Worker backed by WebAssembly. Build the assets as described in
[installation](installation.md), then bundle or serve the `bindings/wasm/`
JavaScript module and generated `pkg/engine/` files together. Your build tool
must support module workers and resolving `new URL(..., import.meta.url)`.

```js
import { createEngineClient } from './bindings/wasm/js/index.mjs';

const engine = createEngineClient();
try {
  const { result } = await engine.run({
    kind: 'run',
    source: modelText,
    request: runRequest, // A pharmflux.run/v0.1 request.
  });
  console.log(result.names, result.times, result.columns);
} finally {
  engine.dispose();
}
```

`modelText` may contain a `.pfx` model or a JSON model document. `runRequest`
uses the same versioned format as the command-line program; the
[synthetic request](../conformance/requests/synthetic-pbpk-24.json) is a
complete example. The returned columns are typed arrays aligned to `times`;
`sides` marks each sample's pre- or post-event side. `run()` returns a request
ID and result, and accepts `timeoutMs` and `onStarted` options. Call
`dispose()` when the client is no longer needed. The worker accepts model and
request payloads of at most one million bytes each.

The public worker handles **simulation requests only**. The underlying WASM
binding also has sensitivity and fit entry points, but this worker does not
expose them. An embedding application would need to integrate those bindings
explicitly before offering browser-based fitting.

## Privacy boundary

For the supplied worker, model source and simulation requests are sent from
the page to its local Worker and processed by WASM on the user's device. The
worker does not upload them or the resulting trajectories to a PharmFlux
service. This allows a browser application to keep a user's model and data
local during simulation.

The browser still fetches the JavaScript and WASM assets. The embedding
application controls where those assets come from and whether it uploads
inputs or results, logs telemetry, or stores anything locally. Review those
application behaviors before promising users that their data never leaves
their device. Browser execution also does not replace application-level
access control or clinical validation.
