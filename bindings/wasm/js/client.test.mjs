import test from 'node:test';
import assert from 'node:assert/strict';
import { WorkerClient } from './client.mjs';
function harness() {
  const workers = [];
  const client = new WorkerClient(() => {
    const worker = { terminated: false, messages: [], postMessage(x) { this.messages.push(x); }, terminate() { this.terminated = true; } };
    workers.push(worker); return worker;
  });
  return { client, workers };
}
test('queued old-worker messages cannot resolve a newer request', async () => {
  const {client, workers} = harness();
  const old = client.run({edit: 1});
  const rejection = assert.rejects(old, {code:'superseded'});
  const stale = workers[0].onmessage;
  const current = client.run({edit: 2});
  assert.equal(workers[0].terminated, true);
  stale({data:{type:'result',requestId:2,result:'stale'}});
  workers[1].onmessage({data:{type:'result',requestId:1,result:'wrong request'}});
  workers[1].onmessage({data:{type:'result',requestId:2,result:'latest'}});
  assert.deepEqual(await current, {requestId:2,result:'latest'});
  await rejection; client.dispose();
});
test('timeout terminates the worker and the next request can recover', async () => {
  const {client, workers} = harness();
  await assert.rejects(client.run({}, {timeoutMs:5}), {code:'time_budget'});
  assert.equal(workers[0].terminated,true);
  const next = client.run({});
  workers[1].onmessage({data:{type:'result',requestId:2,result:42}});
  assert.equal((await next).result,42);client.dispose();
});
test('dispose rejects active work and prevents reuse', async () => {
  const {client,workers} = harness(); const active=client.run({});
  client.dispose();await assert.rejects(active,{code:'disposed'});
  await assert.rejects(client.run({}),{code:'disposed'});
  assert.equal(workers.length,1);assert.equal(workers[0].terminated,true);
});
test('startup and transport failures release the request', async () => {
  const client=new WorkerClient(()=>{throw new Error('creation failed');});
  await assert.rejects(client.run({}),{code:'worker_failure'});client.dispose();
  const h=harness();const active=h.client.run({});
  h.workers[0].onmessageerror();await assert.rejects(active,{code:'worker_protocol'});
  assert.equal(h.workers[0].terminated,true);h.client.dispose();
});
