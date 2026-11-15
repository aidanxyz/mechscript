// Renders every mechanism (reference examples, test/ai, test/cases) into review/sheet.html:
// one static file, no server, no browser-side RDKit. Shows each drawing with the arrow optimiser
// off and on, plus crossing counts.  Usage: node tools/sheet.js [--only name]
import { readFileSync, readdirSync, writeFileSync, mkdirSync } from 'fs';
import initRDKit from '@rdkit/rdkit';
import { createMechScript, formatProblems } from '../src/index.js';
import { examples } from '../src/examples.js';

const afterOnly = process.argv.includes('--after');
const only = process.argv.includes('--only') ? process.argv[process.argv.indexOf('--only') + 1] : null;
const ms = createMechScript(await initRDKit());
const items = [...examples.map((e) => ({ name: `example: ${e.name}`, code: e.code }))];
for (const dir of ['test/ai', 'test/cases']) {
  for (const f of readdirSync(dir).filter((x) => x.endsWith('.mech')).sort()) items.push({ name: `${dir.split('/')[1]}: ${f}`, code: readFileSync(`${dir}/${f}`, 'utf8') });
}
const esc = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;');
const fmt = (r) => `arrow×arrow ${r.arrowCrossings} · arrow×bond ${r.bondCrossings} · through label ${r.labelHits}`;
const totals = { off: { a: 0, b: 0, l: 0 }, on: { a: 0, b: 0, l: 0 } };
let html = '';
for (const it of items) {
  if (only && !it.name.includes(only)) continue;
  const off = ms.toSVG(it.code, { color: '#111', optimizeArrows: false });
  const on = ms.toSVG(it.code, { color: '#111' });
  for (const [k, r] of [['off', off.report], ['on', on.report]]) { totals[k].a += r.arrowCrossings; totals[k].b += r.bondCrossings; totals[k].l += r.labelHits; }
  const bad = (r) => r.arrowCrossings + r.bondCrossings + r.labelHits > 0;
  html += `<section><h2>${esc(it.name)}</h2>
  <div class="pair">${afterOnly ? '' : `<figure><figcaption class="${bad(off.report) ? 'bad' : ''}">before (no search): ${fmt(off.report)}</figcaption>${off.svg}</figure>`}
  <figure><figcaption class="${bad(on.report) ? 'bad' : ''}">after: ${fmt(on.report)}</figcaption>${on.svg}</figure></div>
  ${on.errors.length || on.warnings.length ? `<pre class="err">${esc(formatProblems(on.doc))}</pre>` : ''}
  <details><summary>code</summary><pre>${esc(it.code)}</pre></details></section>`;
  console.log(`${it.name.padEnd(42)} before: ${fmt(off.report)}   after: ${fmt(on.report)}`);
}
const sum = `before: ${totals.off.a} arrow×arrow, ${totals.off.b} arrow×bond, ${totals.off.l} through label — after: ${totals.on.a} / ${totals.on.b} / ${totals.on.l}`;
console.log('\nTOTAL ' + sum);
mkdirSync('review', { recursive: true });
writeFileSync('review/sheet.html', `<!doctype html><html><head><meta charset="utf-8"><title>MechScript review sheet</title>
<style>body{font-family:system-ui;margin:20px;background:#fff;color:#111} section{border-bottom:1px solid #ddd;padding:10px 0 16px}
h2{font-size:14px;margin:4px 0} .pair{display:flex;flex-direction:column;gap:6px} figure{margin:0} figcaption{font-size:12px;color:#555}
figcaption.bad{color:#b42318} svg{max-width:100%;height:auto} pre{font-size:11px;background:#f6f6f6;padding:6px} pre.err{color:#b42318}</style></head>
<body><h1>MechScript review sheet</h1><p>${sum}</p>${html}</body></html>`);
console.log('wrote review/sheet.html');
