// Playground: MechScript code on the left, live drawing on the right.
import { createMechScript } from '../src/index.js';
import { examples } from '../src/examples.js';

const $ = (s) => document.querySelector(s);
const codeEl = $('#code');
const drawingEl = $('#drawing');
const problemsEl = $('#problems');
const statusEl = $('#status');
const ZOOM = 1.5;

let ms = null;
let doc = null;

// ---------- persistence ----------
const b64 = {
  enc: (s) => btoa(String.fromCharCode(...new TextEncoder().encode(s))).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, ''),
  dec: (s) => new TextDecoder().decode(Uint8Array.from(atob(s.replace(/-/g, '+').replace(/_/g, '/')), (c) => c.charCodeAt(0))),
};
function initialText() {
  const m = location.hash.match(/code=([^&]+)/);
  if (m) { try { return b64.dec(m[1]); } catch { /* fall through */ } }
  try { const saved = localStorage.getItem('mechscript.code'); if (saved) return saved; } catch { /* storage blocked */ }
  return examples[0].code;
}
function persist() {
  try { localStorage.setItem('mechscript.code', codeEl.value); } catch { /* storage blocked */ }
}

// ---------- rendering ----------
function options(forExport = false) {
  return {
    showMaps: $('#opt-maps').checked,
    showLonePairs: $('#opt-lp').checked ? 'all' : 'arrows',
    ...(forExport ? { color: '#111', background: '#fff' } : {}),
  };
}

// Draw at once with the plain layout; if that drawing has problems, a moment later redraw with the full
// pass (the renderer turns molecules to clear crossings, then tests fixes for anything left). That pass
// renders several variants, which would make typing feel slow if it ran on every keystroke.
let fullTimer = null;
function draw(res, t0) {
  doc = res.doc;
  drawingEl.innerHTML = res.svg;
  const svg = drawingEl.querySelector('svg');
  if (svg) { svg.setAttribute('width', res.width * ZOOM); svg.setAttribute('height', res.height * ZOOM); }
  showProblems(res);
  statusEl.textContent = `${doc.states.length} step${doc.states.length === 1 ? '' : 's'} · ${Math.round(performance.now() - t0)} ms`;
}
function renderNow() {
  if (!ms) return;
  const t0 = performance.now();
  const res = ms.toSVG(codeEl.value, { ...options(), autoOrient: false, suggestFixes: false });
  clearTimeout(fullTimer);
  draw(res, t0);
  if (res.report.problems.length && !res.errors.length) {
    const text = codeEl.value;
    fullTimer = setTimeout(() => {
      if (codeEl.value !== text) return;
      const t1 = performance.now();
      draw(ms.toSVG(text, options()), t1);
    }, 400);
  }
}

function showProblems(res) {
  problemsEl.innerHTML = '';
  const layoutPending = res.report.problems.length && !res.warnings.some((w) => /^arrow /.test(w.message));
  const warnings = layoutPending
    ? [...res.warnings, ...res.report.problems.map((p) => ({ line: p.ar.line, message: `arrow ${p.what.join(' and ')}; improving the layout…` }))]
    : res.warnings;
  for (const [kind, list] of [['error', res.errors], ['warning', warnings]]) {
    for (const p of list) {
      const div = document.createElement('div');
      div.className = kind;
      div.innerHTML = `<b>line ${p.line}</b> · `;
      div.appendChild(document.createTextNode(p.message));
      div.onclick = () => selectLine(p.line);
      problemsEl.appendChild(div);
    }
  }
}

let timer = null;
function scheduleRender() {
  clearTimeout(timer);
  timer = setTimeout(() => { renderNow(); persist(); }, 120);
}

function selectLine(line) {
  const lines = codeEl.value.split('\n');
  let start = 0;
  for (let i = 0; i < line - 1; i++) start += lines[i].length + 1;
  codeEl.focus();
  codeEl.setSelectionRange(start, start + (lines[line - 1] || '').length);
}

// ---------- controls ----------
$('#examples').innerHTML += examples.map((ex, i) => `<option value="${i}">${ex.name}</option>`).join('');
$('#examples').onchange = (e) => {
  if (e.target.value === '') return;
  codeEl.value = examples[Number(e.target.value)].code;
  e.target.value = '';
  renderNow(); persist();
};
$('#opt-maps').onchange = renderNow;
$('#opt-lp').onchange = renderNow;
codeEl.addEventListener('input', scheduleRender);
codeEl.addEventListener('keydown', (e) => {
  if (e.key === 'Tab') { e.preventDefault(); document.execCommand('insertText', false, '  '); }
});

// ---------- export / share ----------
function download(name, blob) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
function fileBase() { return (doc?.title || 'mechanism').replace(/[^\w-]+/g, '_').slice(0, 60); }
$('#export-svg').onclick = () => {
  const { svg } = ms.toSVG(codeEl.value, options(true));
  download(fileBase() + '.svg', new Blob([svg], { type: 'image/svg+xml' }));
};
$('#export-png').onclick = async () => {
  const { svg, width, height } = ms.toSVG(codeEl.value, options(true));
  const img = new Image();
  img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
  await img.decode();
  const scale = 3;
  const c = document.createElement('canvas');
  c.width = Math.ceil(width * scale); c.height = Math.ceil(height * scale);
  c.getContext('2d').drawImage(img, 0, 0, c.width, c.height);
  c.toBlob((b) => download(fileBase() + '.png', b), 'image/png');
};
$('#share').onclick = async () => {
  const url = location.href.split('#')[0] + '#code=' + b64.enc(codeEl.value);
  history.replaceState(null, '', url);
  try { await navigator.clipboard.writeText(url); statusEl.textContent = 'Link copied'; } catch { statusEl.textContent = 'Link is in the address bar'; }
};
$('#copy-prompt').onclick = async () => {
  try {
    const spec = await (await fetch('../spec/MECHSCRIPT.md')).text();
    await navigator.clipboard.writeText(spec + '\n\nWrite the mechanism I describe next in MechScript, in a ```mechanism code block.');
    statusEl.textContent = 'Spec copied — paste it into Claude';
  } catch { statusEl.textContent = 'Could not copy'; }
};

// ---------- start ----------
codeEl.value = initialText();
statusEl.textContent = 'Loading RDKit…';
window.initRDKitModule({ locateFile: (f) => '../dist/' + f }).then((RDKit) => {
  ms = createMechScript(RDKit);
  renderNow();
});
