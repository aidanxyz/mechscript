// MechScript web page: the mechanism comes from the link (everything after "#", percent-encoded text),
// is drawn, and can be edited in the box below the drawing.
import { createMechScript, formatProblems } from '../../src/index.js';

const $ = (s) => document.querySelector(s);
const codeEl = $('#code');
const drawingEl = $('#drawing');
const problemsEl = $('#problems');
const statusEl = $('#status');
const copyErrBtn = $('#copy-errors');

let ms = null;
let last = null;

// Lenient decoding: AI-written links are not always perfectly encoded. "+" stays "+" (it is a SMILES
// charge or the species separator, never a space), malformed % sequences are left as they are.
function decodeHash(hash) {
  let raw = hash.replace(/^#/, '');
  if (raw.startsWith('code=')) raw = raw.slice(5);
  if (!raw) return null;
  try { return decodeURIComponent(raw); } catch { /* fall through */ }
  return raw.replace(/(%[0-9a-fA-F]{2})+/g, (m) => { try { return decodeURIComponent(m); } catch { return m; } });
}
// Like encodeURIComponent, plus ( ) ! ' * ~ — parentheses (SMILES branches) would break Markdown links.
function encodeText(text) {
  return encodeURIComponent(text).replace(/[()!'*~]/g, (c) => '%' + c.charCodeAt(0).toString(16).toUpperCase());
}
const pageUrl = () => location.href.split('#')[0];

function render() {
  if (!ms) return;
  const text = codeEl.value;
  const res = ms.toSVG(text, { maxWidth: Math.max(480, Math.min(1400, drawingEl.clientWidth - 24)) });
  last = res;
  drawingEl.innerHTML = res.svg;
  const svg = drawingEl.querySelector('svg');
  if (svg) { svg.removeAttribute('width'); svg.removeAttribute('height'); svg.style.width = `${Math.min(res.width * 1.4, drawingEl.clientWidth - 8)}px`; svg.style.height = 'auto'; }
  problemsEl.innerHTML = '';
  for (const [kind, list] of [['error', res.errors], ['warning', res.warnings]]) {
    for (const p of list) {
      const li = document.createElement('li');
      li.className = kind;
      li.textContent = `line ${p.line}: ${kind}: ${p.message}`;
      problemsEl.appendChild(li);
    }
  }
  const n = res.errors.length + res.warnings.length;
  copyErrBtn.hidden = n === 0;
  statusEl.textContent = n ? '' : 'No problems found.';
}

let timer = null;
codeEl.addEventListener('input', () => {
  clearTimeout(timer);
  timer = setTimeout(() => {
    history.replaceState(null, '', pageUrl() + '#' + encodeText(codeEl.value));
    render();
  }, 250);
});
window.addEventListener('hashchange', () => {
  const t = decodeHash(location.hash);
  if (t != null && t !== codeEl.value) { codeEl.value = t; $('#intro').hidden = true; render(); }
});
window.addEventListener('resize', () => { clearTimeout(timer); timer = setTimeout(render, 200); });

async function copy(text, btn, done) {
  try { await navigator.clipboard.writeText(text); } catch {
    const ta = document.createElement('textarea'); ta.value = text; document.body.appendChild(ta); ta.select(); document.execCommand('copy'); ta.remove();
  }
  const old = btn.textContent; btn.textContent = done; setTimeout(() => (btn.textContent = old), 1500);
}
$('#copy-link').onclick = (e) => copy(pageUrl() + '#' + encodeText(codeEl.value), e.target, 'Link copied');
copyErrBtn.onclick = (e) => copy(
  `The MechScript page reports these problems:\n${formatProblems(last.doc)}\n\nMy current code:\n\`\`\`mechanism\n${codeEl.value}\n\`\`\`\nPlease fix it and send a new link.`,
  e.target, 'Copied — paste it to your AI assistant');
$('#download').onclick = async () => {
  const { svg, width, height } = ms.toSVG(codeEl.value, { color: '#111', background: '#fff' });
  const img = new Image();
  img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
  await img.decode();
  const c = document.createElement('canvas');
  c.width = Math.ceil(width * 3); c.height = Math.ceil(height * 3);
  c.getContext('2d').drawImage(img, 0, 0, c.width, c.height);
  c.toBlob((b) => { const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = (last?.doc?.title || 'mechanism').replace(/[^\w-]+/g, '_') + '.png'; a.click(); }, 'image/png');
};

const fromLink = decodeHash(location.hash);
codeEl.value = fromLink ?? $('#example').textContent.trim() + '\n';
$('#intro').hidden = fromLink != null;
statusEl.textContent = 'Loading the chemistry engine (first visit takes a few seconds)…';
window.initRDKitModule({ locateFile: (f) => f }).then((RDKit) => {
  ms = createMechScript(RDKit);
  statusEl.textContent = '';
  render();
}).catch(() => { statusEl.textContent = 'Could not load the chemistry engine. Check your connection and reload.'; });
