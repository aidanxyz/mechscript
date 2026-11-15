// Browser entry: <script src="mechscript.js"></script> replaces every ```mechanism code block with its drawing.
//
// Recognised blocks: <pre><code class="language-mechanism">…</code></pre>, <pre class="mechanism">…</pre>,
// <div class="mechanism">…</div>.
// Options: <script src="mechscript.js" data-rdkit="path/to/RDKit_minimal.js" data-start-on-load="false">
// API: MechScript.initialize(), MechScript.render(text, options) -> Promise<{svg, errors, warnings}>,
//      MechScript.run(root?), MechScript.toPNG(svgString, scale?) -> Promise<Blob>

import { createMechScript, formatProblems } from './index.js';
import { examples } from './examples.js';

const RDKIT_VERSION = '2026.9.1';
const DEFAULT_RDKIT = `https://unpkg.com/@rdkit/rdkit@${RDKIT_VERSION}/dist/RDKit_minimal.js`;

const script = typeof document !== 'undefined' ? document.currentScript : null;
const config = {
  rdkitUrl: script?.dataset.rdkit || (typeof window !== 'undefined' && window.MECHSCRIPT_RDKIT_URL) || DEFAULT_RDKIT,
  startOnLoad: script?.dataset.startOnLoad !== 'false',
  render: {},
};

let ready = null;

function loadScript(src) {
  return new Promise((resolve, reject) => {
    const s = document.createElement('script');
    s.src = src;
    s.onload = resolve;
    s.onerror = () => reject(new Error(`could not load RDKit from ${src}`));
    document.head.appendChild(s);
  });
}

function initialize(opts = {}) {
  Object.assign(config, opts);
  if (ready) return ready;
  ready = (async () => {
    if (typeof window.initRDKitModule !== 'function') await loadScript(config.rdkitUrl);
    const base = config.rdkitUrl.slice(0, config.rdkitUrl.lastIndexOf('/') + 1);
    const RDKit = await window.initRDKitModule({ locateFile: (f) => base + f });
    return createMechScript(RDKit);
  })();
  return ready;
}

async function render(text, options = {}) {
  const ms = await initialize();
  return ms.toSVG(text, { ...config.render, ...options });
}

function blocks(root) {
  return [
    ...root.querySelectorAll('pre > code.language-mechanism, pre.mechanism, div.mechanism'),
  ].filter((el) => !el.dataset.mechscriptDone);
}

async function run(root = document) {
  const els = blocks(root);
  if (!els.length) return;
  await initialize();
  for (const el of els) {
    const text = el.textContent;
    const target = el.tagName === 'CODE' ? el.parentElement : el;
    const res = await render(text);
    const fig = document.createElement('figure');
    fig.className = 'mechscript-figure';
    fig.dataset.mechscriptDone = '1';
    fig.innerHTML = res.svg;
    const svg = fig.querySelector('svg');
    svg.style.maxWidth = '100%';
    svg.style.height = 'auto';
    if (res.errors.length) {
      const pre = document.createElement('pre');
      pre.className = 'mechscript-errors';
      pre.style.cssText = 'color:#b42318;white-space:pre-wrap;font-size:12px;margin:4px 0 0';
      pre.textContent = formatProblems(res.doc);
      fig.appendChild(pre);
    }
    fig.mechscriptSource = text;
    target.replaceWith(fig);
  }
}

// SVG string -> PNG Blob (white background).
async function toPNG(svgString, scale = 2) {
  const doc = new DOMParser().parseFromString(svgString, 'image/svg+xml');
  const svg = doc.documentElement;
  const w = parseFloat(svg.getAttribute('width'));
  const h = parseFloat(svg.getAttribute('height'));
  svg.setAttribute('color', '#111');
  const src = new XMLSerializer().serializeToString(svg).replace(/currentColor/g, '#111');
  const img = new Image();
  img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(src);
  await img.decode();
  const canvas = document.createElement('canvas');
  canvas.width = Math.ceil(w * scale);
  canvas.height = Math.ceil(h * scale);
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
  return new Promise((r) => canvas.toBlob(r, 'image/png'));
}

const api = { initialize, render, run, toPNG, examples, formatProblems, config };
if (typeof window !== 'undefined') {
  window.MechScript = api;
  if (config.startOnLoad) {
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', () => run());
    else run();
  }
}
export default api;
