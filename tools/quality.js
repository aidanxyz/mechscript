// Drawing-quality harness. Renders every mechanism (src/examples.js, test/ai, test/cases, test/quality),
// measures how good each drawing looks, and writes review/quality.html (contact sheet, worst first, issues
// circled on the drawing) and review/quality.json.
//
//   node tools/quality.js                   measure, compare with test/quality/baseline.json
//   node tools/quality.js --save-baseline   measure and store the result as the new baseline
//   node tools/quality.js --check           exit 1 if any mechanism got worse than the baseline
//   node tools/quality.js --only claisen    only mechanisms whose name contains "claisen"
//
// Metrics (L = bond length in px). Each has a threshold; a value past it counts as an issue.
//   ring      worst interior-angle error of a ring vs the regular polygon (°)      > 8°
//   bondLen   worst bond length error (fraction of L)                             > 0.12
//   angle     worst angle error at a non-ring atom or around a ring substituent (°) > 15°
//   clash     closest pair of non-bonded atoms in one molecule (× L)               < 0.6
//   gap       closest pair of atoms in different molecules of one step (× L)       < 0.9
//   drift     molecule shape change between steps for atoms that kept their bonds (× L, RMS) > 0.25
//   turn      rotation of a molecule between steps (°), or mirrored                 > 25°
//   crossings arrow×arrow, arrow×bond, arrow through a label (from the renderer's own report)  > 0
//   detour    arrow length / distance between its ends                            > 2.6
//   miss      gap between an arrowhead and the atom or bond it points at (× L)     > 0.75
//   tangle    closest approach of two arrows away from their ends (px)            < 4
//   clearance closest approach of an arrow to anything it is not about: another atom's label, a bond not
//             at its own atoms, a charge, a "+" sign, a reaction arrow (px)          < 6
//   loop      an arrow coming back close to itself (px between far-apart parts of one curve) < 6
//   snake     an arrow that bends one way, then the other (S shape)               > 0
import { readFileSync, readdirSync, writeFileSync, mkdirSync, existsSync } from 'fs';
import initRDKit from '@rdkit/rdkit';
import { createMechScript } from '../src/index.js';
import { examples } from '../src/examples.js';

const args = process.argv.slice(2);
const only = args.includes('--only') ? args[args.indexOf('--only') + 1] : null;
const BASELINE = 'test/quality/baseline.json';

export const LIMITS = {
  ring: { max: 8, unit: '°', label: 'ring shape' },
  bondLen: { max: 0.12, unit: 'L', label: 'bond length' },
  angle: { max: 15, unit: '°', label: 'bond angle' },
  clash: { min: 0.6, unit: 'L', label: 'atoms too close' },
  gap: { min: 0.9, unit: 'L', label: 'molecules too close' },
  drift: { max: 0.25, unit: 'L', label: 'shape drift between steps' },
  turn: { max: 25, unit: '°', label: 'turn between steps' },
  crossings: { max: 0, unit: '', label: 'arrow crossings' },
  detour: { max: 2.6, unit: '×', label: 'arrow detour' },
  miss: { max: 0.75, unit: 'L', label: 'arrowhead misses target' },
  tangle: { min: 4, unit: 'px', label: 'arrows too close' },
  clearance: { min: 6, unit: 'px', label: 'arrow skims something' },
  loop: { min: 6, unit: 'px', label: 'arrow touches itself' },
  snake: { max: 0, unit: '', label: 'S-shaped arrow' },
};
const bad = (k, v) => v != null && (LIMITS[k].max != null ? v > LIMITS[k].max + 1e-9 : v < LIMITS[k].min - 1e-9);

// ---------- geometry helpers ----------
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const len = (a) => Math.hypot(a.x, a.y);
const ang = (a) => Math.atan2(a.y, a.x);
const deg = (r) => (r * 180) / Math.PI;
const distSeg = (p, a, b) => {
  const ab = sub(b, a); const l2 = ab.x * ab.x + ab.y * ab.y || 1;
  const t = Math.max(0, Math.min(1, ((p.x - a.x) * ab.x + (p.y - a.y) * ab.y) / l2));
  return len(sub(p, { x: a.x + ab.x * t, y: a.y + ab.y * t }));
};
const boxDist = (p, c, hw, hh) => Math.max(0, Math.max(Math.abs(p.x - c.x) - hw, Math.abs(p.y - c.y) - hh));
// sign changes of the turning direction along a polyline, ignoring near-straight stretches
function inflections(pts) {
  let last = 0; let n = 0;
  for (let i = 2; i < pts.length; i++) {
    const a = sub(pts[i - 1], pts[i - 2]); const b = sub(pts[i], pts[i - 1]);
    const t = Math.atan2(a.x * b.y - a.y * b.x, a.x * b.x + a.y * b.y);
    if (Math.abs(t) < 0.035) continue; // ~2°
    const sg = Math.sign(t);
    if (last && sg !== last) n++;
    last = sg;
  }
  return n;
}
// closest approach between parts of one curve that are far apart along it (a third of the curve or more)
function selfGap(pts) {
  let m = Infinity; let at = null;
  const gap = Math.ceil(pts.length / 3);
  for (let i = 0; i < pts.length; i++) for (let j = i + gap; j < pts.length; j++) {
    // skip pairs that are close only because the whole arrow is short (the ends of a small arch)
    const d = len(sub(pts[i], pts[j]));
    if (d < m) { m = d; at = pts[i]; }
  }
  return { d: m, at };
}
const polyDist = (A, B) => {
  let m = Infinity;
  for (let i = 1; i < A.length; i++) for (const q of B) m = Math.min(m, distSeg(q, A[i - 1], A[i]));
  for (let i = 1; i < B.length; i++) for (const q of A) m = Math.min(m, distSeg(q, B[i - 1], B[i]));
  return m;
};

// ring atoms in walking order, or null when the ring is not a simple cycle in the bond graph
function orderRing(ring, nb) {
  const set = new Set(ring);
  const out = [ring[0]];
  let prev = -1;
  while (out.length < ring.length) {
    const cur = out[out.length - 1];
    const next = nb[cur].filter((n) => set.has(n) && n !== prev && !out.includes(n));
    if (!next.length) return null;
    prev = cur; out.push(next[0]);
  }
  return nb[out[out.length - 1]].includes(out[0]) ? out : null;
}

// best 2D rigid fit of B onto A (rotation, optionally mirrored). Returns { rms, turn (°), mirrored }
function rigidFit(A, B) {
  const ca = A.reduce((s, p) => ({ x: s.x + p.x / A.length, y: s.y + p.y / A.length }), { x: 0, y: 0 });
  const cb = B.reduce((s, p) => ({ x: s.x + p.x / B.length, y: s.y + p.y / B.length }), { x: 0, y: 0 });
  const a = A.map((p) => sub(p, ca));
  const fits = [];
  for (const mirror of [false, true]) {
    const b = B.map((p) => { const q = sub(p, cb); return mirror ? { x: -q.x, y: q.y } : q; });
    let sc = 0; let ss = 0;
    for (let i = 0; i < a.length; i++) { sc += b[i].x * a[i].x + b[i].y * a[i].y; ss += b[i].x * a[i].y - b[i].y * a[i].x; }
    const th = Math.atan2(ss, sc);
    const c = Math.cos(th); const s = Math.sin(th);
    let e = 0;
    for (let i = 0; i < a.length; i++) { const x = c * b[i].x - s * b[i].y; const y = s * b[i].x + c * b[i].y; e += (x - a[i].x) ** 2 + (y - a[i].y) ** 2; }
    const rms = Math.sqrt(e / a.length);
    fits.push({ rms, turn: Math.abs(deg(th)), mirrored: mirror });
  }
  // call it mirrored only when the mirror image fits clearly better (near-symmetric atom sets fit both ways)
  const [plain, mirr] = fits;
  return mirr.rms < plain.rms - 0.2 * 34 ? mirr : plain;
}

// ---------- measuring one rendered mechanism ----------
export function measure(res, L) {
  const doc = res.doc;
  const m = { ring: 0, bondLen: 0, angle: 0, clash: null, gap: null, drift: 0, turn: 0, crossings: 0, detour: 0, miss: 0, tangle: null, clearance: null, loop: null, snake: 0 };
  const units = res.units || [];
  const issues = []; // { metric, value, at: {x,y} | null, path?: [{x,y}], note }
  const worst = (k, v, at, note, path) => {
    if (v == null) return;
    const lim = LIMITS[k];
    if (lim.max != null) m[k] = Math.max(m[k] ?? 0, v); else m[k] = m[k] == null ? v : Math.min(m[k], v);
    if (bad(k, v)) issues.push({ metric: k, value: v, at, path, note });
  };

  doc.states.forEach((st, si) => {
    const placed = []; // atoms of every molecule in this step, global coords
    st.species.forEach((sp, k) => {
      if (!sp.mol || !sp.geo || !sp.geo.P) return;
      const { atoms, bonds } = sp.mol;
      const P = sp.geo.P;
      const G = (i) => ({ x: P[i].x + sp.geo.dx, y: P[i].y + sp.geo.dy });
      const nb = atoms.map((a) => a.bonds.map((bi) => (bonds[bi].a === a.index ? bonds[bi].b : bonds[bi].a)));
      const tag = `step ${si + 1}, molecule ${k + 1}`;

      for (const b of bonds) worst('bondLen', Math.abs(len(sub(P[b.a], P[b.b])) / L - 1), G(b.a), `${tag}: bond ${atoms[b.a].symbol}–${atoms[b.b].symbol}`);

      const ringOf = atoms.map(() => []);
      for (const r of sp.mol.rings || []) {
        if (r.length < 3 || r.length > 8) continue;
        const ord = orderRing(r, nb);
        if (!ord) continue;
        const n = ord.length;
        const ideal = 180 - 360 / n;
        let dev = 0;
        for (let i = 0; i < n; i++) {
          const p = P[ord[(i + n - 1) % n]]; const c = P[ord[i]]; const q = P[ord[(i + 1) % n]];
          let t = Math.abs(deg(ang(sub(p, c)) - ang(sub(q, c)))); if (t > 180) t = 360 - t;
          dev = Math.max(dev, Math.abs(t - ideal));
        }
        // fused / bridged rings can't all be regular; still measured, the note says which kind it is
        const shared = ord.filter((i) => (sp.mol.rings || []).filter((rr) => rr.includes(i)).length > 1).length;
        const kind = shared > 2 ? ' (bridged)' : shared === 2 ? ' (fused)' : '';
        const cen = ord.reduce((s, i) => ({ x: s.x + G(i).x / n, y: s.y + G(i).y / n }), { x: 0, y: 0 });
        worst('ring', dev, cen, `${tag}: ${n}-membered ring${kind}`);
        for (const i of ord) ringOf[i].push(ord);
      }

      // a bond this molecule forms within itself in the next step (a ring closure) gets its own slot: the
      // gap left for it is intended, so count it as a neighbour when judging angles
      const idxOfMap = new Map(atoms.map((a, i) => [a.map, i]).filter(([mp]) => mp != null));
      // (only at atoms that keep all their bonds: in a transfer, say of an H, the new bond takes an old one's place)
      const forming = atoms.map(() => []);
      const nextNb = new Map(); // map -> Set of bonded maps in the next step
      for (const nsp of doc.states[si + 1]?.species || []) {
        for (const b of nsp.mol?.bonds || []) {
          const p = nsp.mol.atoms[b.a].map; const q = nsp.mol.atoms[b.b].map;
          if (p == null || q == null) continue;
          if (!nextNb.has(p)) nextNb.set(p, new Set()); if (!nextNb.has(q)) nextNb.set(q, new Set());
          nextNb.get(p).add(q); nextNb.get(q).add(p);
        }
      }
      const keepsBonds = (i) => nb[i].every((j) => atoms[j].map == null || nextNb.get(atoms[i].map)?.has(atoms[j].map));
      for (const [p, set] of nextNb) {
        const i = idxOfMap.get(p);
        if (i == null || atoms[i].symbol === 'H' || !keepsBonds(i)) continue;
        for (const q of set) {
          const j = idxOfMap.get(q);
          if (j == null || nb[i].includes(j)) continue;
          // the renderer only leaves room for rings of 5 to 7 atoms (see reactionLayout in src/chem.js)
          const dist = new Map([[i, 0]]); const queue = [i];
          while (queue.length && !dist.has(j)) { const x = queue.shift(); for (const y of nb[x]) if (!dist.has(y)) { dist.set(y, dist.get(x) + 1); queue.push(y); } }
          const size = (dist.get(j) ?? Infinity) + 1;
          if (size >= 5 && size <= 7) forming[i].push(j);
        }
      }

      atoms.forEach((a, i) => {
        const nbs = [...nb[i], ...forming[i]];
        if (nbs.length < 2) return;
        const dirs = nbs.map((j) => ang(sub(P[j], P[i]))).sort((x, y) => x - y);
        const gaps = dirs.map((d, j) => deg((j + 1 < dirs.length ? dirs[j + 1] : dirs[0] + 2 * Math.PI) - d));
        if (!ringOf[i].length) {
          // linear for sp atoms (a triple bond, or two double bonds), 120° chains otherwise, even spread at branches
          const orders = a.bonds.map((bi) => bonds[bi].order);
          const linear = nbs.length === 2 && (orders.includes(3) || orders.filter((o) => o === 2).length === 2);
          // four bonds: a cross (90° each), "Y plus one" (60, 60, 120, 120) or a zigzag chain with two
          // substituents (60, 90, 90, 120) are all conventional
          const patterns = nbs.length === 4 ? [[90, 90, 90, 90], [60, 60, 120, 120], [60, 90, 90, 120]] : [gaps.map(() => 360 / nbs.length)];
          const sg = [...gaps].sort((x, y) => x - y);
          const dev = nbs.length === 2 ? Math.abs(Math.min(...gaps) - (linear ? 180 : 120)) : Math.min(...patterns.map((pt) => Math.max(...sg.map((g, j) => Math.abs(g - pt[j])))));
          worst('angle', dev, G(i), `${tag}: ${a.symbol}${a.map != null ? ` (atom ${a.map})` : ''}`);
        } else if (ringOf[i].length === 1 && nbs.length === 3) {
          // one substituent on a ring atom should bisect the outside angle
          const ring = ringOf[i][0];
          const exo = nbs.find((j) => !ring.includes(j));
          if (exo == null) return;
          const ex = ang(sub(P[exo], P[i]));
          const others = nbs.filter((j) => j !== exo).map((j) => ang(sub(P[j], P[i])));
          const d = others.map((o) => { let t = Math.abs(deg(o - ex)); return t > 180 ? 360 - t : t; });
          worst('angle', Math.abs(d[0] - d[1]) / 2, G(i), `${tag}: substituent on ring atom ${a.symbol}${a.map != null ? ` (atom ${a.map})` : ''}`);
        }
      });

      for (let i = 0; i < atoms.length; i++) for (let j = i + 1; j < atoms.length; j++) {
        if (nb[i].includes(j)) continue;
        const d = len(sub(P[i], P[j])) / L;
        worst('clash', d, G(i), `${tag}: ${atoms[i].symbol}${i} and ${atoms[j].symbol}${j}`);
      }
      placed.push({ k, pts: atoms.map((_, i) => G(i)) });
    });
    for (let x = 0; x < placed.length; x++) for (let y = x + 1; y < placed.length; y++) {
      let d = Infinity; let at = null;
      for (const p of placed[x].pts) for (const q of placed[y].pts) { const e = len(sub(p, q)); if (e < d) { d = e; at = p; } }
      worst('gap', d / L, at, `step ${si + 1}: molecules ${placed[x].k + 1} and ${placed[y].k + 1}`);
    }

    // step-to-step: molecules should keep their shape and orientation. Compare atoms that keep the same
    // bonded neighbours (by map number) in a molecule of the previous step.
    if (si > 0) {
      const prev = doc.states[si - 1];
      st.species.forEach((sp, k) => {
        if (!sp.mol || !sp.geo?.P) return;
        const mapsOf = (mol, i) => mol.atoms[i].bonds.map((bi) => { const b = mol.bonds[bi]; return mol.atoms[b.a === i ? b.b : b.a].map; });
        let bestPrev = null; let bestShared = [];
        for (const pp of prev.species) {
          if (!pp.mol || !pp.geo?.P) continue;
          const shared = [];
          sp.mol.atoms.forEach((a, i) => {
            if (a.map == null) return;
            const j = pp.mol.atoms.findIndex((b) => b.map === a.map);
            if (j >= 0) shared.push([i, j]);
          });
          if (shared.length > bestShared.length) { bestPrev = pp; bestShared = shared; }
        }
        if (!bestPrev || bestShared.length < 3) return;
        // keep atoms whose mapped neighbours did not change
        const stable = bestShared.filter(([i, j]) => {
          const A = mapsOf(sp.mol, i).filter((x) => x != null).sort().join();
          const B = mapsOf(bestPrev.mol, j).filter((x) => x != null).sort().join();
          return A === B;
        });
        const use = stable.length >= 3 ? stable : bestShared;
        const fit = rigidFit(use.map(([, j]) => bestPrev.geo.P[j]), use.map(([i]) => sp.geo.P[i]));
        const c = sp.geo.centroid();
        const tag = `step ${si} → ${si + 1}, molecule ${k + 1}`;
        if (stable.length >= 3) worst('drift', fit.rms / L, c, `${tag}: ${stable.length} unchanged atoms`);
        worst('turn', fit.mirrored ? 180 : fit.turn, c, `${tag}${fit.mirrored ? ': mirrored' : ''}`);
      });
    }

    // arrows
    const withGeo = st.arrows.filter((a) => a.geo);
    for (const ar of withGeo) {
      const g = ar.geo;
      const chord = len(sub(g.E, g.S));
      worst('detour', g.length / Math.max(chord, L * 0.5), g.samples[Math.floor(g.samples.length / 2)], `line ${ar.line}: ${ar.text?.trim() ?? ''}`, g.samples);
      // target: the atom (or bond midpoint) the arrow is about; head should stop just outside its label
      worst('miss', len(sub(g.E, g.T)) / L, g.E, `line ${ar.line}: ${ar.text?.trim() ?? ''}`, g.samples);
    }
    // what each arrow passes close to (inner part only: away from where it leaves and arrives)
    const fsz = 15;
    const things = []; // { dist(p) -> px, what }
    for (const sp of st.species) {
      if (!sp.mol || !sp.geo?.obstacles) continue;
      const ob = sp.geo.obstacles();
      const off = (p) => ({ x: p.x + sp.geo.dx, y: p.y + sp.geo.dy });
      const mapOf = (i) => sp.mol.atoms[i].map;
      for (const b of ob.bonds) { const a = off(b.a); const c = off(b.b); things.push({ maps: [mapOf(b.atoms[0]), mapOf(b.atoms[1])], d: (p) => distSeg(p, a, c), what: `${sp.mol.atoms[b.atoms[0]].symbol}–${sp.mol.atoms[b.atoms[1]].symbol} bond` }); }
      for (const l of ob.labels) { const c = off(l.c); things.push({ maps: [mapOf(l.atom)], d: (p) => boxDist(p, c, l.hw - 1.5, l.hh - 1.5), what: `${sp.mol.atoms[l.atom].symbol} label` }); }
      for (const ch of ob.charges) { const c = off(ch.c); things.push({ maps: [mapOf(ch.atom)], d: (p) => Math.max(0, len(sub(p, c)) - ch.r + 1), what: 'charge or dots' }); }
    }
    const row = units.filter((u) => u.kind === 'species' && u.si === si);
    const rowMid = row[0]?.mid;
    for (const u of units) {
      if (rowMid == null || Math.abs(u.mid - rowMid) > 1) continue;
      if (u.kind === 'plus') { const c = { x: u.x + u.w / 2, y: u.mid - fsz * 0.05 }; things.push({ maps: [], d: (p) => boxDist(p, c, fsz * 0.33, fsz * 0.33), what: '"+" sign' }); }
      if (u.kind === 'transition') { const a = { x: u.x, y: u.mid }; const b = { x: u.x + u.w, y: u.mid }; things.push({ maps: [], d: (p) => distSeg(p, a, b), what: 'reaction arrow' }); }
    }
    for (const ar of withGeo) {
      const g = ar.geo;
      const own = new Set([...ar.from.atoms, ...ar.to.atoms]);
      let best = Infinity; let at = null; let what = '';
      for (const t of things) {
        if (t.maps.some((mp) => mp != null && own.has(mp))) continue;
        for (const p of g.inner) { const d = t.d(p); if (d < best) { best = d; at = p; what = t.what; } }
      }
      if (isFinite(best)) worst('clearance', best, at, `line ${ar.line}: ${ar.text?.trim() ?? ''} — passes ${best.toFixed(1)} px from a ${what}`, g.samples);
      const sg = selfGap(g.samples);
      // only a curve that bends back on itself can touch itself; a plain arch's ends are just its chord
      const chord = len(sub(g.E, g.S));
      if (sg.d < chord - 1) worst('loop', sg.d, sg.at, `line ${ar.line}: ${ar.text?.trim() ?? ''} (chord ${chord.toFixed(0)} px)`, g.samples);
      else if (chord < 6) worst('loop', chord, g.S, `line ${ar.line}: ${ar.text?.trim() ?? ''} — head ends ${chord.toFixed(1)} px from its own start`, g.samples);
      worst('snake', inflections(g.samples), g.samples[Math.floor(g.samples.length / 2)], `line ${ar.line}: ${ar.text?.trim() ?? ''}`, g.samples);
    }

    for (let i = 0; i < withGeo.length; i++) for (let j = i + 1; j < withGeo.length; j++) {
      const A = withGeo[i].geo.inner; const B = withGeo[j].geo.inner;
      if (A.length < 2 || B.length < 2) continue;
      const d = polyDist(A, B);
      worst('tangle', d, A[Math.floor(A.length / 2)], `lines ${withGeo[i].line} and ${withGeo[j].line}`, [...A]);
    }
  });

  const r = res.report;
  m.crossings = r.arrowCrossings + r.bondCrossings + r.labelHits;
  for (const p of r.problems) {
    const g = p.ar.geo;
    issues.push({ metric: 'crossings', value: 1, at: g ? g.samples[Math.floor(g.samples.length / 2)] : null, path: g?.samples, note: `line ${p.ar.line}: ${p.what.join(', ')}` });
  }
  return { metrics: m, issues };
}

// issue count weighted by how far past the limit each one is; 0 = nothing to complain about
export function score(issues) {
  let s = 0;
  for (const it of issues) {
    const lim = LIMITS[it.metric];
    const over = lim.max != null ? (lim.max ? (it.value - lim.max) / lim.max : it.value) : (lim.min - it.value) / lim.min;
    s += 1 + Math.min(3, Math.max(0, over));
  }
  return Math.round(s * 100) / 100;
}

// ---------- main ----------
function loadItems() {
  const items = examples.map((e) => ({ name: `example/${e.name}`, code: e.code }));
  for (const dir of ['test/ai', 'test/cases', 'test/quality']) {
    if (!existsSync(dir)) continue;
    for (const f of readdirSync(dir).filter((x) => x.endsWith('.mech')).sort()) items.push({ name: `${dir.split('/')[1]}/${f.replace(/\.mech$/, '')}`, code: readFileSync(`${dir}/${f}`, 'utf8') });
  }
  return only ? items.filter((it) => it.name.includes(only)) : items;
}

const fmtV = (k, v) => (v == null ? '–' : LIMITS[k].unit === '°' || LIMITS[k].unit === 'px' ? `${v.toFixed(1)}${LIMITS[k].unit}` : LIMITS[k].unit === '' ? String(v) : `${v.toFixed(2)}${LIMITS[k].unit === 'L' ? 'L' : '×'}`);
const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

function annotate(svg, issues) {
  const marks = [];
  for (const it of issues) {
    const title = `<title>${esc(`${LIMITS[it.metric].label}: ${fmtV(it.metric, it.value)} — ${it.note}`)}</title>`;
    if (it.path && ['detour', 'miss', 'tangle', 'crossings', 'clearance', 'loop', 'snake'].includes(it.metric)) {
      marks.push(`<polyline points="${it.path.map((p) => `${p.x.toFixed(1)},${p.y.toFixed(1)}`).join(' ')}" fill="none" stroke="#f79009" stroke-width="5" stroke-opacity="0.45" stroke-linecap="round">${title}</polyline>`);
    }
    if (it.at) marks.push(`<circle cx="${it.at.x.toFixed(1)}" cy="${it.at.y.toFixed(1)}" r="13" fill="none" stroke="#d92d20" stroke-width="2" stroke-dasharray="4 3">${title}</circle>`);
  }
  return svg.replace(/<\/svg>$/, `<g class="qa-marks">${marks.join('')}</g></svg>`);
}

async function main() {
  const ms = createMechScript(await initRDKit());
  const L = 34;
  const results = [];
  for (const it of loadItems()) {
    const res = ms.toSVG(it.code, { color: '#111' });
    if (res.errors.length) { results.push({ name: it.name, error: res.errors.map((e) => `line ${e.line}: ${e.message}`).join('; '), metrics: {}, issues: [], score: 99 }); continue; }
    const { metrics, issues } = measure(res, L);
    results.push({ name: it.name, code: it.code, svg: annotate(res.svg, issues), metrics, issues, score: score(issues) });
  }

  const baseline = existsSync(BASELINE) ? JSON.parse(readFileSync(BASELINE, 'utf8')) : null;
  const base = new Map((baseline?.mechanisms || []).map((b) => [b.name, b]));
  const keys = Object.keys(LIMITS);
  const counts = Object.fromEntries(keys.map((k) => [k, 0]));
  for (const r of results) for (const it of r.issues) counts[it.metric]++;
  const total = results.reduce((s, r) => s + r.score, 0);

  // terminal table
  const pad = (s, n) => String(s).padEnd(n);
  console.log(pad('mechanism', 40) + pad('score', 8) + pad('was', 8) + 'issues');
  let worse = 0;
  for (const r of [...results].sort((a, b) => b.score - a.score)) {
    const was = base.get(r.name)?.score;
    const delta = was == null ? '' : r.score > was + 0.01 ? ' ▲' : r.score < was - 0.01 ? ' ▼' : '';
    if (was != null && r.score > was + 0.01) worse++;
    const kinds = {};
    for (const it of r.issues) kinds[it.metric] = (kinds[it.metric] || 0) + 1;
    console.log(pad(r.name, 40) + pad(r.score.toFixed(1) + delta, 8) + pad(was == null ? '' : was.toFixed(1), 8) + (r.error || Object.entries(kinds).map(([k, n]) => `${k}×${n}`).join(' ')));
  }
  const issueTotal = Object.values(counts).reduce((a, b) => a + b, 0);
  console.log(`\n${results.length} mechanisms · total score ${total.toFixed(1)}${baseline ? ` (baseline ${baseline.total.toFixed(1)})` : ''} · ${issueTotal} issues: ${keys.filter((k) => counts[k]).map((k) => `${k} ${counts[k]}`).join(', ')}`);

  // report
  const summary = { total: Math.round(total * 100) / 100, counts, mechanisms: results.map((r) => ({ name: r.name, score: r.score, metrics: r.metrics, issues: r.issues.map(({ metric, value, note }) => ({ metric, value: Math.round(value * 1000) / 1000, note })), error: r.error })) };
  mkdirSync('review', { recursive: true });
  writeFileSync('review/quality.json', JSON.stringify(summary, null, 1));
  writeFileSync('review/quality.html', page(results, summary, baseline, base));
  console.log('wrote review/quality.html and review/quality.json');

  if (args.includes('--save-baseline')) {
    if (only) { console.error('refusing to save a baseline from a partial run (--only)'); process.exit(2); }
    writeFileSync(BASELINE, JSON.stringify({ ...summary, mechanisms: summary.mechanisms.map(({ name, score, metrics }) => ({ name, score, metrics })) }, null, 1) + '\n');
    console.log(`saved ${BASELINE}`);
  }
  if (args.includes('--check') && worse) { console.error(`${worse} mechanism(s) scored worse than the baseline`); process.exit(1); }
}

function page(results, summary, baseline, base) {
  const keys = Object.keys(LIMITS);
  const sorted = [...results].sort((a, b) => b.score - a.score);
  const head = `<tr><th>mechanism</th><th>score</th>${keys.map((k) => `<th title="${esc(LIMITS[k].label)} (${LIMITS[k].max != null ? '>' : '<'} ${LIMITS[k].max ?? LIMITS[k].min}${LIMITS[k].unit})">${k}</th>`).join('')}</tr>`;
  const rows = sorted.map((r) => {
    const was = base.get(r.name)?.score;
    const d = was == null ? '' : r.score > was + 0.01 ? ` <span class="up">▲ ${was.toFixed(1)}</span>` : r.score < was - 0.01 ? ` <span class="down">▼ ${was.toFixed(1)}</span>` : '';
    return `<tr><td><a href="#${esc(r.name)}">${esc(r.name)}</a></td><td class="num">${r.score.toFixed(1)}${d}</td>${keys.map((k) => `<td class="num ${bad(k, r.metrics[k]) ? 'bad' : ''}">${r.error ? '' : fmtV(k, r.metrics[k])}</td>`).join('')}</tr>`;
  }).join('');
  const cards = sorted.map((r) => `<section id="${esc(r.name)}"><h2>${esc(r.name)} <span class="score">score ${r.score.toFixed(1)}</span></h2>
${r.error ? `<p class="bad">${esc(r.error)}</p>` : `<div class="fig">${r.svg}</div>`}
${r.issues.length ? `<ul>${r.issues.map((it) => `<li><b>${esc(LIMITS[it.metric].label)}</b> ${fmtV(it.metric, it.value)} — ${esc(it.note)}</li>`).join('')}</ul>` : '<p class="ok">no issues</p>'}
<details><summary>code</summary><pre>${esc(r.code || '')}</pre></details></section>`).join('');
  return `<!doctype html><html><head><meta charset="utf-8"><title>MechScript quality</title>
<style>body{font-family:system-ui,sans-serif;margin:20px;color:#111;background:#fff;max-width:1400px}
table{border-collapse:collapse;font-size:12px;margin:12px 0}th,td{padding:3px 7px;border-bottom:1px solid #eee;text-align:left}
th{position:sticky;top:0;background:#fff}.num{text-align:right;font-variant-numeric:tabular-nums}.bad{color:#b42318;font-weight:600}
.up{color:#b42318}.down{color:#067647}.ok{color:#067647}section{border-top:1px solid #ddd;padding:12px 0}h2{font-size:15px;margin:4px 0}
.score{font-weight:400;color:#666;font-size:13px}.fig svg{max-width:100%;height:auto}ul{font-size:13px;margin:6px 0;padding-left:18px}
pre{font-size:11px;background:#f6f6f6;padding:6px;overflow:auto}</style></head><body>
<h1>MechScript drawing quality</h1>
<p>${results.length} mechanisms · total score <b>${summary.total.toFixed(1)}</b>${baseline ? ` (baseline ${baseline.total.toFixed(1)})` : ''}. Worst first. Red dashed circles mark geometry issues, orange bands mark arrow issues (hover for details). Thresholds are in <code>tools/quality.js</code>.</p>
<table>${head}${rows}</table>${cards}</body></html>`;
}

if (import.meta.url === `file://${process.argv[1]}`) await main();
