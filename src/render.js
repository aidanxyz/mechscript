// MechScript renderer: parsed + analysed document -> SVG string.
// Pure string output (no DOM), so it runs in Node, in workers and in the browser.

import { placeArrows, arrowReport } from './arrows.js';

const DEFAULTS = {
  bondLength: 34,
  fontSize: 15,
  maxWidth: 960,
  arrowColor: '#d0342c',
  color: 'currentColor',
  background: null, // e.g. '#fff' for PNG export
  showTitle: true,
  showMaps: false,
  showLonePairs: 'arrows', // 'arrows' | 'all'
  fontFamily: 'Helvetica, Arial, sans-serif',
};

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const f = (n) => (Math.round(n * 100) / 100).toString();
const V = {
  add: (a, b) => ({ x: a.x + b.x, y: a.y + b.y }),
  sub: (a, b) => ({ x: a.x - b.x, y: a.y - b.y }),
  mul: (a, k) => ({ x: a.x * k, y: a.y * k }),
  len: (a) => Math.hypot(a.x, a.y),
  norm: (a) => { const l = Math.hypot(a.x, a.y) || 1; return { x: a.x / l, y: a.y / l }; },
  perp: (a) => ({ x: -a.y, y: a.x }), // left normal (y down => visually counter-clockwise)
  dot: (a, b) => a.x * b.x + a.y * b.y,
  fromAngle: (t) => ({ x: Math.cos(t), y: Math.sin(t) }),
  lerp: (a, b, t) => ({ x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t }),
};

// Rough Helvetica advance widths (em).
function textWidth(s, size) {
  let w = 0;
  for (const ch of s) {
    if (/[A-Z]/.test(ch)) w += ch === 'I' ? 0.28 : ch === 'M' || ch === 'W' ? 0.83 : 0.68;
    else if (/[a-z]/.test(ch)) w += 'ijl'.includes(ch) ? 0.24 : 'mw'.includes(ch) ? 0.8 : 0.54;
    else if (/[0-9]/.test(ch)) w += 0.556;
    else if (ch === ' ') w += 0.28;
    else w += 0.55;
  }
  return w * size;
}

const angDist = (a, b) => { let d = Math.abs(a - b) % (2 * Math.PI); return d > Math.PI ? 2 * Math.PI - d : d; };

// An arrow from one bond of an atom to another bond of the same atom (e.g. a C-H bond's electrons moving
// into a neighbouring C-C bond) has a clear path only if no other bond of that atom is drawn between the two.
// No stereochemistry is drawn, so end atoms (H, Br, a methyl...) around the shared atom may trade places;
// swap them when that leaves fewer bonds in the arrow's way.
function clearArrowPaths(st) {
  const where = new Map(); // map -> {sp, i}
  for (const sp of st.species) if (sp.mol) sp.mol.atoms.forEach((a, i) => { if (a.map != null) where.set(a.map, { sp, i }); });
  for (const ar of st.arrows) {
    if (ar.from.type !== 'bond' || ar.to.type !== 'bond') continue;
    const shared = ar.from.atoms.find((n) => ar.to.atoms.includes(n));
    if (shared == null) continue;
    const a = ar.from.atoms.find((n) => n !== shared);
    const b = ar.to.atoms.find((n) => n !== shared);
    const X = where.get(shared); const A = where.get(a); const B = where.get(b);
    if (!X || !A || !B || A.sp !== X.sp || B.sp !== X.sp) continue;
    const mol = X.sp.mol;
    const nb = (i) => mol.bonds.filter((bd) => bd.a === i || bd.b === i).map((bd) => (bd.a === i ? bd.b : bd.a));
    const around = nb(X.i);
    if (!around.includes(A.i) || !around.includes(B.i)) continue;
    const ang = (j) => Math.atan2(mol.atoms[j].y - mol.atoms[X.i].y, mol.atoms[j].x - mol.atoms[X.i].x);
    // bonds strictly inside the shorter arc between A and B
    const between = () => {
      const ta = ang(A.i); const tb = ang(B.i);
      let d = tb - ta; while (d <= -Math.PI) d += 2 * Math.PI; while (d > Math.PI) d -= 2 * Math.PI;
      return around.filter((j) => j !== A.i && j !== B.i).filter((j) => {
        let e = ang(j) - ta; while (e <= -Math.PI) e += 2 * Math.PI; while (e > Math.PI) e -= 2 * Math.PI;
        return d > 0 ? e > 0 && e < d : e < 0 && e > d;
      }).length;
    };
    const terminal = (j) => nb(j).length === 1;
    const swap = (i, j) => {
      const p = mol.atoms[i]; const q = mol.atoms[j];
      [p.x, q.x] = [q.x, p.x]; [p.y, q.y] = [q.y, p.y];
    };
    let best = between();
    for (const end of [A.i, B.i]) {
      if (!terminal(end) || best === 0) continue;
      for (const t of around) {
        if (t === A.i || t === B.i || !terminal(t)) continue;
        swap(end, t);
        const c = between();
        if (c < best) best = c; else swap(end, t);
      }
    }
  }
}

function facingInState(st) {
  const kOf = new Map();
  st.species.forEach((sp, k) => sp.atoms.forEach((a) => { if (a.map != null) kOf.set(a.map, k); }));
  const facing = new Map();
  for (const ar of st.arrows) {
    const ends = [...ar.from.atoms, ...ar.to.atoms];
    for (const p of ends) for (const q of ends) {
      const kp = kOf.get(p); const kq = kOf.get(q);
      if (kp == null || kq == null || kp === kq) continue;
      facing.set(p, (facing.get(p) || 0) + (kq > kp ? 1 : -1));
    }
  }
  for (const [mp, v] of facing) facing.set(mp, Math.sign(v));
  return facing;
}

export function render(doc, options = {}) {
  const o = { ...DEFAULTS, ...options };
  const L = o.bondLength;
  const fs = o.fontSize;

  // ---------- 1. per-species geometry ----------
  const lpSources = new Map(); // state index -> Map(map -> [targetMap...])
  doc.states.forEach((st, si) => {
    const m = new Map();
    for (const ar of st.arrows) {
      if (ar.from.type === 'atom') {
        const k = ar.from.atoms[0];
        if (!m.has(k)) m.set(k, []);
        m.get(k).push(ar);
      }
    }
    lpSources.set(si, m);
    // which way should each lone pair face? towards its target
    const speciesOfMap = new Map();
    st.species.forEach((sp, k) => sp.atoms.forEach((a) => { if (a.map != null) speciesOfMap.set(a.map, k); }));
    for (const ar of st.arrows) {
      if (ar.from.type !== 'atom') continue;
      const fromK = speciesOfMap.get(ar.from.atoms[0]);
      const tgt = ar.to.atoms.find((n) => n !== ar.from.atoms[0]) ?? ar.to.atoms[0];
      const toK = speciesOfMap.get(tgt);
      ar._targetMap = tgt;
      ar.preferDir = null;
      if (toK != null && fromK != null && toK !== fromK) ar.preferDir = V.norm({ x: toK > fromK ? 1 : -1, y: -0.8 });
    }
  });

  const blocks = []; // layout items
  doc.states.forEach((st, si) => {
    // molecules keep their orientation from the step before (matched by atom map numbers)
    clearArrowPaths(st);
    const refs = si > 0 ? doc.states[si - 1].species.filter((sp) => sp.mol && sp.geo && sp.geo.P).map((sp) => ({ P: sp.geo.P, atoms: sp.mol.atoms })) : [];
    // which side each reacting atom should face: towards the molecule its arrows connect it to. A molecule
    // drawn for the first time looks ahead to the step where it reacts (it may be a spectator here).
    const facing = new Map();
    for (let sj = si; sj < doc.states.length; sj++) {
      const f = facingInState(doc.states[sj]);
      for (const [mp, v] of f) if (!facing.has(mp)) facing.set(mp, v);
    }
    // author's orientation lines ("face 3 right", "flip 3"), per species
    const layoutOf = new Map();
    for (const lay of st.layout || []) {
      const k = st.species.findIndex((sp) => sp.atoms.some((a) => a.map === lay.atom));
      if (k < 0) continue;
      const cur = layoutOf.get(k) || {};
      if (lay.type === 'flip') cur.flip = true; else cur.face = { map: lay.atom, dir: lay.dir };
      layoutOf.set(k, cur);
    }
    const oState = { ...o, _facing: facing, _layout: layoutOf };
    st.species.forEach((sp, k) => {
      sp.geo = sp.mol ? speciesGeometry(sp.mol, oState, lpSources.get(si), si, k, refs) : errorGeometry(sp, o);
    });
  });

  // ---------- 2. flow layout ----------
  const units = []; // {kind:'species'|'plus'|'transition', w, h, ...}
  doc.states.forEach((st, si) => {
    if (si > 0) {
      const tr = doc.transitions[si - 1];
      const lw = tr && tr.label ? textWidth(tr.label, fs * 0.8) + 16 : 0;
      units.push({ kind: 'transition', tr, w: Math.max(L * 2.6, lw), h: fs * 2.4, breakBefore: true });
    }
    st.species.forEach((sp, k) => {
      if (k > 0) units.push({ kind: 'plus', w: fs * 1.6, h: fs, glueNext: true });
      const g = sp.geo;
      units.push({ kind: 'species', sp, si, k, w: g.box.x1 - g.box.x0, h: g.box.y1 - g.box.y0, breakBefore: k === 0 });
    });
  });

  const titleH = doc.title && o.showTitle ? fs * 2 : 0;
  const gap = L * 0.45;
  // Rows break before a reaction arrow, so each step's species stay together (their curved arrows join them).
  // A single step wider than maxWidth is split between species as a last resort.
  const groups = [];
  for (const u of units) {
    if (u.kind === 'transition' || !groups.length) groups.push([u]);
    else groups[groups.length - 1].push(u);
  }
  const widthOf = (us) => us.reduce((w, u, i) => w + u.w + (i ? gap : 0), 0);
  const rows = [];
  let row = [];
  const flush = () => { if (row.length) rows.push(row); row = []; };
  for (const grp of groups) {
    if (row.length && (grp[0].tr?.rowBreak || widthOf([...row, ...grp]) > o.maxWidth)) flush();
    if (widthOf(grp) <= o.maxWidth) { row.push(...grp); continue; }
    for (let i = 0; i < grp.length; i++) {
      const u = grp[i];
      const piece = u.kind === 'plus' && grp[i + 1] ? [u, grp[i + 1]] : [u];
      if (row.length && widthOf([...row, ...piece]) > o.maxWidth) flush();
      row.push(...piece);
      i += piece.length - 1;
    }
  }
  flush();

  const padArc = L * 0.9; // room above/below molecules for curved arrows
  let y = titleH;
  let maxRowW = 0;
  for (const r of rows) {
    const h = Math.max(...r.map((u) => u.h)) + padArc * 2;
    const mid = y + h / 2;
    let x = 0;
    r.forEach((u, idx) => {
      if (idx > 0) x += gap;
      u.x = x; u.mid = mid;
      if (u.kind === 'species') {
        const g = u.sp.geo;
        g.dx = x - g.box.x0;
        g.dy = mid - (g.box.y0 + g.box.y1) / 2;
      }
      x += u.w;
    });
    maxRowW = Math.max(maxRowW, x);
    y += h;
  }

  // ---------- 3. emit SVG ----------
  const parts = [];
  const bbox = { x0: 0, y0: 0, x1: Math.max(maxRowW, 10), y1: Math.max(y, 10) };
  const grow = (x, yy) => { bbox.x0 = Math.min(bbox.x0, x); bbox.x1 = Math.max(bbox.x1, x); bbox.y0 = Math.min(bbox.y0, yy); bbox.y1 = Math.max(bbox.y1, yy); };

  if (titleH) {
    parts.push(`<text x="0" y="${f(fs * 1.1)}" font-size="${f(fs * 1.05)}" font-weight="600" fill="${o.color}">${esc(doc.title)}</text>`);
    grow(textWidth(doc.title, fs * 1.05), 0);
  }

  // curved arrows: placed per step to avoid crossings (see arrows.js); they claim the lone pairs they start from
  const arrowSvgs = [];
  doc.states.forEach((st, si) => {
    const byMap = new Map();
    st.species.forEach((sp) => {
      if (!sp.mol) return;
      for (const a of sp.mol.atoms) if (a.map != null) byMap.set(a.map, { sp, a });
    });
    const row = rows.find((r) => r.some((u) => u.kind === 'species' && u.si === si)) || [];
    const geos = placeArrows(st, row, byMap, o);
    st.arrows.forEach((ar, ai) => {
      const geo = geos[ai];
      if (!geo) return;
      if (geo.dot) geo.dot.sp.geo.commitDot(geo.dot.i, geo.dot.kind, geo.dot.angle);
      ar.geo = geo;
      arrowSvgs.push(curvedArrowSvg(geo, ar, o, si, ai));
      for (const p of geo.samples) grow(p.x, p.y);
    });
    // move charge symbols out of the way of this step's arrows (incl. a point near each arrowhead's wings)
    const pts = st.arrows.flatMap((ar) => (ar.geo ? ar.geo.samples : []));
    for (const sp of st.species) if (sp.mol && sp.geo.relocateCharges) sp.geo.relocateCharges(pts.map((p) => ({ x: p.x - sp.geo.dx, y: p.y - sp.geo.dy })));
  });

  for (const u of units) {
    if (u.kind === 'species') parts.push(u.sp.geo.svg(u.sp.geo.dx, u.sp.geo.dy));
    else if (u.kind === 'plus') parts.push(`<text x="${f(u.x + u.w / 2)}" y="${f(u.mid + fs * 0.35)}" text-anchor="middle" font-size="${f(fs * 1.2)}" fill="${o.color}">+</text>`);
    else if (u.kind === 'transition') parts.push(transitionSvg(u, o));
  }
  parts.push(...arrowSvgs);

  // tighten bbox to drawn content (species boxes are padded)
  const content = { x0: Infinity, y0: Infinity, x1: -Infinity, y1: -Infinity };
  const cgrow = (x, yy) => { content.x0 = Math.min(content.x0, x); content.x1 = Math.max(content.x1, x); content.y0 = Math.min(content.y0, yy); content.y1 = Math.max(content.y1, yy); };
  if (titleH) { cgrow(0, 0); cgrow(textWidth(doc.title, fs * 1.05), titleH); }
  for (const u of units) {
    if (u.kind === 'species') { const g = u.sp.geo; const ik = g.inkFinal ? g.inkFinal() : g.ink; cgrow(ik.x0 + g.dx, ik.y0 + g.dy); cgrow(ik.x1 + g.dx, ik.y1 + g.dy); }
    else { cgrow(u.x, u.mid - u.h / 2); cgrow(u.x + u.w, u.mid + u.h / 2); }
  }
  doc.states.forEach((st) => st.arrows.forEach((ar) => { if (ar.geo) for (const p of ar.geo.samples) cgrow(p.x, p.y); }));
  const pad = 8;
  const vb = isFinite(content.x0) ? content : bbox;
  const W = vb.x1 - vb.x0 + pad * 2;
  const H = vb.y1 - vb.y0 + pad * 2;
  const bg = o.background ? `<rect x="${f(vb.x0 - pad)}" y="${f(vb.y0 - pad)}" width="${f(W)}" height="${f(H)}" fill="${o.background}"/>` : '';
  const defs = `<defs><marker id="ms-head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="${o.color}"/></marker></defs>`;
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" class="mechscript" viewBox="${f(vb.x0 - pad)} ${f(vb.y0 - pad)} ${f(W)} ${f(H)}" width="${f(W)}" height="${f(H)}" font-family="${esc(o.fontFamily)}" role="img" aria-label="${esc(doc.title || 'reaction mechanism')}">${defs}${bg}${parts.join('')}</svg>`;
  return { svg, width: W, height: H, report: arrowReport(doc, o) };
}

// ---------------- species ----------------

function errorGeometry(sp, o) {
  const fs = o.fontSize;
  const label = sp.smiles || '(empty)';
  const w = textWidth(label, fs * 0.85) + 16;
  const h = fs * 2;
  const box = { x0: 0, y0: 0, x1: w, y1: h };
  return {
    box, ink: box, atomPos: () => null,
    svg: (dx, dy) => `<g class="ms-error"><rect x="${f(dx)}" y="${f(dy)}" width="${f(w)}" height="${f(h)}" rx="4" fill="none" stroke="#d0342c" stroke-dasharray="4 3"/><text x="${f(dx + w / 2)}" y="${f(dy + h / 2 + fs * 0.3)}" text-anchor="middle" font-size="${f(fs * 0.85)}" fill="#d0342c">${esc(label)}</text></g>`,
  };
}

// Rotates (and if better, mirrors) P in place so atoms shared with a reference species line up with it.
function alignToReference(P, atoms, refs) {
  let best = null;
  for (const ref of refs) {
    const pairs = [];
    atoms.forEach((a, i) => {
      if (a.map == null) return;
      const j = ref.atoms.findIndex((b) => b.map === a.map);
      if (j >= 0) pairs.push([P[i], ref.P[j]]);
    });
    if (pairs.length >= 2 && (!best || pairs.length > best.length)) best = pairs;
  }
  if (!best) return false;
  const cA = best.reduce((c, [p]) => ({ x: c.x + p.x / best.length, y: c.y + p.y / best.length }), { x: 0, y: 0 });
  const cB = best.reduce((c, [, q]) => ({ x: c.x + q.x / best.length, y: c.y + q.y / best.length }), { x: 0, y: 0 });
  const fit = (mirror) => {
    let sxx = 0; let sxy = 0;
    for (const [p, q] of best) {
      const ax = p.x - cA.x; const ay = (p.y - cA.y) * (mirror ? -1 : 1);
      const bx = q.x - cB.x; const by = q.y - cB.y;
      sxx += ax * bx + ay * by; sxy += ax * by - ay * bx;
    }
    const th = Math.atan2(sxy, sxx);
    const c = Math.cos(th); const sn = Math.sin(th);
    let err = 0;
    for (const [p, q] of best) {
      const ax = p.x - cA.x; const ay = (p.y - cA.y) * (mirror ? -1 : 1);
      err += (c * ax - sn * ay - (q.x - cB.x)) ** 2 + (sn * ax + c * ay - (q.y - cB.y)) ** 2;
    }
    return { c, sn, mirror, err };
  };
  const a = fit(false); const b = fit(true);
  const t = b.err < a.err - 1e-6 ? b : a;
  for (const p of P) {
    const ax = p.x - cA.x; const ay = (p.y - cA.y) * (t.mirror ? -1 : 1);
    p.x = t.c * ax - t.sn * ay + cA.x;
    p.y = t.sn * ax + t.c * ay + cA.y;
  }
  return true;
}

// Turns (and if needed mirrors) a molecule so atoms that react with a molecule to its right face right and
// atoms that react with one to its left face left. facing[i] = +1 / -1 / 0. Small preference for RDKit's
// own orientation so unreactive molecules keep their usual look.
function orientForReaction(P, facing) {
  if (!facing.some((v) => v)) return;
  const cx = P.reduce((t, p) => t + p.x, 0) / P.length;
  const cy = P.reduce((t, p) => t + p.y, 0) / P.length;
  const R = Math.max(1, ...P.map((p) => Math.hypot(p.x - cx, p.y - cy)));
  let best = null;
  for (const mirror of [false, true]) {
    for (let k = 0; k < 24; k++) {
      const th = (k / 24) * 2 * Math.PI;
      const c = Math.cos(th); const sn = Math.sin(th);
      let score = 0;
      P.forEach((p, i) => {
        if (!facing[i]) return;
        const ax = p.x - cx; const ay = (p.y - cy) * (mirror ? -1 : 1);
        score += facing[i] * (c * ax - sn * ay) / R;
      });
      score -= 0.08 * Math.min(th, 2 * Math.PI - th) + (mirror ? 0.05 : 0);
      if (!best || score > best.score + 1e-9) best = { score, c, sn, mirror };
    }
  }
  for (const p of P) {
    const ax = p.x - cx; const ay = (p.y - cy) * (best.mirror ? -1 : 1);
    p.x = best.c * ax - best.sn * ay + cx;
    p.y = best.sn * ax + best.c * ay + cy;
  }
}

// "flip N": mirror the molecule; "face N dir": turn it so atom N points right/left/up/down from its centre.
function applyAuthorLayout(P, atoms, lay) {
  const cx = P.reduce((t, p) => t + p.x, 0) / P.length;
  const cy = P.reduce((t, p) => t + p.y, 0) / P.length;
  if (lay.flip) for (const p of P) p.y = 2 * cy - p.y;
  if (lay.face) {
    const i = atoms.findIndex((a) => a.map === lay.face.map);
    if (i < 0) return;
    const vx = P[i].x - cx; const vy = P[i].y - cy;
    if (Math.hypot(vx, vy) < 1e-6) return;
    const target = { right: 0, down: Math.PI / 2, left: Math.PI, up: -Math.PI / 2 }[lay.face.dir];
    const th = target - Math.atan2(vy, vx);
    const c = Math.cos(th); const sn = Math.sin(th);
    for (const p of P) { const x = p.x - cx; const y = p.y - cy; p.x = cx + c * x - sn * y; p.y = cy + sn * x + c * y; }
  }
}

function speciesGeometry(mol, o, lpSourceMap, si, k, refs = []) {
  const L = o.bondLength;
  const fs = o.fontSize;
  const { atoms, bonds } = mol;

  // scale RDKit coordinates to bond length L
  const lens = bonds.map((b) => Math.hypot(atoms[b.a].x - atoms[b.b].x, atoms[b.a].y - atoms[b.b].y)).sort((p, q) => p - q);
  const s = lens.length ? L / (lens[Math.floor(lens.length / 2)] || 1) : L;
  const P = atoms.map((a) => ({ x: a.x * s, y: a.y * s }));
  const aligned = refs.length ? alignToReference(P, atoms, refs) : false;
  if (!aligned && o._facing) orientForReaction(P, atoms.map((a) => (a.map != null ? o._facing.get(a.map) || 0 : 0)));
  const lay = o._layout && o._layout.get(k);
  if (lay) applyAuthorLayout(P, atoms, lay);

  const neighbors = atoms.map((a) => a.bonds.map((bi) => (bonds[bi].a === a.index ? bonds[bi].b : bonds[bi].a)));
  const heavyCount = atoms.filter((a) => a.symbol !== 'H').length;

  // ---- labels ----
  const labels = atoms.map((a, i) => {
    const deg = neighbors[i].length;
    let show = a.symbol !== 'C' || deg === 0 || (deg === 1 && heavyCount <= 3) || (a.radicals > 0 && deg <= 1);
    if (a.symbol === 'C' && a.charge !== 0 && deg <= 1) show = true;
    if (!show) return null;
    let hSide = 'right';
    if (a.hCount > 0) {
      if (deg === 0) hSide = (a.symbol === 'O' || a.symbol === 'S') && a.hCount === 2 ? 'left' : 'right';
      else {
        // put the H part on the side with more room (away from bonds)
        const angs = neighbors[i].map((j) => Math.atan2(P[j].y - P[i].y, P[j].x - P[i].x));
        const dR = Math.min(...angs.map((q) => angDist(q, 0)));
        const dL = Math.min(...angs.map((q) => angDist(q, Math.PI)));
        if (Math.abs(dR - dL) > 0.15) hSide = dR > dL ? 'right' : 'left';
        else {
          const sx = neighbors[i].reduce((acc, j) => acc + (P[j].x - P[i].x), 0) / deg;
          hSide = sx > L * 0.1 ? 'left' : 'right';
        }
      }
    }
    const symW = textWidth(a.symbol, fs);
    const hW = a.hCount > 0 ? textWidth('H', fs) + (a.hCount > 1 ? textWidth(String(a.hCount), fs * 0.7) : 0) : 0;
    return { hSide, symW, hW, r: Math.max(symW / 2, fs * 0.42) + 1.5 };
  });

  // occupied directions per atom (angles)
  const occupied = atoms.map((a, i) => {
    const angs = neighbors[i].map((j) => Math.atan2(P[j].y - P[i].y, P[j].x - P[i].x));
    if (labels[i] && a.hCount > 0) angs.push(labels[i].hSide === 'right' ? 0 : Math.PI);
    return angs;
  });

  const freeDir = (i, prefer = null, extra = []) => {
    const occ = [...occupied[i], ...extra];
    let best = null;
    for (let n = 0; n < 36; n++) {
      const t = (n / 36) * 2 * Math.PI - Math.PI;
      const minD = occ.length ? Math.min(...occ.map((q) => angDist(q, t))) : Math.PI;
      let score = Math.min(minD, 1.0);
      if (prefer) score += 0.8 * V.dot(V.fromAngle(t), prefer);
      else score += 0.05 * Math.cos(t + Math.PI / 4); // mild preference for upper right
      if (!best || score > best.score) best = { t, score, minD };
    }
    return best.t;
  };

  // ---- charges ----
  const labelCornerPos = (i, t) => {
    const a = atoms[i]; const lab = labels[i];
    const c = Math.cos(t); const sn = Math.sin(t);
    let half = lab.symW / 2;
    if (a.hCount && ((c > 0.1 && lab.hSide === 'right') || (c < -0.1 && lab.hSide === 'left'))) half += lab.hW;
    return { x: P[i].x + (Math.abs(c) < 0.1 ? 0 : Math.sign(c) * (half + fs * 0.3)), y: P[i].y + Math.sign(sn) * fs * (Math.abs(c) < 0.1 ? 0.95 : 0.55) };
  };
  const bondAngles = (i) => neighbors[i].map((j) => Math.atan2(P[j].y - P[i].y, P[j].x - P[i].x));
  const chargeAt = atoms.map((a, i) => {
    if (!a.charge) return null;
    let pos;
    if (labels[i]) {
      const lab = labels[i];
      // try the four label corners, then top/bottom; take the one furthest from bonds and the H label
      const cands = [-Math.PI / 4, -3 * Math.PI / 4, Math.PI / 4, 3 * Math.PI / 4, -Math.PI / 2, Math.PI / 2];
      let best = null;
      cands.forEach((t, rank) => {
        const minD = occupied[i].length ? Math.min(...occupied[i].map((q) => angDist(q, t))) : Math.PI;
        const score = Math.min(minD, 1.2) - rank * 0.01;
        if (!best || score > best.score) best = { t, score };
      });
      const t = best.t;
      const c = Math.cos(t); const sn = Math.sin(t);
      let half = lab.symW / 2;
      if (a.hCount && ((c > 0.1 && lab.hSide === 'right') || (c < -0.1 && lab.hSide === 'left'))) half += lab.hW;
      pos = { x: P[i].x + (Math.abs(c) < 0.1 ? 0 : Math.sign(c) * (half + fs * 0.3)), y: P[i].y + Math.sign(sn) * fs * (Math.abs(c) < 0.1 ? 0.95 : 0.55) };
      occupied[i].push(t);
    } else {
      const t = freeDir(i);
      occupied[i].push(t);
      pos = V.add(P[i], V.mul(V.fromAngle(t), L * 0.42));
    }
    return pos;
  });

  // After the curved arrows are placed: a charge symbol an arrow runs into (often the only free gap next to
  // a carbocation is also where its arrow arrives) moves to another free spot around its atom.
  const relocateCharges = (points) => {
    const avoid = [...points, ...lpDots.map((d) => d.center)];
    const clearance = (p) => (avoid.length ? Math.min(...avoid.map((q) => V.len(V.sub(p, q)))) : Infinity);
    atoms.forEach((a, i) => {
      if (!chargeAt[i] || clearance(chargeAt[i]) > 9) return;
      const cands = [];
      if (labels[i]) {
        for (const t of [-Math.PI / 4, -3 * Math.PI / 4, Math.PI / 4, 3 * Math.PI / 4, -Math.PI / 2, Math.PI / 2]) cands.push({ p: labelCornerPos(i, t), t });
      } else {
        for (let n = 0; n < 36; n++) for (const r of [0.42, 0.55]) { const t = (n / 36) * 2 * Math.PI - Math.PI; cands.push({ p: V.add(P[i], V.mul(V.fromAngle(t), L * r)), t, r }); }
      }
      const bonds = bondAngles(i);
      let best = null;
      for (const c of cands) {
        const bondGap = bonds.length ? Math.min(...bonds.map((q) => angDist(q, c.t))) : Math.PI;
        if (bondGap < 0.45) continue; // never on top of a bond
        const score = Math.min(clearance(c.p), 14) + 3 * Math.min(bondGap, 1) - (c.r === 0.55 ? 1 : 0);
        if (!best || score > best.score) best = { ...c, score };
      }
      if (best && clearance(best.p) > clearance(chargeAt[i])) chargeAt[i] = best.p;
    });
  };

  // ---- lone pairs + radicals ----
  // Lone pairs that start arrows are claimed by the arrow stage (after layout), so they face the arrow's curve.
  // finalizeDots() then places any remaining lone pairs (showLonePairs: 'all') and radicals.
  const lpDots = []; // {i, center, dir, kind}
  const claimed = atoms.map(() => ({ lp: 0, rad: 0 }));
  const dotBase = (i) => (labels[i] ? labels[i].r + 4 : L * 0.28);
  const placeDot = (i, kind, prefer) => {
    const t = freeDir(i, prefer);
    occupied[i].push(t);
    const dir = V.fromAngle(t);
    const d = { i, center: V.add(P[i], V.mul(dir, dotBase(i))), dir, kind };
    lpDots.push(d);
    return d;
  };
  // pure: where a dot facing `prefer` would go (nothing is reserved)
  const dotFor = (i, prefer) => {
    const t = freeDir(i, prefer);
    const dir = V.fromAngle(t);
    return { t, dir, center: V.add(P[i], V.mul(dir, dotBase(i))) };
  };
  const commitDot = (i, kind, t) => {
    claimed[i][kind]++;
    occupied[i].push(t);
    const dir = V.fromAngle(t);
    lpDots.push({ i, center: V.add(P[i], V.mul(dir, dotBase(i))), dir, kind });
  };
  let finalized = false;
  const finalizeDots = () => {
    if (finalized) return;
    finalized = true;
    atoms.forEach((a, i) => {
      const wantLp = o.showLonePairs === 'all' ? a.lonePairs : 0;
      for (let q = claimed[i].lp; q < wantLp; q++) placeDot(i, 'lp', null);
      for (let q = claimed[i].rad; q < a.radicals; q++) placeDot(i, 'rad', null);
    });
  };
  // atoms that may get dots, for the layout box
  const mayDot = atoms.map((a, i) => a.radicals > 0 || (o.showLonePairs === 'all' && a.lonePairs > 0) ||
    (a.map != null && lpSourceMap && (lpSourceMap.get(a.map) || []).length > 0));

  // ---- ink box ----
  const ink = { x0: Infinity, y0: Infinity, x1: -Infinity, y1: -Infinity };
  const g = (x, y) => { ink.x0 = Math.min(ink.x0, x); ink.x1 = Math.max(ink.x1, x); ink.y0 = Math.min(ink.y0, y); ink.y1 = Math.max(ink.y1, y); };
  atoms.forEach((a, i) => {
    const lab = labels[i];
    if (lab) {
      const left = lab.hSide === 'left' && a.hCount ? lab.symW / 2 + lab.hW : lab.symW / 2;
      const right = lab.hSide === 'right' && a.hCount ? lab.symW / 2 + lab.hW : lab.symW / 2;
      g(P[i].x - left, P[i].y - fs * 0.6); g(P[i].x + right, P[i].y + fs * 0.6);
    } else { g(P[i].x - 2, P[i].y - 2); g(P[i].x + 2, P[i].y + 2); }
    if (chargeAt[i]) { g(chargeAt[i].x - 6, chargeAt[i].y - 6); g(chargeAt[i].x + 6, chargeAt[i].y + 6); }
  });
  const box = { ...ink };
  atoms.forEach((a, i) => {
    if (!mayDot[i]) return;
    const r = dotBase(i) + 5;
    box.x0 = Math.min(box.x0, P[i].x - r); box.x1 = Math.max(box.x1, P[i].x + r);
    box.y0 = Math.min(box.y0, P[i].y - r); box.y1 = Math.max(box.y1, P[i].y + r);
  });
  const inkFinal = () => {
    const r = { ...ink };
    for (const d of lpDots) {
      r.x0 = Math.min(r.x0, d.center.x - 5); r.x1 = Math.max(r.x1, d.center.x + 5);
      r.y0 = Math.min(r.y0, d.center.y - 5); r.y1 = Math.max(r.y1, d.center.y + 5);
    }
    return r;
  };

  // ---- svg ----
  const trim = (i) => (labels[i] ? labels[i].r : 0);
  const svgFn = (dx, dy) => {
    finalizeDots();
    const out = [];
    const T = (p) => ({ x: p.x + dx, y: p.y + dy });
    const line = (p, q, extra = '') => `<line x1="${f(p.x)}" y1="${f(p.y)}" x2="${f(q.x)}" y2="${f(q.y)}"${extra}/>`;
    out.push(`<g class="ms-species" data-state="${si}" data-species="${k}">`);
    out.push(`<g stroke="${o.color}" stroke-width="1.6" stroke-linecap="round">`);
    for (const b of bonds) {
      const pa = T(P[b.a]); const pb = T(P[b.b]);
      const dir = V.norm(V.sub(pb, pa));
      const A = V.add(pa, V.mul(dir, trim(b.a)));
      const B = V.sub(pb, V.mul(dir, trim(b.b)));
      const nrm = V.perp(dir);
      if (b.order === 1) out.push(line(A, B));
      else if (b.order === 3) {
        out.push(line(A, B));
        out.push(line(V.add(A, V.mul(nrm, 4)), V.add(B, V.mul(nrm, 4))));
        out.push(line(V.sub(A, V.mul(nrm, 4)), V.sub(B, V.mul(nrm, 4))));
      } else {
        const side = doubleBondSide(mol, b, P, neighbors, labels);
        if (side === 0) {
          out.push(line(V.add(A, V.mul(nrm, 2.6)), V.add(B, V.mul(nrm, 2.6))));
          out.push(line(V.sub(A, V.mul(nrm, 2.6)), V.sub(B, V.mul(nrm, 2.6))));
        } else {
          out.push(line(A, B));
          const off = V.mul(nrm, 5.5 * side);
          const sh = L * 0.14;
          const A2 = V.add(V.add(A, off), V.mul(dir, labels[b.a] ? 0 : sh));
          const B2 = V.sub(V.add(B, off), V.mul(dir, labels[b.b] ? 0 : sh));
          out.push(line(A2, B2));
        }
      }
    }
    out.push('</g>');
    atoms.forEach((a, i) => {
      const p = T(P[i]);
      const lab = labels[i];
      const attrs = `class="ms-atom" data-state="${si}" data-species="${k}" data-atom="${i}"${a.map != null ? ` data-map="${a.map}"` : ''}`;
      out.push(`<g ${attrs}>`);
      if (lab) {
        const ty = p.y + fs * 0.36;
        out.push(`<text x="${f(p.x)}" y="${f(ty)}" text-anchor="middle" font-size="${fs}" fill="${o.color}">${esc(a.symbol)}</text>`);
        if (a.hCount > 0) {
          const sub = a.hCount > 1 ? `<tspan font-size="${f(fs * 0.7)}" dy="${f(fs * 0.22)}">${a.hCount}</tspan>` : '';
          if (lab.hSide === 'right') out.push(`<text x="${f(p.x + lab.symW / 2)}" y="${f(ty)}" font-size="${fs}" fill="${o.color}">H${sub}</text>`);
          else out.push(`<text x="${f(p.x - lab.symW / 2)}" y="${f(ty)}" text-anchor="end" font-size="${fs}" fill="${o.color}">H${sub}</text>`);
        }
      }
      if (chargeAt[i]) {
        const c = T(chargeAt[i]);
        const mag = Math.abs(a.charge);
        const sign = a.charge > 0 ? '+' : '−';
        if (mag === 1) {
          const r = 5.2;
          out.push(`<circle cx="${f(c.x)}" cy="${f(c.y)}" r="${r}" fill="none" stroke="${o.color}" stroke-width="1"/>`);
          out.push(`<line x1="${f(c.x - 3)}" y1="${f(c.y)}" x2="${f(c.x + 3)}" y2="${f(c.y)}" stroke="${o.color}" stroke-width="1.2"/>`);
          if (a.charge > 0) out.push(`<line x1="${f(c.x)}" y1="${f(c.y - 3)}" x2="${f(c.x)}" y2="${f(c.y + 3)}" stroke="${o.color}" stroke-width="1.2"/>`);
        } else out.push(`<text x="${f(c.x)}" y="${f(c.y + 4)}" text-anchor="middle" font-size="${f(fs * 0.7)}" fill="${o.color}">${mag}${sign}</text>`);
      }
      if (o.showMaps && a.map != null) {
        out.push(`<text x="${f(p.x + 7)}" y="${f(p.y + 13)}" font-size="${f(fs * 0.55)}" fill="#3b7dd8">${a.map}</text>`);
      }
      out.push('</g>');
    });
    for (const d of lpDots) {
      const c = T(d.center);
      if (d.kind === 'rad') { out.push(`<circle cx="${f(c.x)}" cy="${f(c.y)}" r="1.9" fill="${o.color}"/>`); continue; }
      const t = V.perp(d.dir);
      const p1 = V.add(c, V.mul(t, 3.3));
      const p2 = V.sub(c, V.mul(t, 3.3));
      out.push(`<circle cx="${f(p1.x)}" cy="${f(p1.y)}" r="1.7" fill="${o.color}"/><circle cx="${f(p2.x)}" cy="${f(p2.y)}" r="1.7" fill="${o.color}"/>`);
    }
    out.push('</g>');
    return out.join('');
  };

  const geo = {
    box, ink, inkFinal, svg: svgFn, dx: 0, dy: 0, labels, P, dotFor, commitDot, dotBase, relocateCharges,
    // things a curved arrow should not run through, in species coordinates
    obstacles: () => {
      const out = { bonds: [], labels: [], charges: [], vertices: [] };
      for (const b of bonds) {
        const dir = V.norm(V.sub(P[b.b], P[b.a]));
        out.bonds.push({ a: V.add(P[b.a], V.mul(dir, trim(b.a) + 1)), b: V.sub(P[b.b], V.mul(dir, trim(b.b) + 1)), atoms: [b.a, b.b] });
      }
      atoms.forEach((a, i) => {
        const lab = labels[i];
        if (lab) {
          // text boxes: half width / half height of the glyphs (cap height is about 0.72 em)
          out.labels.push({ c: P[i], r: lab.r, hw: lab.symW / 2 + 1.5, hh: fs * 0.38 + 1.5, atom: i });
          if (a.hCount) {
            const sgn = lab.hSide === 'right' ? 1 : -1;
            out.labels.push({ c: { x: P[i].x + sgn * (lab.symW / 2 + lab.hW / 2), y: P[i].y + fs * 0.05 }, r: lab.hW / 2 + 2, hw: lab.hW / 2 + 1.5, hh: fs * 0.45 + 1.5, atom: i });
          }
        } else out.vertices.push({ c: P[i], atom: i });
        if (chargeAt[i]) out.charges.push({ c: chargeAt[i], r: 6, atom: i });
      });
      for (const d of lpDots) out.charges.push({ c: d.center, r: 5, atom: d.i });
      out.rings = (mol.rings || []).filter((r) => r.length >= 3 && r.length <= 8).map((r) => ({ atoms: r, poly: r.map((i) => P[i]) }));
      return out;
    },
    atomPos: (i) => ({ x: P[i].x + geo.dx, y: P[i].y + geo.dy }),
    labelRadius: (i) => (labels[i] ? labels[i].r : 0),
    centroid: () => {
      const c = P.reduce((acc, p) => ({ x: acc.x + p.x, y: acc.y + p.y }), { x: 0, y: 0 });
      return { x: c.x / P.length + geo.dx, y: c.y / P.length + geo.dy };
    },
  };
  return geo;
}

// +1 / -1: which side the second line of a double bond goes on (in the bond's left-normal frame); 0 = centred.
function doubleBondSide(mol, b, P, neighbors, labels) {
  const ring = mol.rings.find((r) => r.includes(b.a) && r.includes(b.b));
  const dir = V.norm(V.sub(P[b.b], P[b.a]));
  const nrm = V.perp(dir);
  if (ring) {
    const c = ring.reduce((acc, i) => ({ x: acc.x + P[i].x / ring.length, y: acc.y + P[i].y / ring.length }), { x: 0, y: 0 });
    return V.dot(V.sub(c, P[b.a]), nrm) >= 0 ? 1 : -1;
  }
  const others = [...neighbors[b.a].filter((j) => j !== b.b), ...neighbors[b.b].filter((j) => j !== b.a)];
  if (neighbors[b.a].length === 1 || neighbors[b.b].length === 1 || labels[b.a] || labels[b.b]) {
    if (!(labels[b.a] && labels[b.b]) && (neighbors[b.a].length === 1 && neighbors[b.b].length >= 3 || neighbors[b.b].length === 1 && neighbors[b.a].length >= 3)) return 0;
    if (neighbors[b.a].length === 1 || neighbors[b.b].length === 1) return 0;
  }
  if (!others.length) return 0;
  const sum = others.reduce((acc, j) => acc + V.dot(V.sub(P[j], P[b.a]), nrm), 0);
  return sum >= 0 ? 1 : -1;
}

// ---------------- reaction arrows ----------------

function transitionSvg(u, o) {
  const fs = o.fontSize;
  const x0 = u.x + 4;
  const x1 = u.x + u.w - 4;
  const y = u.mid;
  const st = `stroke="${o.color}" stroke-width="1.6" fill="none"`;
  let s = '';
  if (u.tr?.type === 'equilibrium') {
    s += `<path d="M${f(x0)},${f(y - 3)} H${f(x1)} l-8,-5" ${st}/>`;
    s += `<path d="M${f(x1)},${f(y + 3)} H${f(x0)} l8,5" ${st}/>`;
  } else if (u.tr?.type === 'resonance') {
    s += `<line x1="${f(x0 + 2)}" y1="${f(y)}" x2="${f(x1 - 2)}" y2="${f(y)}" ${st} marker-start="url(#ms-head)" marker-end="url(#ms-head)"/>`;
  } else {
    s += `<line x1="${f(x0)}" y1="${f(y)}" x2="${f(x1 - 2)}" y2="${f(y)}" ${st} marker-end="url(#ms-head)"/>`;
  }
  if (u.tr?.label) s += `<text x="${f((x0 + x1) / 2)}" y="${f(y - 9)}" text-anchor="middle" font-size="${f(fs * 0.8)}" fill="${o.color}">${esc(u.tr.label)}</text>`;
  return `<g class="ms-transition" data-line="${u.tr?.line ?? ''}">${s}</g>`;
}

// ---------------- curved arrows ----------------

function curvedArrowSvg(geo, ar, o, si, ai) {
  const { S, E } = geo;
  const C1 = geo.C1 || V.lerp(S, geo.C, 2 / 3);
  const C2 = geo.C2 || V.lerp(E, geo.C, 2 / 3);
  const dir = V.norm(V.sub(E, C2));
  const head = 9;
  const base = V.sub(E, V.mul(dir, head));
  // pull the curve end back to the head's base so the line doesn't poke through the tip
  const Eline = V.sub(E, V.mul(dir, head * 0.6));
  const w = V.perp(dir);
  const col = o.arrowColor;
  const mid = { x: (S.x + 3 * C1.x + 3 * C2.x + E.x) / 8, y: (S.y + 3 * C1.y + 3 * C2.y + E.y) / 8 };
  let headSvg;
  if (ar.kind === 'fishhook') {
    // single barb on the outside of the curve
    const outward = V.dot(V.sub(mid, V.lerp(S, E, 0.5)), w) >= 0 ? 1 : -1;
    const barb = V.add(base, V.mul(w, 4.5 * outward));
    headSvg = `<path d="M${f(E.x)},${f(E.y)} L${f(barb.x)},${f(barb.y)} L${f(Eline.x)},${f(Eline.y)} z" fill="${col}"/>`;
  } else {
    const l = V.add(base, V.mul(w, 4.2));
    const r = V.sub(base, V.mul(w, 4.2));
    headSvg = `<path d="M${f(E.x)},${f(E.y)} L${f(l.x)},${f(l.y)} L${f(Eline.x)},${f(Eline.y)} L${f(r.x)},${f(r.y)} z" fill="${col}"/>`;
  }
  // end the line at the head's base, keeping the final tangent
  const C2line = V.add(Eline, V.sub(C2, E));
  const d = `M${f(S.x)},${f(S.y)} C${f(C1.x)},${f(C1.y)} ${f(C2line.x)},${f(C2line.y)} ${f(Eline.x)},${f(Eline.y)}`;
  let s = `<g class="ms-arrow" data-state="${si}" data-arrow="${ai}" data-line="${ar.line}">`;
  s += `<path d="${d}" stroke="${col}" stroke-width="1.5" fill="none"/>${headSvg}`;
  return s + '</g>';
}
