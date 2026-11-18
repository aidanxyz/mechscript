// Curved-arrow placement. For each step, every arrow gets a set of candidate curves (bow side, bow size,
// where the peak sits); each candidate is scored against the drawing (bonds crossed, labels/charges hit,
// length) and every pair of candidates against each other (crossings, near-overlaps, heads touching).
// The cheapest combination wins: exhaustive search when small, iterated improvement otherwise.
// Authors can steer an arrow with placement words ("above", "below", "left", "right", "tight", "wide").

import { findBond } from './chem.js';

const V = {
  add: (a, b) => ({ x: a.x + b.x, y: a.y + b.y }),
  sub: (a, b) => ({ x: a.x - b.x, y: a.y - b.y }),
  mul: (a, k) => ({ x: a.x * k, y: a.y * k }),
  len: (a) => Math.hypot(a.x, a.y),
  norm: (a) => { const l = Math.hypot(a.x, a.y) || 1; return { x: a.x / l, y: a.y / l }; },
  perp: (a) => ({ x: -a.y, y: a.x }),
  dot: (a, b) => a.x * b.x + a.y * b.y,
  lerp: (a, b, t) => ({ x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t }),
};

export const COST = {
  arrowCross: 120, // two curved arrows crossing
  arrowNear: 25, // two curved arrows running within a few px of each other
  headsTouch: 40, // arrowheads (or a head and a tail) on top of each other
  bondCross: 70, // a curve crossing a bond
  labelHit: 70, // a curve running through an atom label
  chargeHit: 35, // through a charge symbol or lone-pair dots
  atomGraze: 12, // passing right over an unlabelled atom (a skeletal vertex)
  furniture: 80, // through a "+" sign or a reaction arrow or its label
  length: 0.03, // per px of curve
  offSide: 4, // bowing towards the molecule instead of away from it
  dotClash: 90, // two arrows starting from the same lone-pair slot
  ringFace: 60, // cutting across a ring that the arrow does not move electrons around
  degenerate: 300, // curve too short to arch, or arrowhead not following the curve
  cramped: 150, // the curve comes back close to itself, or its head lands next to its own tail
  snake: 150, // the curve bends one way and then the other
  near: 45, // passing closer than CLEAR px to something the arrow is not about (full cost when touching)
};
const CLEAR = 8;

const SCALES = [0.55, 0.85, 1.25, 1.75, 2.6];
const PEAKS = [0.5, 0.32, 0.68, 0.18, 0.82];
const N = 20; // polyline samples per curve

// ---------- geometry helpers ----------

function segIntersect(p1, p2, p3, p4) {
  const d = (p2.x - p1.x) * (p4.y - p3.y) - (p2.y - p1.y) * (p4.x - p3.x);
  if (Math.abs(d) < 1e-9) return false;
  const u = ((p3.x - p1.x) * (p4.y - p3.y) - (p3.y - p1.y) * (p4.x - p3.x)) / d;
  const v = ((p3.x - p1.x) * (p2.y - p1.y) - (p3.y - p1.y) * (p2.x - p1.x)) / d;
  return u > 0 && u < 1 && v > 0 && v < 1;
}
// inside a label's text box, shrunk (slack > 0) or grown (slack < 0) by `slack` px
export function inLabel(p, l, slack = 0) {
  return Math.abs(p.x - l.c.x) < l.hw - slack && Math.abs(p.y - l.c.y) < l.hh - slack;
}
function inPoly(p, poly) {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const a = poly[i]; const b = poly[j];
    if ((a.y > p.y) !== (b.y > p.y) && p.x < ((b.x - a.x) * (p.y - a.y)) / (b.y - a.y) + a.x) inside = !inside;
  }
  return inside;
}
function distPointSeg(p, a, b) {
  const ab = V.sub(b, a);
  const l2 = ab.x * ab.x + ab.y * ab.y || 1;
  const t = Math.max(0, Math.min(1, V.dot(V.sub(p, a), ab) / l2));
  return V.len(V.sub(p, V.lerp(a, b, t)));
}
function cubic(S, C1, C2, E, u) {
  const v = 1 - u;
  return {
    x: v * v * v * S.x + 3 * v * v * u * C1.x + 3 * v * u * u * C2.x + u * u * u * E.x,
    y: v * v * v * S.y + 3 * v * v * u * C1.y + 3 * v * u * u * C2.y + u * u * u * E.y,
  };
}
function quad(S, C, E, u) {
  return {
    x: (1 - u) * (1 - u) * S.x + 2 * (1 - u) * u * C.x + u * u * E.x,
    y: (1 - u) * (1 - u) * S.y + 2 * (1 - u) * u * C.y + u * u * E.y,
  };
}

// ---------- endpoints (independent of the curve shape) ----------

export function arrowEndpoints(ar, byMap, o) {
  const L = o.bondLength;
  const get = (n) => byMap.get(n);
  if ([...ar.from.atoms, ...ar.to.atoms].some((n) => !get(n))) return null;
  const pos = (n) => { const { sp, a } = get(n); return sp.geo.atomPos(a.index); };

  const ep = { ar, involved: new Set() };
  if (ar.from.type === 'atom') {
    const { sp, a } = get(ar.from.atoms[0]);
    ep.S0 = pos(ar.from.atoms[0]);
    ep.dot = { sp, i: a.index, kind: ar.from.source === 'rad' ? 'rad' : 'lp' };
    ep.involved.add(sp);
  } else {
    const [p, q] = ar.from.atoms.map(pos);
    ep.S0 = V.lerp(p, q, 0.5);
    ep.bondEnds = [p, q];
    ep.bondNormal = V.perp(V.norm(V.sub(q, p)));
    ar.from.atoms.forEach((n) => ep.involved.add(get(n).sp));
  }
  ep.endRadius = 3;
  if (ar.to.type === 'atom') {
    const t = get(ar.to.atoms[0]);
    ep.T = pos(ar.to.atoms[0]);
    ep.endRadius = Math.max(t.sp.geo.labelRadius(t.a.index) + 3, L * 0.22);
    ep.targetAtom = t;
    ep.involved.add(t.sp);
  } else {
    const [x, y] = ar.to.atoms;
    const gx = get(x); const gy = get(y);
    ep.involved.add(gx.sp); ep.involved.add(gy.sp);
    const exists = gx.sp === gy.sp && findBond(gx.sp.mol, gx.a.index, gy.a.index);
    const px = pos(x); const py = pos(y);
    if (exists) ep.T = V.lerp(px, py, 0.5);
    else if (ar.kind === 'fishhook') ep.T = V.lerp(px, py, 0.5); // fishhooks forming a bond meet in the middle
    else {
      // a new bond (to another molecule, or within one: a 1,2-shift, a ring closure): point at the atom the
      // electrons are heading for
      const from = ar.from.atoms;
      const n = from.includes(x) && !from.includes(y) ? y : from.includes(y) && !from.includes(x) ? x : y;
      const t = get(n);
      ep.T = pos(n);
      ep.endRadius = Math.max(t.sp.geo.labelRadius(t.a.index) + 3, L * 0.22);
      ep.targetAtom = t;
    }
  }
  if (V.len(V.sub(ep.T, ep.S0)) < 1) ep.T = V.add(ep.T, { x: L * 0.3, y: 0 });
  // text boxes of the atom the arrow points at, so the head can stop just outside them
  ep.targetBoxes = [];
  if (ep.targetAtom) {
    const { sp, a } = ep.targetAtom;
    const g = sp.geo.obstacles();
    for (const l of g.labels) if (l.atom === a.index) ep.targetBoxes.push({ c: { x: l.c.x + sp.geo.dx, y: l.c.y + sp.geo.dy }, hw: l.hw, hh: l.hh });
  }

  // the "away from the molecule" side, used as a mild preference
  let cx = 0; let cy = 0; let cnt = 0;
  for (const sp of ep.involved) for (const a of sp.mol.atoms) { const p = sp.geo.atomPos(a.index); cx += p.x; cy += p.y; cnt++; }
  const centroid = { x: cx / cnt, y: cy / cnt };
  const nrm = V.perp(V.norm(V.sub(ep.T, ep.S0)));
  let side = V.dot(V.sub(V.lerp(ep.S0, ep.T, 0.5), centroid), nrm);
  if (Math.abs(side) < L * 0.15) side = -nrm.y;
  ep.awaySide = side >= 0 ? 1 : -1;
  return ep;
}

// Placement words -> allowed sides/scales.
function placementSides(ep, place) {
  const nrm = V.perp(V.norm(V.sub(ep.T, ep.S0)));
  const dirs = { above: { x: 0, y: -1 }, below: { x: 0, y: 1 }, left: { x: -1, y: 0 }, right: { x: 1, y: 0 } };
  if (!place || !place.side) return [1, -1];
  const d = dirs[place.side];
  return [V.dot(nrm, d) >= 0 ? 1 : -1];
}

// Where the arrow should stop when arriving at T from direction `from`: at least endRadius out, and
// further if needed to clear the target's label text (e.g. the "H2" of a CH2) by 3 px.
function endPoint(ep, from) {
  const u = V.norm(V.sub(ep.T, from));
  let d = ep.endRadius;
  const inside = (p) => ep.targetBoxes.some((b) => Math.abs(p.x - b.c.x) < b.hw + 3 && Math.abs(p.y - b.c.y) < b.hh + 3);
  while (d < ep.endRadius + 40 && inside(V.sub(ep.T, V.mul(u, d)))) d += 1;
  return V.sub(ep.T, V.mul(u, d));
}

// Build one candidate curve. side: +1/-1 (left normal of S0->T), scale: bow size, t: where the peak sits.
function buildCandidate(ep0, o, { side, scale, t, perp, along }) {
  const L = o.bondLength;
  // optionally start part-way along the bond instead of at its midpoint
  const ep = along != null && ep0.bondEnds ? { ...ep0, S0: V.lerp(ep0.bondEnds[0], ep0.bondEnds[1], along) } : ep0;
  let S = ep.S0;
  const T = ep.T;
  let dot = null;
  const shape = (S1) => {
    const chord = V.sub(T, S1);
    const clen = Math.max(V.len(chord), 1);
    const nrm = V.perp(V.norm(chord));
    const bow = Math.max(L * 0.9, Math.min(clen * 0.5, L * 1.6)) * scale;
    return { C: V.add(V.add(S1, V.mul(chord, t)), V.mul(nrm, side * bow)), t, h: (side * bow) / clen, clen };
  };
  let sh = shape(S);
  if (ep.bondNormal) {
    // leave the bond from the side the curve heads to, not from on top of the bond line
    const n = ep.bondNormal;
    const sgn = V.dot(n, V.sub(sh.C, S)) >= 0 ? 1 : -1;
    S = V.add(S, V.mul(n, 4 * sgn));
    sh = shape(S);
  }
  if (ep.dot) {
    // the lone pair faces the curve; the curve then starts just outside it
    const { sp, i } = ep.dot;
    const prefer = V.norm(V.add(V.norm(V.sub(sh.C, S)), V.mul(V.norm(V.sub(T, S)), 0.3)));
    const d = sp.geo.dotFor(i, prefer);
    dot = { ...ep.dot, angle: d.t };
    S = V.add({ x: d.center.x + sp.geo.dx, y: d.center.y + sp.geo.dy }, V.mul(d.dir, 3));
    // keep the absolute bow when re-aiming from the dot
    const chord = V.sub(T, S);
    const clen = Math.max(V.len(chord), 1);
    const nrm = V.perp(V.norm(chord));
    const absBow = sh.h * sh.clen;
    sh = { C: V.add(V.add(S, V.mul(chord, sh.t)), V.mul(nrm, absBow)), t: sh.t, h: absBow / clen, clen };
  }
  // Curves are cubic Beziers S-C1-C2-E. Most candidates are quadratic shapes written as cubics; "perp"
  // candidates leave a bond at right angles and then sweep round to the target (two independent tangents).
  let C = sh.C;
  let C1; let C2; let E;
  if (perp && ep.bondNormal) {
    const n = V.mul(ep.bondNormal, perp.sgn);
    S = V.add(ep.S0, V.mul(n, 4));
    C1 = V.add(S, V.mul(n, perp.d));
    E = endPoint(ep, C1);
    for (let it = 0; it < 3; it++) {
      const dir = V.norm(V.sub(E, C1));
      C2 = V.add(V.lerp(C1, E, 0.55), V.mul(V.perp(dir), perp.side * perp.bow));
      E = endPoint(ep, C2);
    }
    const chord = V.sub(T, S);
    C = { x: (3 * (C1.x + C2.x) - S.x - E.x) / 4, y: (3 * (C1.y + C2.y) - S.y - E.y) / 4 };
    sh = { C, t: V.dot(V.sub(C, S), chord) / Math.max(1, V.dot(chord, chord)), h: V.dot(V.sub(C, S), V.perp(V.norm(chord))) / Math.max(1, V.len(chord)), clen: Math.max(1, V.len(chord)) };
  } else {
    // The curve really ends at E (outside the target label), not at the atom centre T. Re-shape it between
    // S and E with the same side and bow, so short arrows (a bond's electrons moving onto its own atom)
    // still arch, and the head follows the curve instead of an almost-zero end tangent.
    E = endPoint(ep, C);
    const absBow = sh.h * sh.clen;
    for (let it = 0; it < 3; it++) {
      const chord = V.sub(E, S);
      const nrm = V.perp(V.norm(chord));
      C = V.add(V.add(S, V.mul(chord, sh.t)), V.mul(nrm, absBow));
      E = endPoint(ep, C);
    }
    C1 = V.lerp(S, C, 2 / 3);
    C2 = V.lerp(E, C, 2 / 3);
  }
  return finishCandidate(ep, o, { S, C, C1, C2, E, T, t: sh.t, h: sh.h, chordLen: sh.clen, dot, side: perp ? perp.side : side, scale });
}

// Hook: the shape for an arrow whose ends are close together (a bond's electrons moving onto one of its own
// atoms, a lone pair becoming the bond next to it). It leaves its source at right angles, and arrives at
// the target from the side (phi: angle of the arrival point around the target, measured from the source
// side), so the tail and head stay apart and the curve is a single open arch that cannot touch itself.
function buildHook(ep, o, { side, phi, scale }) {
  const T = ep.T;
  const d = V.norm(V.sub(T, ep.S0));
  const nrm = V.mul(V.perp(d), side);
  let S = ep.S0;
  let leave = nrm;
  let dot = null;
  if (ep.bondNormal) {
    const out = V.mul(ep.bondNormal, V.dot(ep.bondNormal, nrm) >= 0 ? 1 : -1);
    S = V.add(ep.S0, V.mul(out, 4));
    // leave the bond outwards but already leaning towards the target, so the curve never swings away first
    leave = V.norm(V.add(out, V.mul(d, 0.7)));
  } else if (ep.dot) {
    const { sp, i } = ep.dot;
    const dd = sp.geo.dotFor(i, V.norm(V.add(nrm, V.mul(d, 0.35))));
    dot = { ...ep.dot, angle: dd.t };
    leave = dd.dir;
    S = V.add({ x: dd.center.x + sp.geo.dx, y: dd.center.y + sp.geo.dy }, V.mul(dd.dir, 3));
  }
  const u = V.norm(V.add(V.mul(d, -Math.cos(phi)), V.mul(nrm, Math.sin(phi))));
  const E = ep.targetAtom ? endPoint(ep, V.add(T, V.mul(u, 60))) : V.add(T, V.mul(u, 5));
  const chord = Math.max(V.len(V.sub(E, S)), 1);
  const k = chord * scale;
  const C1 = V.add(S, V.mul(leave, k));
  const C2 = V.add(E, V.mul(V.norm(V.sub(E, T)), k));
  const C = { x: (3 * (C1.x + C2.x) - S.x - E.x) / 4, y: (3 * (C1.y + C2.y) - S.y - E.y) / 4 };
  return finishCandidate(ep, o, { S, C, C1, C2, E, T, t: 0.5, h: 0, chordLen: chord, dot, side, scale, hook: true });
}

function finishCandidate(ep, o, c) {
  const L = o.bondLength;
  const { S, C1, C2, E, T } = c;
  const pts = [];
  for (let k = 0; k <= N; k++) pts.push(cubic(S, C1, C2, E, k / N));
  let length = 0;
  for (let k = 1; k <= N; k++) length += V.len(V.sub(pts[k], pts[k - 1]));
  // inner part: skip where the curve leaves its source and arrives at its target
  const inner = pts.filter((p) => V.len(V.sub(p, S)) > 7 && V.len(V.sub(p, E)) > 6 && V.len(V.sub(p, T)) > ep.endRadius + 4);
  // degenerate shapes: the head has no clear direction (end on top of its control point), the arrow is a
  // stub, or it is too flat to read as a curved arrow (its middle barely leaves the straight line S-E)
  const se = V.sub(E, S);
  const mid = cubic(S, C1, C2, E, 0.5);
  const sag = Math.abs(V.dot(V.sub(mid, S), V.perp(V.norm(se))));
  // S-shaped (control points on opposite sides of S-E): an electron-pushing arrow is a single arch
  const cross = (q) => se.x * (q.y - S.y) - se.y * (q.x - S.x);
  const sCurve = cross(C1) * cross(C2) < 0 && Math.min(Math.abs(cross(C1)), Math.abs(cross(C2))) / Math.max(1, V.len(se)) > 1.5;
  const degenerate = sCurve || V.len(V.sub(E, C2)) < 6 || length < L * 0.6 || sag < Math.max(5, 0.12 * V.len(se));
  // cramped: far-apart parts of the curve (a third of it or more) come within 6 px, or the head lands
  // next to the tail; snake: the curve's bend changes direction
  let selfGap = V.len(se);
  const far = Math.ceil(pts.length / 3);
  for (let i = 0; i < pts.length; i++) for (let j = i + far; j < pts.length; j++) selfGap = Math.min(selfGap, V.len(V.sub(pts[i], pts[j])));
  const cramped = selfGap < 6;
  let turns = 0; let last = 0;
  for (let k = 2; k < pts.length; k++) {
    const a = V.sub(pts[k - 1], pts[k - 2]); const b = V.sub(pts[k], pts[k - 1]);
    const tt = Math.atan2(a.x * b.y - a.y * b.x, V.dot(a, b));
    if (Math.abs(tt) < 0.035) continue;
    if (last && Math.sign(tt) !== last) turns++;
    last = Math.sign(tt);
  }
  return { ...c, nrm: V.perp(V.norm(V.sub(T, S))), samples: pts, inner, length, degenerate, cramped, snake: turns > 0 };
}

// ---------- obstacles of one step, in page coordinates ----------

export function stepObstacles(state, units, o) {
  const fs = o.fontSize;
  const ob = { bonds: [], labels: [], charges: [], vertices: [], rings: [], furniture: [], furnitureSegs: [] };
  for (const sp of state.species) {
    if (!sp.mol || !sp.geo.obstacles) continue;
    const g = sp.geo.obstacles();
    const off = (p) => ({ x: p.x + sp.geo.dx, y: p.y + sp.geo.dy });
    for (const b of g.bonds) ob.bonds.push({ a: off(b.a), b: off(b.b), atoms: b.atoms, sp });
    for (const l of g.labels) ob.labels.push({ c: off(l.c), r: l.r, hw: l.hw, hh: l.hh, atom: l.atom, sp });
    for (const c of g.charges) ob.charges.push({ c: off(c.c), r: c.r, atom: c.atom, sp });
    for (const v of g.vertices) ob.vertices.push({ c: off(v.c), atom: v.atom, sp });
    for (const r of g.rings || []) ob.rings.push({ poly: r.poly.map(off), maps: r.atoms.map((i) => sp.mol.atoms[i].map), sp });
  }
  for (const u of units) {
    if (u.kind === 'plus') ob.furniture.push({ c: { x: u.x + u.w / 2, y: u.mid }, r: fs * 0.55 });
    if (u.kind === 'transition') {
      ob.furnitureSegs.push({ a: { x: u.x, y: u.mid }, b: { x: u.x + u.w, y: u.mid } });
      if (u.tr?.label) ob.furnitureSegs.push({ a: { x: u.x, y: u.mid - 13 }, b: { x: u.x + u.w, y: u.mid - 13 } });
    }
  }
  return ob;
}

// ---------- costs ----------

// Does the whole curve cross bond b? The bonds the arrow starts from or ends on are exempt; every other bond
// counts, including ones right next to the start (an arrow leaving a bond must not cut its neighbours).
function ownBond(ar, b) {
  const m = [b.sp.mol.atoms[b.atoms[0]].map, b.sp.mol.atoms[b.atoms[1]].map];
  const same = (r) => r.type === 'bond' && m.includes(r.atoms[0]) && m.includes(r.atoms[1]);
  return same(ar.from) || same(ar.to);
}
function crossesBond(samples, b) {
  for (let k = 1; k < samples.length; k++) if (segIntersect(samples[k - 1], samples[k], b.a, b.b)) return true;
  return false;
}

// Closest approach (px) of curve points to anything the arrow is not about.
function clearance(pts, ep, ob) {
  const own = new Set([...ep.ar.from.atoms, ...ep.ar.to.atoms]);
  const mine = (sp, i) => { const m = sp.mol.atoms[i].map; return m != null && own.has(m); };
  let d = Infinity;
  const box = (p, c, hw, hh) => Math.max(0, Math.max(Math.abs(p.x - c.x) - hw, Math.abs(p.y - c.y) - hh));
  for (const p of pts) {
    for (const b of ob.bonds) if (!mine(b.sp, b.atoms[0]) && !mine(b.sp, b.atoms[1])) d = Math.min(d, distPointSeg(p, b.a, b.b));
    for (const l of ob.labels) if (!mine(l.sp, l.atom)) d = Math.min(d, box(p, l.c, l.hw - 1.5, l.hh - 1.5));
    for (const c of ob.charges) if (!mine(c.sp, c.atom)) d = Math.min(d, V.len(V.sub(p, c.c)) - c.r + 1);
    for (const fz of ob.furniture) d = Math.min(d, box(p, { x: fz.c.x, y: fz.c.y - fz.r * 0.09 }, fz.r * 0.6, fz.r * 0.6));
    for (const sg of ob.furnitureSegs) d = Math.min(d, distPointSeg(p, sg.a, sg.b));
  }
  return d;
}

function unaryCost(cand, ep, ob) {
  const pts = cand.inner;
  let cost = cand.length * COST.length + (cand.degenerate ? COST.degenerate : 0) + (cand.cramped ? COST.cramped : 0) + (cand.snake ? COST.snake : 0);
  if (cand.side && cand.side !== ep.awaySide) cost += COST.offSide;
  if (pts.length < 2) return cost;
  for (const b of ob.bonds) if (!ownBond(ep.ar, b) && crossesBond(cand.samples, b)) cost += COST.bondCross;
  // labels: the whole curve counts (an arrow leaving a bond must not cut through the atom at its end);
  // only a lone-pair start sits right next to its own label, so that one gets a little slack
  const away = cand.samples.filter((p) => V.len(V.sub(p, cand.S)) > 3);
  for (const l of ob.labels) {
    const map = l.sp.mol.atoms[l.atom].map;
    const slack = ep.dot && map != null && map === ep.ar.from.atoms[0] ? 2 : -2;
    if (away.some((p) => inLabel(p, l, slack))) cost += COST.labelHit;
  }
  for (const c of ob.charges) if (away.some((p) => V.len(V.sub(p, c.c)) < c.r + 2)) cost += COST.chargeHit;
  for (const v of ob.vertices) if (pts.some((p) => V.len(V.sub(p, v.c)) < 4)) cost += COST.atomGraze;
  for (const r of ob.rings) {
    // fine inside a ring when both ends sit on it (a pi shift around the ring); otherwise stay outside
    const allOn = [...ep.ar.from.atoms, ...ep.ar.to.atoms].every((n) => r.maps.includes(n));
    if (!allOn && cand.samples.slice(1, -1).some((p) => inPoly(p, r.poly))) cost += COST.ringFace;
  }
  for (const fz of ob.furniture) if (pts.some((p) => V.len(V.sub(p, fz.c)) < fz.r)) cost += COST.furniture;
  // clearance: an arrow that skims a label, bond, charge, "+" sign or reaction arrow reads as touching it
  const near = clearance(pts, ep, ob);
  if (near < CLEAR) cost += COST.near * (1 - Math.max(0, near) / CLEAR);
  for (const s of ob.furnitureSegs) {
    for (let k = 1; k < pts.length; k++) if (segIntersect(pts[k - 1], pts[k], s.a, s.b)) { cost += COST.furniture; break; }
  }
  return cost;
}

function pairCost(a, b) {
  let cost = 0;
  const A = a.samples; const B = b.samples;
  let crossings = 0;
  for (let i = 1; i < A.length; i++) {
    for (let j = 1; j < B.length; j++) if (segIntersect(A[i - 1], A[i], B[j - 1], B[j])) crossings++;
  }
  cost += crossings * COST.arrowCross;
  // running alongside each other (ignoring shared endpoints)
  // running alongside each other (ignoring shared endpoints): cost grows with how much of the curve is close
  const ai = a.inner; const bi = b.inner;
  const closeFrac = (P, Q) => {
    if (P.length < 2 || Q.length < 2) return 0;
    let near = 0;
    for (const p of P) for (let j = 1; j < Q.length; j++) if (distPointSeg(p, Q[j - 1], Q[j]) < 11) { near++; break; }
    return near / P.length;
  };
  cost += COST.arrowNear * 4 * Math.max(closeFrac(ai, bi), closeFrac(bi, ai));
  // heads/tails on top of each other
  for (const [p, q] of [[a.E, b.E], [a.E, b.S], [a.S, b.E], [a.S, b.S]]) if (V.len(V.sub(p, q)) < 9) cost += COST.headsTouch;
  // both start from the same lone-pair slot
  if (a.dot && b.dot && a.dot.sp === b.dot.sp && a.dot.i === b.dot.i && Math.abs(Math.atan2(Math.sin(a.dot.angle - b.dot.angle), Math.cos(a.dot.angle - b.dot.angle))) < 0.7) cost += COST.dotClash;
  return cost;
}

// ---------- optimisation ----------

export function placeArrows(state, units, byMap, o) {
  const eps = state.arrows.map((ar) => arrowEndpoints(ar, byMap, o));
  const ob = stepObstacles(state, units, o);
  const cands = eps.map((ep) => {
    if (!ep) return [];
    const ar = ep.ar;
    if (o.optimizeArrows === false) return [buildCandidate(ep, o, { side: ep.awaySide, scale: 1, t: 0.5 })];
    const sides = placementSides(ep, ar.place);
    // o.fastArrows: a coarser search, used when comparing whole-drawing variants (auto-orientation)
    const scales = ar.place?.size === 'tight' ? [0.45, 0.65] : ar.place?.size === 'wide' ? [1.5, 2.0] : o.fastArrows ? [0.85, 1.25, 1.75] : SCALES;
    const peaks = o.fastArrows ? [0.5, 0.32, 0.68] : PEAKS;
    const list = [];
    for (const side of sides) for (const scale of scales) for (const t of peaks) list.push(buildCandidate(ep, o, { side, scale, t }));
    // ends close together: open hooks that arrive at the target from the side
    if (V.len(V.sub(ep.T, ep.S0)) < 1.4 * o.bondLength && ar.place?.size !== 'wide') {
      const phis = o.fastArrows ? [2.0] : [1.6, 2.0, 2.4];
      for (const side of sides) for (const phi of phis) for (const scale of [0.45, 0.65, 0.85]) list.push(buildHook(ep, o, { side, phi, scale }));
    }
    if (o.fastArrows) return list;
    if (ep.bondEnds) {
      for (const along of [0.36, 0.64]) for (const side of sides) for (const scale of [0.85, 1.25]) list.push(buildCandidate(ep, o, { side, scale, t: 0.5, along }));
    }
    if (ep.bondNormal && !ar.place) {
      const L = o.bondLength;
      for (const sgn of [1, -1]) for (const d of [0.6 * L, 1.0 * L]) for (const side of [1, -1]) for (const bow of [0.5 * L, 0.9 * L]) list.push(buildCandidate(ep, o, { perp: { sgn, d, side, bow } }));
    }
    return list;
  });
  const unary = cands.map((list, i) => list.map((c) => unaryCost(c, eps[i], ob)));
  const n = cands.length;
  const pair = new Map();
  const pc = (i, a, j, b) => {
    const key = i < j ? `${i}:${a}|${j}:${b}` : `${j}:${b}|${i}:${a}`;
    if (!pair.has(key)) pair.set(key, pairCost(cands[i][a], cands[j][b]));
    return pair.get(key);
  };
  const total = (choice) => {
    let c = 0;
    for (let i = 0; i < n; i++) if (cands[i].length) c += unary[i][choice[i]];
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (cands[i].length && cands[j].length) c += pc(i, choice[i], j, choice[j]);
    return c;
  };

  let best = cands.map((list, i) => (list.length ? unary[i].indexOf(Math.min(...unary[i])) : -1));
  const combos = cands.reduce((p, l) => p * Math.max(1, l.length), 1);
  if (combos <= 100000) {
    let bestCost = Infinity;
    const choice = new Array(n).fill(0);
    const rec = (i) => {
      if (i === n) { const c = total(choice); if (c < bestCost) { bestCost = c; best = [...choice]; } return; }
      if (!cands[i].length) { choice[i] = -1; rec(i + 1); return; }
      for (let a = 0; a < cands[i].length; a++) { choice[i] = a; rec(i + 1); }
    };
    rec(0);
  } else {
    // iterated conditional modes from the best individual choices
    let improved = true;
    let rounds = 0;
    while (improved && rounds++ < 20) {
      improved = false;
      for (let i = 0; i < n; i++) {
        if (!cands[i].length) continue;
        let bi = best[i]; let bc = total(best);
        for (let a = 0; a < cands[i].length; a++) {
          const trial = [...best]; trial[i] = a;
          const c = total(trial);
          if (c < bc - 1e-9) { bc = c; bi = a; }
        }
        if (bi !== best[i]) { best[i] = bi; improved = true; }
      }
    }
  }
  return eps.map((ep, i) => (ep && best[i] >= 0 ? cands[i][best[i]] : null));
}

// Quality numbers for a rendered document (used by tests, the review sheet and fix-it warnings).
// problems: one entry per arrow that still crosses something, with what it crosses.
export function arrowReport(doc, o) {
  const rep = { arrowCrossings: 0, bondCrossings: 0, labelHits: 0, problems: [] };
  for (const st of doc.states) {
    const issues = new Map(); // arrow -> Set of descriptions
    const note = (ar, what) => { if (!issues.has(ar)) issues.set(ar, new Set()); issues.get(ar).add(what); };
    const withGeo = st.arrows.filter((a) => a.geo);
    for (let i = 0; i < withGeo.length; i++) for (let j = i + 1; j < withGeo.length; j++) {
      const A = withGeo[i].geo.samples; const B = withGeo[j].geo.samples;
      let n = 0;
      for (let x = 1; x < A.length; x++) for (let y = 1; y < B.length; y++) if (segIntersect(A[x - 1], A[x], B[y - 1], B[y])) n++;
      if (n) {
        rep.arrowCrossings += n;
        note(withGeo[i], `crosses the arrow on line ${withGeo[j].line}`);
        note(withGeo[j], `crosses the arrow on line ${withGeo[i].line}`);
      }
    }
    const ob = stepObstacles(st, [], o);
    for (const ar of withGeo) {
      const g = ar.geo;
      const away = g.samples.filter((p) => V.len(V.sub(p, g.S)) > 3);
      for (const b of ob.bonds) {
        if (ownBond(ar, b) || !crossesBond(g.samples, b)) continue;
        rep.bondCrossings++;
        note(ar, 'crosses a bond');
      }
      for (const l of ob.labels) {
        const a = l.sp.mol.atoms[l.atom];
        // the report measures overlap with the text itself (boxes carry 1.5 px of padding for the optimiser)
        const slack = (ar.from.type === 'atom' && a.map === ar.from.atoms[0] ? 2 : 0) + 1.5;
        if (!away.some((p) => inLabel(p, l, slack))) continue;
        rep.labelHits++;
        note(ar, `runs through the ${a.symbol}${a.map != null ? ` (atom ${a.map})` : ''} label`);
      }
    }
    for (const [ar, set] of issues) rep.problems.push({ ar, state: st, what: [...set] });
  }
  return rep;
}

// Candidate fixes for an arrow the renderer could not route cleanly: the arrow line with the other side
// words, and orientation lines for the molecule the arrow starts from. Each is an edit to the source text;
// index.js tries them and only suggests the ones that actually remove the problem.
export function fixCandidates(problem) {
  const { ar, state } = problem;
  const base = ar.text.replace(/\s+(above|below|left|right)\b/, '');
  const out = [];
  for (const side of ['above', 'below', 'left', 'right']) {
    if (ar.place?.side === side) continue;
    out.push({ label: `"${base} ${side}"`, edit: { type: 'replace', line: ar.line, text: `  ${base} ${side}` } });
  }
  const src = ar.from.atoms[0];
  for (const dir of ['right', 'left', 'up', 'down']) out.push({ label: `"face ${src} ${dir}" under the species line`, edit: { type: 'insert', after: state.line, text: `  face ${src} ${dir}` } });
  out.push({ label: `"flip ${src}" under the species line`, edit: { type: 'insert', after: state.line, text: `  flip ${src}` } });
  return out;
}

export function describeProblem(problem) {
  return `arrow ${problem.what.join(' and ')}`;
}
