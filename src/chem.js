// Chemistry layer: turns parsed species into molecule graphs with 2D coordinates (via RDKit.js)
// and runs the checks that need the molecule (bonds exist, atoms have lone pairs, SMILES valid).

const SYMBOLS = ['', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
  'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe'];
const VALENCE_ELECTRONS = { H: 1, Li: 1, Na: 1, K: 1, Be: 2, Mg: 2, B: 3, Al: 3, C: 4, Si: 4, Sn: 4, N: 5, P: 5, As: 5,
  O: 6, S: 6, Se: 6, F: 7, Cl: 7, Br: 7, I: 7 };

export const symbolOf = (z) => SYMBOLS[z] || '?';
export const zOf = (sym) => SYMBOLS.indexOf(sym);

const MOL_OPTS = JSON.stringify({ removeHs: false, setAromaticity: false, kekulize: true });

// Builds a molecule graph from a SMILES string. Returns { mol } or { error }.
export function buildMolecule(RDKit, smiles, atomsFromText) {
  let m = null;
  try { m = RDKit.get_mol(smiles, MOL_OPTS); } catch (e) { m = null; }
  if (!m || !m.is_valid()) {
    if (m) m.delete();
    return { error: `RDKit cannot read SMILES "${smiles}" (check brackets, ring numbers and valences)` };
  }
  try {
    m.set_new_coords(true);
    const json = JSON.parse(m.get_json());
    const d = json.defaults;
    const raw = json.molecules[0];
    const coords = raw.conformers[0].coords;
    const atoms = raw.atoms.map((a, i) => {
      const z = a.z ?? d.atom.z;
      return {
        index: i,
        z,
        symbol: symbolOf(z),
        charge: a.chg ?? d.atom.chg,
        hCount: a.impHs ?? d.atom.impHs,
        radicals: a.nRad ?? d.atom.nRad,
        map: atomsFromText && atomsFromText[i] ? atomsFromText[i].map : null,
        x: coords[i][0],
        y: -coords[i][1], // SVG y points down
        bonds: [],
      };
    });
    const bonds = (raw.bonds || []).map((b, i) => {
      const bond = { index: i, a: b.atoms[0], b: b.atoms[1], order: b.bo ?? d.bond.bo };
      atoms[bond.a].bonds.push(i);
      atoms[bond.b].bonds.push(i);
      return bond;
    });
    const rings = [];
    for (const ext of raw.extensions || []) if (ext.atomRings) rings.push(...ext.atomRings);
    // Rings from RDKit's extension are only present when aromaticity is perceived; find small rings ourselves.
    const mol = { smiles, atoms, bonds, rings: rings.length ? rings : findRings(atoms, bonds) };
    for (const a of atoms) a.lonePairs = lonePairs(mol, a);
    return { mol };
  } finally {
    m.delete();
  }
}

export function bondOrderSum(mol, atom) {
  return atom.bonds.reduce((s, bi) => s + mol.bonds[bi].order, 0);
}

export function lonePairs(mol, atom) {
  const ve = VALENCE_ELECTRONS[atom.symbol];
  if (ve == null) return 0;
  const left = ve - atom.charge - bondOrderSum(mol, atom) - atom.hCount - atom.radicals;
  return Math.max(0, Math.floor(left / 2));
}

export function findBond(mol, i, j) {
  return mol.bonds.find((b) => (b.a === i && b.b === j) || (b.a === j && b.b === i)) || null;
}

// Smallest rings up to size 8, by BFS from each bond.
function findRings(atoms, bonds) {
  const adj = atoms.map(() => []);
  bonds.forEach((b) => { adj[b.a].push(b.b); adj[b.b].push(b.a); });
  const seen = new Set();
  const rings = [];
  for (const b of bonds) {
    // shortest path from b.a to b.b not using the bond itself
    const prev = new Map([[b.a, -1]]);
    const queue = [b.a];
    let found = false;
    while (queue.length && !found) {
      const u = queue.shift();
      for (const v of adj[u]) {
        if (u === b.a && v === b.b) continue;
        if (prev.has(v)) continue;
        prev.set(v, u);
        if (v === b.b) { found = true; break; }
        queue.push(v);
      }
    }
    if (!found) continue;
    const ring = [];
    for (let v = b.b; v !== -1; v = prev.get(v)) ring.push(v);
    if (ring.length > 8) continue;
    const key = [...ring].sort((x, y) => x - y).join(',');
    if (!seen.has(key)) { seen.add(key); rings.push(ring); }
  }
  return rings;
}

// Builds molecules for every species and checks arrows against them.
// Adds to doc.errors / doc.warnings. Returns the same doc with sp.mol set where possible.
// Reaction-aware layout: when the bonds a step forms (present in the next step) or the bonds the previous
// step broke close a ring, lay the molecules out as part of that ring, the way textbooks draw pericyclic and
// intramolecular reactions (s-cis diene for a Diels-Alder, a six-membered ring for a Cope or an
// intramolecular proton transfer). Only changes the shape of molecules inside such a ring.
function mappedBonds(state) {
  const set = new Set();
  for (const sp of state.species) {
    if (!sp.mol) continue;
    for (const b of sp.mol.bonds) {
      const p = sp.mol.atoms[b.a].map; const q = sp.mol.atoms[b.b].map;
      if (p != null && q != null) set.add(p < q ? `${p}-${q}` : `${q}-${p}`);
    }
  }
  return set;
}

function reactionLayout(RDKit, doc) {
  const bondSets = doc.states.map(mappedBonds);
  doc.states.forEach((state, si) => {
    const nodes = []; // {sp, i}
    const nodeOf = new Map(); // map number -> node index
    for (const sp of state.species) {
      if (!sp.mol) continue;
      sp._node0 = nodes.length;
      sp.mol.atoms.forEach((a, i) => { if (a.map != null) nodeOf.set(a.map, nodes.length); nodes.push({ sp, i }); });
    }
    if (!nodes.length) return;
    const extra = new Set();
    const forming = new Set(); // bonds the next step forms (as opposed to ones the previous step broke)
    const addDiff = (other, next) => {
      for (const key of other) {
        if (bondSets[si].has(key)) continue;
        const [p, q] = key.split('-').map(Number);
        if (nodeOf.has(p) && nodeOf.has(q)) { extra.add(key); if (next) forming.add(key); }
      }
    };
    if (si + 1 < doc.states.length) addDiff(bondSets[si + 1], true);
    if (si > 0) addDiff(bondSets[si - 1], false);
    if (!extra.size) return;
    // union-find over the step's own bonds, then see which extra bonds close a ring
    const parent = nodes.map((_, i) => i);
    const find = (x) => (parent[x] === x ? x : (parent[x] = find(parent[x])));
    const edges = [];
    for (const sp of state.species) {
      if (!sp.mol) continue;
      for (const b of sp.mol.bonds) {
        const u = sp._node0 + b.a; const v = sp._node0 + b.b;
        edges.push([u, v]);
        parent[find(u)] = find(v);
      }
    }
    // Only rings of 5 to 7 atoms get this treatment (six-membered pericyclic transition states, five- and
    // six-membered ring closures), plus 4-membered ones that two molecules are about to form together (the
    // four-centre transition state of a hydroboration). Laying a molecule out around a small ring that
    // exists only for an instant (a 1,3-proton shift, a bond that just broke out of an epoxide or bromonium
    // ion, a Wittig betaine closing) bends its real bonds out of shape.
    const keys = [...extra];
    const allExtra = keys.map((k) => k.split('-').map((n) => nodeOf.get(Number(n))));
    const adj = nodes.map(() => []);
    for (const [u, v] of [...edges, ...allExtra]) { adj[u].push(v); adj[v].push(u); }
    const ringSize = (u, v) => {
      const dist = new Map([[u, 0]]);
      const queue = [u];
      while (queue.length) {
        const x = queue.shift();
        for (const y of adj[x]) {
          if ((x === u && y === v) || (x === v && y === u) || dist.has(y)) continue;
          dist.set(y, dist.get(x) + 1);
          if (y === v) return dist.get(y) + 1;
          queue.push(y);
        }
      }
      return Infinity;
    };
    const extraEdges = allExtra.filter(([u, v], k) => {
      const n = ringSize(u, v);
      if (n >= 5 && n <= 7) return true;
      return n === 4 && forming.has(keys[k]) && nodes[u].sp !== nodes[v].sp;
    });
    const closing = [];
    for (const [u, v] of extraEdges) {
      if (find(u) === find(v)) closing.push([u, v]);
      else parent[find(u)] = find(v);
    }
    if (!closing.length) return;
    const ringRoots = new Set(closing.map(([u]) => find(u)));
    // lay out each ring-containing component as a skeleton of dummy atoms
    for (const root of ringRoots) {
      const members = nodes.map((_, i) => i).filter((i) => find(i) === root);
      const local = new Map(members.map((n, k) => [n, k]));
      const g = { atoms: members.map((n, k) => ({ map: k + 1 })), bonds: [] };
      for (const [u, v] of [...edges, ...extraEdges]) if (local.has(u) && local.has(v)) g.bonds.push({ a: local.get(u), b: local.get(v) });
      const smiles = graphToSmiles({ atoms: g.atoms.map((x) => ({ symbol: '*', charge: 0, hCount: 0, map: x.map })), bonds: g.bonds.map((x) => ({ ...x, order: 1 })) }).join('.');
      let m = null;
      try { m = RDKit.get_mol(smiles, JSON.stringify({ removeHs: false })); } catch (e) { m = null; }
      if (!m) continue;
      try {
        m.set_new_coords(true);
        const coords = JSON.parse(m.get_json()).molecules[0].conformers[0].coords;
        // dummy atoms come back in SMILES order; their map numbers say which member they are
        const order = [...smiles.matchAll(/\[\*:(\d+)\]/g)].map((x) => Number(x[1]) - 1);
        const placed = order.map((k, idx) => ({ node: nodes[members[k]], x: coords[idx][0], y: -coords[idx][1] }));
        // several molecules in one ring: turn it so they sit left to right in the order they are written
        const spOrder = [...new Set(placed.map((p) => p.node.sp))].sort((a, b) => state.species.indexOf(a) - state.species.indexOf(b));
        if (spOrder.length > 1) {
          const cen = (sp) => { const ps = placed.filter((p) => p.node.sp === sp); return { x: ps.reduce((t, p) => t + p.x, 0) / ps.length, y: ps.reduce((t, p) => t + p.y, 0) / ps.length }; };
          const a = cen(spOrder[0]);
          const rest = placed.filter((p) => p.node.sp !== spOrder[0]);
          const b = { x: rest.reduce((t, p) => t + p.x, 0) / rest.length, y: rest.reduce((t, p) => t + p.y, 0) / rest.length };
          const th = -Math.atan2(b.y - a.y, b.x - a.x);
          const c = Math.cos(th); const sn = Math.sin(th);
          for (const p of placed) { const x = p.x; p.x = c * x - sn * p.y; p.y = sn * x + c * p.y; }
        }
        for (const p of placed) { p.node.sp.mol.atoms[p.node.i].x = p.x; p.node.sp.mol.atoms[p.node.i].y = p.y; p.node.sp._ringLayout = true; }
      } finally { m.delete(); }
    }
  });
}

// ---------- SMILES writer (bracket atoms, Kekule bonds, map numbers kept) ----------
// Used to write the skeleton graphs of the reaction-aware layout.
// Graph: { atoms: [{symbol, charge, hCount, map}], bonds: [{a, b, order}] }
function atomText(a) {
  const h = a.hCount > 0 ? 'H' + (a.hCount > 1 ? a.hCount : '') : '';
  const c = a.charge === 0 ? '' : (a.charge > 0 ? '+' : '-') + (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '');
  const m = a.map != null ? ':' + a.map : '';
  return `[${a.symbol}${h}${c}${m}]`;
}

const BOND = { 1: '', 2: '=', 3: '#' };

// Returns one SMILES string per connected fragment.
export function graphToSmiles(g) {
  const n = g.atoms.length;
  const adj = g.atoms.map(() => []);
  g.bonds.forEach((b, bi) => { adj[b.a].push(bi); adj[b.b].push(bi); });
  const other = (bi, u) => (g.bonds[bi].a === u ? g.bonds[bi].b : g.bonds[bi].a);
  const visited = new Array(n).fill(false);
  const out = [];
  for (let s = 0; s < n; s++) {
    if (visited[s]) continue;
    const children = new Map();
    const closures = new Map();
    const seenBond = new Set();
    const dfs = (u, parentBond) => {
      visited[u] = true;
      children.set(u, []);
      if (!closures.has(u)) closures.set(u, []);
      for (const bi of adj[u]) {
        if (bi === parentBond || seenBond.has(bi)) continue;
        const v = other(bi, u);
        seenBond.add(bi);
        if (!visited[v]) { children.get(u).push({ v, bi }); dfs(v, bi); }
        else {
          closures.get(u).push(bi);
          if (!closures.has(v)) closures.set(v, []);
          closures.get(v).push(bi);
        }
      }
    };
    dfs(s, -1);
    const digitOf = new Map();
    const free = [];
    let next = 1;
    const emit = (u) => {
      let str = atomText(g.atoms[u]);
      for (const bi of closures.get(u)) {
        if (digitOf.has(bi)) {
          const d = digitOf.get(bi);
          str += d > 9 ? `%${d}` : d;
          digitOf.delete(bi);
          free.push(d);
          free.sort((p, q) => p - q);
        } else {
          const d = free.length ? free.shift() : next++;
          digitOf.set(bi, d);
          str += BOND[g.bonds[bi].order] + (d > 9 ? `%${d}` : d);
        }
      }
      const ch = children.get(u);
      ch.forEach(({ v, bi }, idx) => {
        const part = BOND[g.bonds[bi].order] + emit(v);
        str += idx < ch.length - 1 ? `(${part})` : part;
      });
      return str;
    };
    out.push(emit(s));
  }
  return out;
}

export function analyse(RDKit, doc) {
  for (const [si, state] of doc.states.entries()) {
    const byMap = new Map(); // map -> {sp, atom}
    for (const sp of state.species) {
      if (!sp.smiles) continue;
      const res = buildMolecule(RDKit, sp.smiles, sp.atoms);
      if (res.error) { doc.errors.push({ line: sp.line, message: res.error }); continue; }
      sp.mol = res.mol;
      if (sp.atoms.length !== res.mol.atoms.length) {
        // our tokenizer and RDKit disagree; fall back to no map numbers rather than wrong ones
        doc.warnings.push({ line: sp.line, message: `could not match atom numbers in "${sp.smiles}"` });
        res.mol.atoms.forEach((a) => (a.map = null));
      }
      for (const a of res.mol.atoms) if (a.map != null) byMap.set(a.map, { sp, atom: a });
    }
    for (const ar of state.arrows) {
      const missing = [...ar.from.atoms, ...ar.to.atoms].some((n) => !byMap.has(n));
      if (missing) continue; // parser already reported it
      const get = (n) => byMap.get(n);
      if (ar.from.type === 'bond') {
        const [p, q] = ar.from.atoms.map(get);
        if (p.sp !== q.sp || !findBond(p.sp.mol, p.atom.index, q.atom.index)) {
          doc.errors.push({ line: ar.line, message: `there is no bond between atoms ${ar.from.atoms[0]} and ${ar.from.atoms[1]} in this step, so electrons cannot move out of it` });
        }
      } else if (ar.from.source === 'lp') {
        const a = get(ar.from.atoms[0]).atom;
        if (a.lonePairs === 0) {
          doc.warnings.push({ line: ar.line, message: `atom ${ar.from.atoms[0]} (${a.symbol}${a.charge ? (a.charge > 0 ? '+' : '−') : ''}) has no lone pair to start an arrow from` });
        }
      } else if (ar.from.source === 'rad') {
        const a = get(ar.from.atoms[0]).atom;
        if (a.radicals === 0) {
          doc.warnings.push({ line: ar.line, message: `atom ${ar.from.atoms[0]} has no unpaired electron; write it as a radical in SMILES, e.g. [CH3:${ar.from.atoms[0]}]` });
        }
      }
      if (ar.to.type === 'bond' && ar.from.type === 'bond' &&
          ar.from.atoms.every((n) => ar.to.atoms.includes(n))) {
        doc.warnings.push({ line: ar.line, message: 'arrow starts and ends on the same bond' });
      }
    }
    void si;
  }
  reactionLayout(RDKit, doc);
  stepLayout(RDKit, doc);
  return doc;
}

// ---------- step-to-step layout ----------
// Each step's molecules come from RDKit fresh, so groups can trade places between steps (the OH and OMe of a
// tetrahedral intermediate swapping sides) even after the renderer turns the molecule to match. Instead, lay
// each molecule out with the atoms it shares with the previous step pinned where they were, and let RDKit
// place only what is new. Atoms are matched by map number, then unnumbered neighbours are matched outwards
// from them. Atoms that both lose and gain a bond (a migrating methyl, a transferred H) are not pinned.
// Molecules shaped by the ring-closure layout keep that shape.
function stepLayout(RDKit, doc) {
  for (let si = 1; si < doc.states.length; si++) {
    const prevSpecies = doc.states[si - 1].species.filter((sp) => sp.mol);
    for (const sp of doc.states[si].species) {
      if (!sp.mol || sp._ringLayout || sp.mol.atoms.length < 3) continue;
      // a molecule about to bond to itself (the next step closes a ring) needs room to curl up: leave it free
      if (closesOnItself(sp, doc.states[si + 1])) continue;
      let best = null;
      for (const pp of prevSpecies) {
        const pairs = pairAtoms(sp.mol, pp.mol);
        if (!best || pairs.size > best.pairs.size) best = { pp, pairs };
      }
      if (!best || best.pairs.size < 2) continue;
      const pinned = pinnedAtoms(sp.mol, best.pp.mol, best.pairs);
      if (pinned.length < 2 || pinned.length === sp.mol.atoms.length && samePositions(sp.mol, best.pp.mol, best.pairs)) continue;
      alignToPrevious(RDKit, sp, best.pp.mol, best.pairs, pinned);
    }
  }
}

function closesOnItself(sp, next) {
  if (!next) return false;
  const maps = new Set(sp.mol.atoms.map((a) => a.map).filter((m) => m != null));
  const bonded = mappedBonds({ species: [sp] });
  for (const key of mappedBonds(next)) {
    const [p, q] = key.split('-').map(Number);
    if (maps.has(p) && maps.has(q) && !bonded.has(key)) return true;
  }
  return false;
}

// cur atom index -> prev atom index
function pairAtoms(cur, prev) {
  const pairs = new Map();
  const prevByMap = new Map(prev.atoms.filter((a) => a.map != null).map((a) => [a.map, a.index]));
  for (const a of cur.atoms) if (a.map != null && prevByMap.has(a.map)) pairs.set(a.index, prevByMap.get(a.map));
  const nbrs = (mol, i) => mol.atoms[i].bonds.map((bi) => (mol.bonds[bi].a === i ? mol.bonds[bi].b : mol.bonds[bi].a));
  const used = new Set(pairs.values());
  const queue = [...pairs.keys()];
  while (queue.length) {
    const i = queue.shift(); const j = pairs.get(i);
    const freeCur = nbrs(cur, i).filter((x) => !pairs.has(x) && cur.atoms[x].map == null);
    const freePrev = nbrs(prev, j).filter((y) => !used.has(y) && prev.atoms[y].map == null);
    for (const x of freeCur) {
      // same element and same number of bonds; equivalent groups (two methyls) can pair either way
      const y = freePrev.find((q) => !used.has(q) && prev.atoms[q].symbol === cur.atoms[x].symbol && prev.atoms[q].bonds.length === cur.atoms[x].bonds.length);
      if (y == null) continue;
      pairs.set(x, y); used.add(y); queue.push(x);
    }
  }
  return pairs;
}

function pinnedAtoms(cur, prev, pairs) {
  const nbrs = (mol, i) => mol.atoms[i].bonds.map((bi) => (mol.bonds[bi].a === i ? mol.bonds[bi].b : mol.bonds[bi].a));
  // An atom that lost a bond (a leaving group going, a tetrahedral intermediate collapsing) needs its
  // remaining bonds spread out again, so it and its neighbours are left free. An atom that only gained a
  // bond keeps its old bonds where they were and takes the new one in the gap.
  const free = new Set();
  for (const [i, j] of pairs) {
    const now = nbrs(cur, i).map((x) => (pairs.has(x) ? pairs.get(x) : -1 - x));
    const before = nbrs(prev, j);
    const gained = now.some((y) => !before.includes(y));
    const lost = before.some((y) => !now.includes(y));
    if (lost) { free.add(i); for (const x of nbrs(cur, i)) free.add(x); } else if (gained && lost) free.add(i);
  }
  const pinned = [...pairs.keys()].filter((i) => !free.has(i));
  // a bond between two pinned atoms that did not exist before would be drawn at whatever length the old
  // positions give: leave the molecule to RDKit
  for (const b of cur.bonds) {
    if (!pinned.includes(b.a) || !pinned.includes(b.b)) continue;
    const pa = pairs.get(b.a); const pb = pairs.get(b.b);
    if (!findBond(prev, pa, pb)) return [];
  }
  return pinned;
}

function samePositions(cur, prev, pairs) {
  for (const [i, j] of pairs) if (Math.hypot(cur.atoms[i].x - prev.atoms[j].x, cur.atoms[i].y - prev.atoms[j].y) > 1e-6) return false;
  return true;
}

function isoLines(labels) {
  const out = [];
  for (let k = 0; k < labels.length; k += 8) {
    const chunk = labels.slice(k, k + 8);
    out.push(`M  ISO${String(chunk.length).padStart(3)}${chunk.map(([i, v]) => `${String(i).padStart(4)}${String(v).padStart(4)}`).join('')}`);
  }
  return out;
}

function alignToPrevious(RDKit, sp, prev, pairs, pinned) {
  const mol = sp.mol;
  const label = (i) => 100 + i; // isotope labels tie each template atom to one atom of the molecule
  const plens = prev.bonds.map((b) => Math.hypot(prev.atoms[b.a].x - prev.atoms[b.b].x, prev.atoms[b.a].y - prev.atoms[b.b].y)).sort((x, y) => x - y);
  const L = plens.length ? plens[Math.floor(plens.length / 2)] : 1;
  let m = null;
  try {
    m = RDKit.get_mol(sp.smiles, MOL_OPTS);
    if (!m || !m.is_valid() || m.get_num_atoms() !== mol.atoms.length) return;
    const block = m.get_molblock().replace(/^M {2}END$/m, [...isoLines(mol.atoms.map((a) => [a.index + 1, label(a.index)])), 'M  END'].join('\n'));
    m.delete();
    m = RDKit.get_mol(block, JSON.stringify({ removeHs: false, kekulize: false, setAromaticity: false }));
    if (!m || !m.is_valid()) return;
    // RDKit's own depictor keeps template atoms exactly in place (it works at a bond length of 1.5);
    // CoordGen (bond length 1) only follows the template loosely but copes better when the pinned atoms
    // leave awkward room. Try the exact one first and keep the first layout that passes the checks.
    for (const { coordGen, scale } of [{ coordGen: false, scale: 1.5 }, { coordGen: true, scale: 1 }]) {
      const coords = pinnedLayout(RDKit, m, mol, prev, pairs, pinned, label, scale / L, coordGen);
      if (coords && layoutOk(mol, coords, L)) {
        mol.atoms.forEach((a, i) => { a.x = coords[i][0]; a.y = -coords[i][1]; });
        return;
      }
    }
  } catch (e) {
    // keep the molecule's own layout
  } finally {
    if (m) m.delete();
  }
}

// Coordinates (molecule's own scale, RDKit frame) with the pinned atoms at their previous positions, or null.
function pinnedLayout(RDKit, m, mol, prev, pairs, pinned, label, K, coordGen) {
  const idx = new Map(pinned.map((i, k) => [i, k + 1]));
  const atomLines = pinned.map((i) => {
    const p = prev.atoms[pairs.get(i)];
    return `${(p.x * K).toFixed(4).padStart(10)}${(-p.y * K).toFixed(4).padStart(10)}${'0.0000'.padStart(10)} ${mol.atoms[i].symbol.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`;
  });
  const bondLines = mol.bonds.filter((b) => idx.has(b.a) && idx.has(b.b)).map((b) => `${String(idx.get(b.a)).padStart(3)}${String(idx.get(b.b)).padStart(3)}${String(b.order).padStart(3)}  0`);
  const tblock = ['template', '     MechScript', '', `${String(atomLines.length).padStart(3)}${String(bondLines.length).padStart(3)}  0  0  0  0  0  0  0  0999 V2000`,
    ...atomLines, ...bondLines, ...isoLines(pinned.map((i) => [idx.get(i), label(i)])), 'M  END'].join('\n');
  let tpl = null;
  try {
    tpl = RDKit.get_mol(tblock, JSON.stringify({ removeHs: false, sanitize: false }));
    if (!tpl || !tpl.is_valid()) return null;
    const res = m.generate_aligned_coords(tpl, JSON.stringify({ acceptFailure: false, allowRGroups: false, useCoordGen: coordGen }));
    if (!res) return null;
    return JSON.parse(m.get_json()).molecules[0].conformers[0].coords.map(([x, y]) => [x / K, y / K]);
  } finally {
    if (tpl) tpl.delete();
  }
}

// A pinned layout is only used when it is a clean drawing: bonds of normal length, no atoms on top of each
// other, and no crowded bond angles (outside small rings).
function layoutOk(mol, coords, L) {
  const d = (i, j) => Math.hypot(coords[i][0] - coords[j][0], coords[i][1] - coords[j][1]);
  for (const b of mol.bonds) if (d(b.a, b.b) < 0.75 * L || d(b.a, b.b) > 1.3 * L) return false;
  for (let i = 0; i < coords.length; i++) for (let j = i + 1; j < coords.length; j++) if (!findBond(mol, i, j) && d(i, j) < 0.5 * L) return false;
  const small = new Set((mol.rings || []).filter((r) => r.length <= 4).flat());
  for (const a of mol.atoms) {
    if (small.has(a.index) || a.bonds.length < 2) continue;
    const nb = a.bonds.map((bi) => (mol.bonds[bi].a === a.index ? mol.bonds[bi].b : mol.bonds[bi].a));
    const dirs = nb.map((j) => Math.atan2(coords[j][1] - coords[a.index][1], coords[j][0] - coords[a.index][0])).sort((x, y) => x - y);
    const gaps = dirs.map((t, k) => ((k + 1 < dirs.length ? dirs[k + 1] : dirs[0] + 2 * Math.PI) - t) * 180 / Math.PI);
    const minGap = Math.min(...gaps);
    if (nb.length === 2 ? minGap < 100 : minGap < 52) return false;
    if (nb.length === 3 && Math.max(...gaps) > 170) return false; // a T shape: one bond squeezed to the side
  }
  return true;
}

