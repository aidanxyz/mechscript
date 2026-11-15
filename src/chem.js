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
    const addDiff = (other) => {
      for (const key of other) {
        if (bondSets[si].has(key)) continue;
        const [p, q] = key.split('-').map(Number);
        if (nodeOf.has(p) && nodeOf.has(q)) extra.add(key);
      }
    };
    if (si + 1 < doc.states.length) addDiff(bondSets[si + 1]);
    if (si > 0) addDiff(bondSets[si - 1]);
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
    const extraEdges = [...extra].map((k) => k.split('-').map((n) => nodeOf.get(Number(n))));
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
        for (const p of placed) { p.node.sp.mol.atoms[p.node.i].x = p.x; p.node.sp.mol.atoms[p.node.i].y = p.y; }
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
  return doc;
}
