// MechScript parser: text -> data structure. No dependencies, runs in Node and the browser.
//
// Document shape:
// {
//   title: string|null, titleLine: number|null,
//   states: [{ line, layout: [{ type: 'face'|'flip', atom, dir?, line }], species: [{ smiles, line, start, end, atoms: [{index, symbol, map}] }],
//              arrows: [{ line, kind: 'full'|'fishhook', from: Ref, to: Ref, place: {side?, size?}|null, text }] }],
//   transitions: [{ line, type: 'forward'|'equilibrium'|'resonance', label: string|null, rowBreak?: true }],
//   errors: [{ line, message }], warnings: [{ line, message }]
// }
// Ref: { type: 'atom', atoms: [a], source?: 'lp'|'rad' } or { type: 'bond', atoms: [a, b] }
//
// states[i] and states[i+1] are joined by transitions[i].

const TRANSITION_RE = /^(-->|<=>|<->)\s*(?:"((?:[^"\\]|\\.)*)")?\s*$/;
const ARROW_RE = /^(lp|rad|bond)\s+(\d+(?:\s*-\s*\d+)?)\s*(->|~>)\s*(\d+(?:\s*-\s*\d+)?)(.*)$/;
const SIDE_WORDS = ['above', 'below', 'left', 'right'];
const FACE_RE = /^face\s+(\d+)\s+(right|left|up|down)$/;
const FLIP_RE = /^flip\s+(\d+)$/;
const SIZE_WORDS = ['tight', 'wide'];
const TITLE_RE = /^mechanism\s+"((?:[^"\\]|\\.)*)"\s*$/;

const ORGANIC = ['Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I', 'b', 'c', 'n', 'o', 'p', 's'];

// Tokenises the atoms of a SMILES string in the order RDKit numbers them.
// Returns { atoms: [{index, symbol, map}], error: string|null }.
export function smilesAtoms(smiles) {
  const atoms = [];
  let i = 0;
  while (i < smiles.length) {
    const ch = smiles[i];
    if (ch === '[') {
      const close = smiles.indexOf(']', i);
      if (close < 0) return { atoms, error: `unclosed "[" in "${smiles}"` };
      const body = smiles.slice(i + 1, close);
      const m = body.match(/^(\d*)(Cl|Br|[A-Z][a-z]?|[a-z][a-z]?|\*)/);
      if (!m) return { atoms, error: `cannot read atom "[${body}]"` };
      const mapM = body.match(/:(\d+)$/);
      atoms.push({ index: atoms.length, symbol: m[2], map: mapM ? Number(mapM[1]) : null });
      i = close + 1;
      continue;
    }
    const org = ORGANIC.find((s) => smiles.startsWith(s, i));
    if (org) {
      atoms.push({ index: atoms.length, symbol: org, map: null });
      i += org.length;
      continue;
    }
    if (ch === '%') { i += 3; continue; }
    if ('()=#-:/\\.$~+0123456789@'.includes(ch)) { i++; continue; }
    if (ch === '*') { atoms.push({ index: atoms.length, symbol: '*', map: null }); i++; continue; }
    return { atoms, error: `unexpected character "${ch}" in SMILES "${smiles}"` };
  }
  return { atoms, error: null };
}

function parseRef(text, sourceKind) {
  const parts = text.split('-').map((s) => Number(s.trim()));
  if (parts.length === 2) return { type: 'bond', atoms: parts };
  const ref = { type: 'atom', atoms: parts };
  if (sourceKind) ref.source = sourceKind;
  return ref;
}

// Splits a species line on " + " separators (a "+" with whitespace on both sides;
// charges inside brackets like [OH-] or [NH4+] never have spaces around them).
function splitSpecies(content, offset) {
  const out = [];
  const re = /\s\+\s/g;
  let last = 0;
  let m;
  const push = (s, e) => {
    const raw = content.slice(s, e);
    const lead = raw.length - raw.trimStart().length;
    const smiles = raw.trim();
    if (smiles) out.push({ smiles, start: offset + s + lead, end: offset + s + lead + smiles.length });
    else out.push({ smiles: '', start: offset + s, end: offset + e });
  };
  while ((m = re.exec(content))) {
    push(last, m.index);
    last = m.index + m[0].length;
    re.lastIndex = last - 1;
  }
  push(last, content.length);
  return out;
}

export function parse(text) {
  const doc = { title: null, titleLine: null, states: [], transitions: [], errors: [], warnings: [] };
  const lines = text.split('\n');
  let current = null; // current state
  let pendingTransition = null;
  let pendingBreak = false; // a "---" line: start a new row at the next step
  let offset = 0;

  const err = (line, message) => doc.errors.push({ line, message });
  const warn = (line, message) => doc.warnings.push({ line, message });

  lines.forEach((rawLine, idx) => {
    const lineNo = idx + 1;
    const lineStart = offset;
    offset += rawLine.length + 1;
    // "#" starts a comment only at line start or after whitespace ("#" inside SMILES is a triple bond)
    const cm = rawLine.match(/(^|\s)#/);
    const commentAt = cm ? cm.index + cm[1].length : -1;
    const content = commentAt >= 0 ? rawLine.slice(0, commentAt) : rawLine;
    const line = content.trim();
    if (!line) return;

    let m;
    if ((m = line.match(TITLE_RE))) {
      if (doc.title !== null) err(lineNo, 'only one "mechanism" title line is allowed');
      else if (doc.states.length) err(lineNo, 'the "mechanism" title line must come first');
      doc.title = m[1].replace(/\\"/g, '"');
      doc.titleLine = lineNo;
      return;
    }
    if (/^mechanism\b/.test(line)) {
      err(lineNo, 'title line must look like: mechanism "Title in quotes"');
      return;
    }
    if (line === '---') {
      if (!doc.states.length) { err(lineNo, '"---" starts a new row before a step; put it between steps'); return; }
      pendingBreak = true;
      return;
    }
    if (/^(face|flip)\b/.test(line)) {
      if (!current) { err(lineNo, `"${line.split(/\s+/)[0]}" must come after a line of species (it orients a molecule in that step)`); return; }
      if ((m = line.match(FACE_RE))) current.layout.push({ type: 'face', atom: Number(m[1]), dir: m[2], line: lineNo });
      else if ((m = line.match(FLIP_RE))) current.layout.push({ type: 'flip', atom: Number(m[1]), line: lineNo });
      else err(lineNo, 'orientation lines are "face 3 right" (right, left, up or down) and "flip 3"; the number is an atom of the molecule to turn');
      return;
    }
    if ((m = line.match(TRANSITION_RE))) {
      if (!current) { err(lineNo, 'a reaction arrow needs species on a line before it'); return; }
      if (pendingTransition) { err(lineNo, 'two reaction arrows in a row; put the species of the intermediate step between them'); return; }
      const type = m[1] === '-->' ? 'forward' : m[1] === '<=>' ? 'equilibrium' : 'resonance';
      pendingTransition = { line: lineNo, type, label: m[2] != null ? m[2].replace(/\\"/g, '"') : null };
      if (pendingBreak) { pendingTransition.rowBreak = true; pendingBreak = false; }
      current = null;
      return;
    }
    if (/^(-->|<=>|<->|=>|->)/.test(line)) {
      err(lineNo, 'reaction arrows are "-->", "<=>" or "<->", optionally followed by a "label" in double quotes');
      return;
    }
    if (/^(lp|rad|bond)\b/.test(line)) {
      if (!current) { err(lineNo, 'a curved arrow must come after a line of species (it belongs to the species above it)'); return; }
      if (!(m = line.match(ARROW_RE))) {
        err(lineNo, 'cannot read curved arrow. Forms: "lp 1 -> 2", "lp 1 -> 1-2", "bond 2-3 -> 3", "bond 1-2 -> 2-3", "rad 4 ~> 4-5"; "~>" for a one-electron fishhook');
        return;
      }
      const kindWord = m[1];
      const from = parseRef(m[2], kindWord === 'bond' ? null : kindWord);
      const to = parseRef(m[4]);
      const kind = m[3] === '~>' ? 'fishhook' : 'full';
      if (kindWord === 'bond' && from.type !== 'bond') { err(lineNo, `"bond" needs two atoms joined by "-", e.g. "bond 2-3 -> 3"`); return; }
      if (kindWord !== 'bond' && from.type !== 'atom') { err(lineNo, `"${kindWord}" starts from one atom, e.g. "${kindWord} 1 -> 2"`); return; }
      if (kindWord === 'rad' && kind !== 'fishhook') warn(lineNo, 'a "rad" arrow moves one electron; you probably want "~>"');
      if (from.type === 'bond' && from.atoms[0] === from.atoms[1]) { err(lineNo, `bond ${from.atoms[0]}-${from.atoms[1]} joins an atom to itself`); return; }
      if (to.type === 'bond' && to.atoms[0] === to.atoms[1]) { err(lineNo, `bond ${to.atoms[0]}-${to.atoms[1]} joins an atom to itself`); return; }
      // optional placement: a side word and/or a size word
      const place = {};
      for (const word of m[5].trim().split(/\s+/).filter(Boolean)) {
        if (SIDE_WORDS.includes(word)) {
          if (place.side) { err(lineNo, `arrow has two sides ("${place.side}" and "${word}"); use one of above, below, left, right`); return; }
          place.side = word;
        } else if (SIZE_WORDS.includes(word)) place.size = word;
        else { err(lineNo, `unknown arrow option "${word}"; options are above, below, left, right, tight, wide`); return; }
      }
      current.arrows.push({ line: lineNo, kind, from, to, place: place.side || place.size ? place : null, text: line });
      return;
    }

    // Otherwise: a line of species.
    if (current) {
      err(lineNo, 'two species lines in a row; put all species of one step on one line, joined by " + ", or add a reaction arrow "-->" between steps');
      return;
    }
    const contentStart = lineStart + (content.length - content.trimStart().length);
    const species = splitSpecies(content.trim(), contentStart).map((s) => ({ ...s, line: lineNo, atoms: [] }));
    for (const sp of species) {
      if (!sp.smiles) { err(lineNo, 'empty species: check for a stray "+"'); continue; }
      if (/\s/.test(sp.smiles)) { err(lineNo, `"${sp.smiles}" contains a space; separate species with " + " and arrows go on their own lines`); continue; }
      const { atoms, error } = smilesAtoms(sp.smiles);
      if (error) err(lineNo, error);
      sp.atoms = atoms;
    }
    current = { line: lineNo, species, arrows: [], layout: [] };
    if (pendingTransition) { doc.transitions.push(pendingTransition); pendingTransition = null; }
    doc.states.push(current);
  });

  if (pendingTransition) err(pendingTransition.line, 'the reaction arrow has nothing after it; add the products on the next line');
  if (!doc.states.length && !doc.errors.length) err(1, 'no species found; start with a line of SMILES such as "[OH-:1] + [CH3:2][Br:3]"');

  checkMaps(doc);
  return doc;
}

// Checks that only depend on the text (map numbers), not on chemistry.
function checkMaps(doc) {
  const elementOf = new Map(); // map number -> {symbol, line}
  doc.states.forEach((state, si) => {
    const inState = new Map();
    for (const sp of state.species) {
      for (const a of sp.atoms) {
        if (a.map == null) continue;
        if (inState.has(a.map)) {
          doc.errors.push({ line: sp.line, message: `atom map number ${a.map} is used twice in this step; each atom needs its own number` });
          continue;
        }
        inState.set(a.map, a);
        const sym = a.symbol.charAt(0).toUpperCase() + a.symbol.slice(1);
        const prev = elementOf.get(a.map);
        if (prev && prev.symbol !== sym) {
          doc.errors.push({ line: sp.line, message: `atom ${a.map} is ${sym} here but was ${prev.symbol} on line ${prev.line}; a map number must stay on the same atom through the whole mechanism` });
        } else if (!prev) elementOf.set(a.map, { symbol: sym, line: sp.line });
      }
    }
    for (const lay of state.layout || []) {
      if (!inState.has(lay.atom)) doc.errors.push({ line: lay.line, message: `"${lay.type}" names atom ${lay.atom}, which does not exist in this step; use the map number of an atom in the molecule to turn` });
    }
    for (const ar of state.arrows) {
      const refs = [...ar.from.atoms, ...ar.to.atoms];
      for (const n of refs) {
        if (!inState.has(n)) {
          doc.errors.push({ line: ar.line, message: `arrow points at atom ${n}, which does not exist in this step (step ${si + 1}, line ${state.line}); give the atom a map number like [C:${n}]` });
        }
      }
    }
  });
}
