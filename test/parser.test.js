import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parse, smilesAtoms } from '../src/parser.js';

const SN2 = `mechanism "SN2"

[OH-:1] + [CH3:2][Br:3]
  lp 1 -> 2
  bond 2-3 -> 3
--> "backside attack"
[OH:1][CH3:2] + [Br-:3]`;

test('parses the SN2 example', () => {
  const d = parse(SN2);
  assert.deepEqual(d.errors, []);
  assert.equal(d.title, 'SN2');
  assert.equal(d.states.length, 2);
  assert.deepEqual(d.states[0].species.map((s) => s.smiles), ['[OH-:1]', '[CH3:2][Br:3]']);
  assert.equal(d.states[0].arrows.length, 2);
  assert.deepEqual(d.states[0].arrows[0].from, { type: 'atom', atoms: [1], source: 'lp' });
  assert.deepEqual(d.states[0].arrows[1].from, { type: 'bond', atoms: [2, 3] });
  assert.deepEqual(d.transitions, [{ line: 6, type: 'forward', label: 'backside attack' }]);
});

test('species source positions point back into the text', () => {
  const d = parse(SN2);
  for (const sp of d.states.flatMap((s) => s.species)) assert.equal(SN2.slice(sp.start, sp.end), sp.smiles);
});

test('tokenises SMILES atoms in RDKit order with map numbers', () => {
  const { atoms, error } = smilesAtoms('C[C:1](=[O:2])[O-:3]Cl');
  assert.equal(error, null);
  assert.deepEqual(atoms.map((a) => [a.symbol, a.map]), [['C', null], ['C', 1], ['O', 2], ['O', 3], ['Cl', null]]);
  assert.equal(smilesAtoms('c1ccccc1Br').atoms.length, 7);
  assert.equal(smilesAtoms('C%10CC%10').atoms.length, 3);
});

test('arrow forms: bond to bond, fishhook, placement, equilibrium and resonance', () => {
  const d = parse(`[CH2:1]=[CH:2][CH2+:3]
  bond 1-2 -> 2-3 below
<->
[CH2+:1][CH:2]=[CH2:3]
<=> "K"
[CH2+:1][CH:2]=[CH2:3]
  rad 1 ~> 1-2`);
  const ar = d.states[0].arrows[0];
  assert.deepEqual(ar.to, { type: 'bond', atoms: [2, 3] });
  assert.deepEqual(ar.place, { side: 'below' });
  assert.deepEqual(d.transitions.map((t) => t.type), ['resonance', 'equilibrium']);
  assert.equal(d.states[2].arrows[0].kind, 'fishhook');
  assert.deepEqual(d.errors, []);
});

test('reports an arrow to a missing atom with a helpful message', () => {
  const d = parse(`[OH-:1] + [CH3:2][Br:3]
  lp 1 -> 7`);
  assert.equal(d.errors.length, 1);
  assert.match(d.errors[0].message, /arrow points at atom 7, which does not exist in this step/);
  assert.equal(d.errors[0].line, 2);
});

test('reports a map number that changes element', () => {
  const d = parse(`[OH-:1] + [CH3:2][Br:3]
-->
[OH:1][CH3:2] + [Cl-:3]`);
  assert.match(d.errors[0].message, /atom 3 is Cl here but was Br on line 1/);
});

test('reports duplicate map numbers and structural mistakes', () => {
  assert.match(parse('[OH-:1] + [CH3:1]').errors[0].message, /used twice/);
  assert.match(parse('lp 1 -> 2').errors[0].message, /must come after a line of species/);
  assert.match(parse('[OH-:1]\n-->').errors[0].message, /nothing after it/);
  assert.match(parse('[OH-:1]\n  lp 1 => 2').errors[0].message, /cannot read curved arrow/);
  assert.match(parse('[OH-:1]\n[OH-:1]').errors[0].message, /two species lines in a row/);
  assert.match(parse('[OH-:1] + + C').errors[0].message, /empty species/);
});

test('comments and blank lines are ignored', () => {
  const d = parse('# a comment\n\n[OH-:1] + C   # trailing\n');
  assert.deepEqual(d.errors, []);
  assert.deepEqual(d.states[0].species.map((s) => s.smiles), ['[OH-:1]', 'C']);
});

test('a triple bond "#" in SMILES is not a comment', () => {
  const d = parse('N#[C-:1] + [CH3:2][Br:3]   # cyanide\n  lp 1 -> 2');
  assert.deepEqual(d.errors, []);
  assert.equal(d.states[0].species[0].smiles, 'N#[C-:1]');
});

test('arrow placement words', () => {
  const d = parse('[OH-:1] + [CH3:2][Br:3]\n  lp 1 -> 2 below wide\n  bond 2-3 -> 3 above');
  assert.deepEqual(d.errors, []);
  assert.deepEqual(d.states[0].arrows[0].place, { side: 'below', size: 'wide' });
  assert.deepEqual(d.states[0].arrows[1].place, { side: 'above' });
  assert.match(parse('[OH-:1] + [CH3:2]\n  lp 1 -> 2 sideways').errors[0].message, /unknown arrow option "sideways"/);
  assert.match(parse('[OH-:1] + [CH3:2]\n  lp 1 -> 2 above below').errors[0].message, /two sides/);
});

test('orientation lines and row breaks', () => {
  const d = parse('[OH-:1] + [CH3:2][Br:3]\n  face 2 left\n  flip 3\n  lp 1 -> 2\n---\n-->\n[OH:1][CH3:2] + [Br-:3]');
  assert.deepEqual(d.errors, []);
  assert.deepEqual(d.states[0].layout.map((l) => [l.type, l.atom, l.dir]), [['face', 2, 'left'], ['flip', 3, undefined]]);
  assert.equal(d.transitions[0].rowBreak, true);
  assert.match(parse('[OH-:1]\n  face 7 up').errors[0].message, /names atom 7, which does not exist/);
  assert.match(parse('[OH-:1]\n  face 1 sideways').errors[0].message, /orientation lines are/);
  assert.match(parse('face 1 up').errors[0].message, /must come after a line of species/);
});
