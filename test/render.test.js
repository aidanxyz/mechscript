import { test, before } from 'node:test';
import assert from 'node:assert/strict';
import initRDKit from '@rdkit/rdkit';
import { createMechScript } from '../src/index.js';
import { examples } from '../src/examples.js';
import { readdirSync, readFileSync } from 'fs';

let ms;
before(async () => { ms = createMechScript(await initRDKit()); });

for (const ex of examples) {
  test(`reference mechanism renders without problems: ${ex.name}`, () => {
    const r = ms.toSVG(ex.code);
    assert.deepEqual(r.errors, [], JSON.stringify(r.errors));
    assert.deepEqual(r.warnings, [], JSON.stringify(r.warnings));
    assert.match(r.svg, /^<svg /);
    const arrows = ex.code.split('\n').filter((l) => /^\s*(lp|bond|rad)\b/.test(l)).length;
    assert.equal((r.svg.match(/class="ms-arrow"/g) || []).length, arrows, 'one curved arrow per arrow line');
    assert.ok(!/NaN|Infinity/.test(r.svg), 'no NaN coordinates');
  });
}

test('chemistry checks: bond that does not exist', () => {
  const d = ms.check('[OH-:1] + [CH3:2][Br:3]\n  bond 1-2 -> 2');
  assert.match(d.errors[0].message, /no bond between atoms 1 and 2/);
});

test('chemistry checks: lone pair arrow from an atom without lone pairs', () => {
  const d = ms.check('[CH4:1] + [Br:2][Br:3]\n  lp 1 -> 2');
  assert.match(d.warnings[0].message, /atom 1 \(C\) has no lone pair/);
});

test('chemistry checks: invalid SMILES', () => {
  const d = ms.check('C1CC + [OH-:1]');
  assert.match(d.errors[0].message, /RDKit cannot read SMILES "C1CC"/);
});

test('rendering is deterministic', () => {
  assert.equal(ms.toSVG(examples[2].code).svg, ms.toSVG(examples[2].code).svg);
});

test('automatic placement: no curved arrows cross, none run through labels (all test mechanisms)', () => {
  const files = ['test/ai', 'test/cases'].flatMap((d) => readdirSync(d).filter((f) => f.endsWith('.mech')).map((f) => readFileSync(`${d}/${f}`, 'utf8')));
  let bondCrossings = 0;
  let labelHits = 0;
  for (const code of [...examples.map((e) => e.code), ...files]) {
    const r = ms.toSVG(code);
    assert.deepEqual(r.errors, []);
    assert.equal(r.report.arrowCrossings, 0, code.split('\n')[0]);
    bondCrossings += r.report.bondCrossings;
    labelHits += r.report.labelHits;
  }
  assert.equal(bondCrossings, 0, `bond crossings: ${bondCrossings}`);
  assert.equal(labelHits, 0, `label hits: ${labelHits}`);
});

test('placement words choose the side of the curve', () => {
  const base = '[OH-:1] + [CH3:2][Br:3]\n  lp 1 -> 2 PLACE\n  bond 2-3 -> 3';
  const below = ms.toSVG(base.replace('PLACE', 'below')).doc.states[0].arrows[0].geo;
  const above = ms.toSVG(base.replace('PLACE', 'above')).doc.states[0].arrows[0].geo;
  const midY = (g) => (g.S.y + 2 * g.C.y + g.E.y) / 4;
  assert.ok(midY(below) > Math.max(below.S.y, below.E.y), 'below bows down');
  assert.ok(midY(above) < Math.min(above.S.y, above.E.y), 'above bows up');
});

test('a Diels-Alder diene is drawn curled towards the dienophile (ring-aware layout)', () => {
  const code = readFileSync('test/cases/diels_alder.mech', 'utf8');
  const doc = ms.toSVG(code).doc;
  const diene = doc.states[0].species[0];
  const byMap = (m) => diene.geo.P[diene.mol.atoms.findIndex((a) => a.map === m)];
  const d14 = Math.hypot(byMap(1).x - byMap(4).x, byMap(1).y - byMap(4).y);
  assert.ok(d14 < 34 * 2.2, `C1 and C4 close (s-cis): ${d14.toFixed(1)}`);
});

test('face turns a molecule so the atom points that way; flip mirrors it; both carry to later steps', () => {
  const code = (extra) => `C[C:1](=[O:2])[O:3]C + [OH-:4]\n${extra}  lp 4 -> 1\n  bond 1-2 -> 2\n-->\nC[C:1]([O-:2])([OH:4])[O:3]C`;
  const pos = (doc, si, m) => { const sp = doc.states[si].species.find((x) => x.mol.atoms.some((a) => a.map === m)); const i = sp.mol.atoms.findIndex((a) => a.map === m); const P = sp.geo.P; const c = P.reduce((t, p) => ({ x: t.x + p.x / P.length, y: t.y + p.y / P.length }), { x: 0, y: 0 }); return { x: P[i].x - c.x, y: P[i].y - c.y }; };
  const up = ms.toSVG(code('  face 2 up\n')).doc;
  const v = pos(up, 0, 2);
  assert.ok(v.y < 0 && Math.abs(v.x) < 1e-6 * 1000, `O2 points up: ${JSON.stringify(v)}`);
  const w = pos(up, 1, 2);
  assert.ok(w.y < 0 && Math.abs(w.y) > Math.abs(w.x), 'orientation carries to the next step');
  const plain = pos(ms.toSVG(code('')).doc, 0, 3);
  const flipped = pos(ms.toSVG(code('  flip 3\n')).doc, 0, 3);
  assert.ok(Math.abs(plain.y + flipped.y) < 1e-6 || Math.sign(plain.y) !== Math.sign(flipped.y), 'flip mirrors');
});

test('"---" puts the next step on a new row', () => {
  const base = '[OH-:1] + [CH3:2][Br:3]\n  lp 1 -> 2\n  bond 2-3 -> 3\nBREAK-->\n[OH:1][CH3:2] + [Br-:3]';
  const one = ms.toSVG(base.replace('BREAK', ''));
  const two = ms.toSVG(base.replace('BREAK', '---\n'));
  assert.ok(two.height > one.height * 1.5, `${two.height} vs ${one.height}`);
});

test('the renderer turns molecules itself when that clears a crossing (bromonium)', () => {
  const r = ms.toSVG(readFileSync('test/cases/bromonium.mech', 'utf8'));
  assert.equal(r.report.problems.length, 0);
  assert.deepEqual(r.warnings, []);
});

test('crossings that remain become warnings with fixes that were verified', () => {
  const code = readFileSync('test/cases/bromonium.mech', 'utf8');
  const r = ms.toSVG(code, { autoOrient: false });
  const w = r.warnings.find((x) => /runs through/.test(x.message));
  assert.ok(w, JSON.stringify(r.warnings));
  const m = w.message.match(/try "face (\d+) (\w+)"/);
  assert.ok(m, w.message);
  const fixed = code.replace('[CH2:1]=[CH2:2] + [Br:3][Br:4]', `[CH2:1]=[CH2:2] + [Br:3][Br:4]\n  face ${m[1]} ${m[2]}`);
  assert.deepEqual(ms.toSVG(fixed, { autoOrient: false }).warnings, []);
});
