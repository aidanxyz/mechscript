// Public API. Pass in an initialised RDKit.js module.
import { parse } from './parser.js';
import { analyse } from './chem.js';
import { render } from './render.js';
import { fixCandidates, describeProblem } from './arrows.js';

export { parse } from './parser.js';
export { smilesAtoms } from './parser.js';

export function createMechScript(RDKit) {
  // Parse + chemistry checks, no drawing. Errors are sorted by line.
  function check(text) {
    const doc = parse(text);
    analyse(RDKit, doc);
    doc.errors.sort((a, b) => a.line - b.line);
    doc.warnings.sort((a, b) => a.line - b.line);
    return doc;
  }
  // Text -> { svg, width, height, doc, errors, warnings }
  // Crossings the arrow optimiser could not avoid become warnings on the arrow's line, with fixes that were
  // tried on a copy of the text and shown to remove the problem (options.suggestFixes = false to skip).
  function toSVG(text, options = {}) {
    let doc = check(text);
    let out = render(doc, options);
    if (options.autoOrient !== false && !doc.errors.length && out.report.problems.length) ({ doc, out } = autoOrient(text, doc, out, options));
    if (options.suggestFixes !== false && !doc.errors.length) {
      for (const pr of out.report.problems) {
        const verified = [];
        for (const c of fixCandidates(pr)) {
          if (verified.length >= 2) break;
          const trial = applyEdit(text, c.edit);
          const tdoc = check(trial);
          if (tdoc.errors.length) continue;
          const trep = render(tdoc, options).report;
          if (trep.problems.length < out.report.problems.length) verified.push(c.label);
        }
        const fix = verified.length ? `; try ${verified.join(' or ')}` : '; no single placement word or orientation line fixes it';
        doc.warnings.push({ line: pr.ar.line, message: describeProblem(pr) + fix });
      }
      doc.warnings.sort((a, b) => a.line - b.line);
    }
    return { ...out, doc, errors: doc.errors, warnings: doc.warnings };
  }

  // When arrows still cross something, try turning or mirroring the molecules at either end of each problem
  // arrow (the same "face"/"flip" an author could write), keep the variant with the fewest problems.
  // Molecules the author oriented themselves are left alone. Two rounds at most.
  function autoOrient(text, doc, out, options) {
    const score = (r) => r.report.problems.length * 1000 + r.report.arrowCrossings + r.report.bondCrossings + r.report.labelHits;
    const applied = []; // [{si, lay}]
    let bestScore = score(render(doc, { ...options, fastArrows: true }));
    for (let round = 0; round < 2 && out.report.problems.length; round++) {
      let best = null;
      const tried = new Set();
      for (const pr of out.report.problems) {
        const si = doc.states.indexOf(pr.state);
        const st = doc.states[si];
        const speciesOf = (n) => st.species.findIndex((sp) => sp.atoms.some((a) => a.map === n));
        const authored = new Set(st.layout.filter((l) => !l.auto).map((l) => speciesOf(l.atom)));
        for (const atom of new Set([pr.ar.from.atoms[0], ...pr.ar.to.atoms])) {
          const k = speciesOf(atom);
          if (k < 0 || authored.has(k)) continue;
          const variants = [...['right', 'left', 'up', 'down'].map((dir) => ({ type: 'face', atom, dir })), { type: 'flip', atom }];
          for (const lay of variants) {
            const key = `${si}:${k}:${lay.type}:${lay.dir || ''}`;
            if (tried.has(key)) continue;
            tried.add(key);
            // re-render the already analysed document with the extra orientation (no RDKit work)
            const extra = [...applied, { si, lay }].map((a) => ({ st: doc.states[a.si], entry: { ...a.lay, line: null, auto: true } }));
            for (const x of extra) x.st.layout.push(x.entry);
            const r2 = render(doc, { ...options, fastArrows: true });
            for (const x of extra) x.st.layout.splice(x.st.layout.indexOf(x.entry), 1);
            const sc = score(r2);
            if (sc < bestScore && (!best || sc < best.sc)) best = { sc, si, lay };
          }
        }
      }
      if (!best) break;
      applied.push({ si: best.si, lay: best.lay });
      bestScore = best.sc;
    }
    // render once more with the chosen orientations so the document's geometry matches the output
    for (const a of applied) doc.states[a.si].layout.push({ ...a.lay, line: null, auto: true });
    out = render(doc, options);
    return { doc, out };
  }
  return { check, toSVG, parse, RDKit };
}

export function formatProblems(doc) {
  return [
    ...doc.errors.map((e) => `line ${e.line}: error: ${e.message}`),
    ...doc.warnings.map((e) => `line ${e.line}: warning: ${e.message}`),
  ].join('\n');
}

function applyEdit(text, edit) {
  const lines = text.split('\n');
  if (edit.type === 'replace') lines[edit.line - 1] = edit.text;
  else if (edit.type === 'insert') lines.splice(edit.after, 0, edit.text);
  return lines.join('\n');
}
