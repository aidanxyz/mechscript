#!/usr/bin/env node
// mechscript check <file|->            print errors and warnings; exit 1 on errors
// mechscript render <file|-> [-o out.svg] [--maps] [--all-lone-pairs]
import { readFileSync, writeFileSync } from 'fs';
import initRDKit from '@rdkit/rdkit';
import { createMechScript, formatProblems } from '../src/index.js';

const [cmd, file, ...rest] = process.argv.slice(2);
if (!cmd || !['check', 'render'].includes(cmd) || !file) {
  console.error('usage: mechscript check <file|->\n       mechscript render <file|-> [-o out.svg] [--maps] [--all-lone-pairs]');
  process.exit(2);
}
const text = readFileSync(file === '-' ? 0 : file, 'utf8');
const ms = createMechScript(await initRDKit());
if (cmd === 'check') {
  // render too, so layout problems (arrows the renderer could not route cleanly) are reported with fixes
  const doc = ms.toSVG(text).doc;
  const out = formatProblems(doc);
  console.log(out || 'ok: no problems found');
  process.exit(doc.errors.length ? 1 : 0);
}
const oi = rest.indexOf('-o');
const res = ms.toSVG(text, { color: '#111', showMaps: rest.includes('--maps'), showLonePairs: rest.includes('--all-lone-pairs') ? 'all' : 'arrows' });
const problems = formatProblems(res.doc);
if (problems) console.error(problems);
if (oi >= 0 && rest[oi + 1]) writeFileSync(rest[oi + 1], res.svg);
else process.stdout.write(res.svg + '\n');
process.exit(res.errors.length ? 1 : 0);
