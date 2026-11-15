// Bundles the browser script and copies RDKit next to it (for offline use / the playground).
import { build } from 'esbuild';
import { cpSync, mkdirSync } from 'fs';

mkdirSync('dist', { recursive: true });
await build({ entryPoints: ['src/browser.js'], bundle: true, format: 'iife', outfile: 'dist/mechscript.js', target: 'es2020', minify: false });
await build({ entryPoints: ['src/browser.js'], bundle: true, format: 'iife', outfile: 'dist/mechscript.min.js', target: 'es2020', minify: true });
for (const f of ['RDKit_minimal.js', 'RDKit_minimal.wasm']) cpSync(`node_modules/@rdkit/rdkit/dist/${f}`, `dist/${f}`);
console.log('built dist/');
