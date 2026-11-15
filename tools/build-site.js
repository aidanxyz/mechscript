// Builds the public web page into site/: index.html (instructions + language reference baked in as static
// text, so an AI assistant that fetches the page can read them), app.js, RDKit, llms.txt.
// Usage: SITE_URL=https://example.github.io/mechscript/ node tools/build-site.js
import { build } from 'esbuild';
import { readFileSync, writeFileSync, mkdirSync, cpSync, rmSync } from 'fs';
import { examples } from '../src/examples.js';

const SITE = (process.env.SITE_URL || 'http://localhost:8377/site/').replace(/\/?$/, '/');
const out = 'site';
rmSync(out, { recursive: true, force: true });
mkdirSync(out, { recursive: true });

const esc = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const encodeText = (t) => encodeURIComponent(t).replace(/[()!'*~]/g, (c) => '%' + c.charCodeAt(0).toString(16).toUpperCase());

const spec = readFileSync('spec/MECHSCRIPT.md', 'utf8');
const example = examples.find((e) => e.name === 'E2').code.trim();
const exampleLink = SITE + '#' + encodeText(example + '\n');

const html = readFileSync('tools/site/index.html', 'utf8')
  .replaceAll('{{SITE}}', SITE)
  .replaceAll('{{EXAMPLE_LINK}}', esc(exampleLink))
  .replaceAll('{{EXAMPLE}}', esc(example))
  .replaceAll('{{SPEC}}', esc(spec));
writeFileSync(`${out}/index.html`, html);

const llms = `# MechScript

> Draws organic reaction mechanisms (curved electron-pushing arrows) from a small text language.
> Page: ${SITE}

## How to draw a mechanism for a user

1. Write the mechanism in MechScript (reference below).
2. Link: ${SITE}# followed by the MechScript text percent-encoded like encodeURIComponent (newline %0A, space %20, # %23, " %22), plus ( %28 and ) %29.
3. Reply with the link as a clickable Markdown link and the code in a \`\`\`mechanism block.
4. The page lists problems by line, with suggested fixes; if the user pastes them back, fix the code and send a new link.

Example code:

\`\`\`mechanism
${example}
\`\`\`

Example link: ${exampleLink}

---

${spec}`;
writeFileSync(`${out}/llms.txt`, llms);
writeFileSync(`${out}/spec.md`, spec);

await build({ entryPoints: ['tools/site/app.js'], bundle: true, format: 'iife', minify: true, outfile: `${out}/app.js`, target: 'es2020' });
for (const f of ['RDKit_minimal.js', 'RDKit_minimal.wasm']) cpSync(`node_modules/@rdkit/rdkit/dist/${f}`, `${out}/${f}`);
writeFileSync(`${out}/.nojekyll`, '');
console.log(`built ${out}/ for ${SITE}`);
