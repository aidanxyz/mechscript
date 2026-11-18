# MechScript

A small text language for reaction mechanisms with curved arrows, plus a renderer that draws it as SVG. AI models write it well (it is SMILES plus one line per arrow), and code can be drawn exactly.

```mechanism
mechanism "SN2: hydroxide + bromomethane"

[OH-:1] + [CH3:2][Br:3]
  lp 1 -> 2
  bond 2-3 -> 3
--> "backside attack"
[OH:1][CH3:2] + [Br-:3]
```

- **Live page:** https://aidanxyz.github.io/mechscript/ (send the link to an AI assistant and ask it to draw a mechanism there; it replies with a link to the drawing)
- **Language spec (one page, also the prompt for Claude):** [spec/MECHSCRIPT.md](spec/MECHSCRIPT.md)
- **Playground:** `playground/index.html` (code on the left, live drawing on the right)
- **Reference mechanisms:** `src/examples.js` (SN2, E2, acid ester hydrolysis, aldol addition, aromatic bromination, radical bromination, resonance)

## Run it

```bash
npm install
npm run build      # bundles dist/mechscript.js and copies RDKit next to it
npm run serve      # http://localhost:8377/playground/  (also /test/gallery.html, /test/autorun.html)
npm test
```

## The public page

`tools/build-site.js` builds `site/`: one page that reads the mechanism from the link (everything after `#`, percent-encoded), draws it, and lists problems. The instructions for AI assistants and the full spec are baked into the page as static text (also at `/llms.txt`), so an assistant that fetches the link can read how to use it. Every push to `main` tests, builds and deploys it to GitHub Pages (`.github/workflows/pages.yml`).

## Use it on a page (like Mermaid)

```html
<script src="dist/mechscript.js"></script>
```

Every `<pre><code class="language-mechanism">` (what Markdown renders from a ```` ```mechanism ```` block), `<pre class="mechanism">` or `<div class="mechanism">` is replaced by an SVG drawing. RDKit.js loads from unpkg by default. To use a local copy, add `data-rdkit="dist/RDKit_minimal.js"` to the script tag. The API is `MechScript.render(text, options)`, which returns `{svg, errors, warnings}`, plus `MechScript.run()` and `MechScript.toPNG(svg)`.

## Command line (for an AI fixing its own output)

```bash
node bin/mechscript.js check mechanism.mech       # line-numbered errors and drawing warnings with tested fixes; exit 1 on errors
node bin/mechscript.js render mechanism.mech -o out.svg
```

## Reviewing renders without a server

```bash
node tools/sheet.js            # renders every mechanism into review/sheet.html (optimiser off vs on, crossing counts)
node tools/sheet.js --after --only cases
```

Open `review/sheet.html` straight from disk. It contains only SVG, with no scripts and no RDKit in the browser.

## Measuring drawing quality

```bash
npm run quality          # scores every mechanism, writes review/quality.html (worst first, issues circled)
npm run quality:check    # exits 1 if any mechanism scores worse than test/quality/baseline.json
node tools/quality.js --save-baseline   # accept the current scores as the new baseline
```

It covers the reference examples, `test/ai`, `test/cases` and the extra fixtures in `test/quality` (rings, fused and bridged systems, ring closures, pericyclic steps). Each drawing is measured for ring shape, bond lengths, bond angles, atoms or molecules too close together, molecules changing shape or turning between steps, and arrow crossings, detours, missed targets, arrows running too close to each other or skimming labels, bonds and "+" signs, arrows that bend back to touch themselves, and S-shaped arrows. Thresholds are at the top of `tools/quality.js`.

## How it works

| File | Role |
| --- | --- |
| `src/parser.js` | Text → steps, species (SMILES), arrows, reaction arrows. No dependencies. Checks map numbers. |
| `src/chem.js` | RDKit.js builds each molecule (Kekulé, explicit H kept) with 2D coordinates. Checks that bonds exist and that atoms have lone pairs. Lays molecules out as part of the ring a pericyclic or intramolecular step closes. |
| `src/arrows.js` | Curved-arrow placement: candidate curves per arrow, scored for crossings (arrows, bonds, labels, ring faces), cheapest combination chosen. If crossings remain, `src/index.js` first tries turning and mirroring the molecules involved and keeps the best drawing; anything still left becomes a warning with fixes that were tested on a copy. |
| `src/render.js` | Layout and SVG: species left to right with wrapping, molecules aligned to the previous step, lone pairs on the side facing their arrow, Bézier arrows bowing away from the molecules. |
| `playground/` | Code on the left, live drawing on the right; examples, error list, SVG/PNG export, share link. |

## Decisions on the spec's open questions

- **Partial charges and stereo wedges:** not in v1. Bonds are drawn flat.
- **Visual editor:** dropped for now; MechScript is a text → drawing tool.
- **Publishing:** MIT-licensed and laid out as an npm package (`bin`, `main`). Nothing has been published.
