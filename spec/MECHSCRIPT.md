# MechScript

MechScript is a plain-text language for reaction mechanisms with curved (electron-pushing) arrows. Write *what* happens; the renderer decides where everything goes. Put it in a fenced code block tagged `mechanism`.

```mechanism
mechanism "SN2: hydroxide + bromomethane"

[OH-:1] + [CH3:2][Br:3]
  lp 1 -> 2
  bond 2-3 -> 3
--> "backside attack"
[OH:1][CH3:2] + [Br-:3]
```

## Structure

A mechanism is a sequence of **steps**. Each step is one line of species, then (indented, optional) the curved arrows acting on those species. Steps are joined by reaction-arrow lines. The products of one step are the starting species of the next.

| Line | Meaning |
| --- | --- |
| `mechanism "Title"` | Optional, first line only. |
| `A + B + C` | Species of one step, as SMILES, separated by ` + ` (spaces around the `+`). |
| `lp 1 -> 2` | A lone pair on atom 1 forms a bond to atom 2. |
| `lp 1 -> 1-2` | A lone pair on atom 1 becomes a (new π) bond between 1 and 2. |
| `bond 2-3 -> 3` | The electrons of bond 2–3 move onto atom 3 (the bond breaks). |
| `bond 1-2 -> 2-3` | The electrons of bond 1–2 form a bond 2–3 (π bond shift, or attack on atom 3). |
| `rad 1 ~> 1-2` | One unpaired electron on atom 1 moves (fishhook). |
| `~>` instead of `->` | Fishhook arrow: one electron (radical chemistry). |
| `--> "label"` | Reaction arrow; the optional label (reagents, conditions) sits above it. |
| `<=> "label"` | Equilibrium arrow. |
| `<->` | Resonance arrow (the next line is a resonance form, not a new step). |
| `# text` | Comment. |

Arrow sources: `lp N` (lone pair), `rad N` (single electron), `bond N-M`. Arrow targets: an atom `N` or a bond/position between two atoms `N-M` (the bond may be new, even between two different species).

## Rules

1. **SMILES, Kekulé form.** Write rings with explicit double bonds (`C1=CC=CC=C1`, not `c1ccccc1`) so arrows can name a specific π bond.
2. **Number every atom an arrow touches** with an atom-map number: `[CH3:2]`, `[O-:1]`, `[H:7]`. Unnumbered atoms are fine elsewhere (`CC[O-:1]`).
3. **Numbers are permanent.** An atom keeps its number in every step of the mechanism, and two atoms in one step never share a number.
4. **Hydrogens that move must be explicit atoms** with a number, e.g. `[H:2][CH2:3]C`, and appear as `[H:2]` in later steps too.
5. **Charges and radicals go in brackets**: `[OH-:1]`, `[OH2+:5]`, `[C+:3]`, `[CH3:1]` (a methyl radical: bracket atom with too few H). Write H counts inside brackets explicitly.
6. **Every step must balance**: the next step's species line must show exactly the result of the arrows (new bonds, broken bonds, changed charges).
7. Arrows belong to the species line directly above them. Never put two species lines in a row; use `-->` between steps.
8. **Species may come and go between steps.** A new reagent (a second water, a base) can first appear in the step where it reacts, and spectators can be dropped. Only the atoms that keep a number must stay consistent.
9. **An implicit H may become explicit later**: `[OH2:6]` in one step and `[H:12][OH:6]` in the next is fine, as long as the new number is unused.
10. **Bracket H counts exclude explicit H atoms**: `[H:2][CH2:3]C` is a CH3 group with one numbered H.

## Arrow details

- `bond N-M` on a double or triple bond means its π electrons; the bond order drops by one.
- Atom order inside `N-M` does not matter: `bond 1-2 -> 1-8` and `bond 2-1 -> 8-1` are the same arrow.
- Put an arrow under the form or step where the electrons start. For resonance, the arrows converting one form into the next go under the earlier form; the next line after `<->` repeats every species (spectators too). The last resonance form may carry the next step's arrows and continue with `-->`.
- `-->` is one elementary step; use `<=>` for a reversible step (proton transfers, equilibria). Labels are optional on both.
- Comments (`#` at the start of a line or after a space) may go anywhere. A `#` inside SMILES is a triple bond.
- The checker checks map numbers, bonds and lone pairs, not formal-charge bookkeeping or metal valences. Draw the conventional mechanism (e.g. a lone pair on Br of FeBr4- attacking H).

## Steering the drawing (optional)

Leave these out. The renderer turns molecules to face their partners and routes arrows so they do not cross each other, bonds or labels. Add them only when the checker warns about the drawing, using the fix it suggests:

| Line or words | Effect |
| --- | --- |
| `above`, `below`, `left`, `right` after an arrow | The side the curve bows towards. |
| `tight`, `wide` after an arrow | A flatter or a rounder curve. |
| `face 3 right` under a species line | Turn the molecule containing atom 3 so atom 3 points right (also `left`, `up`, `down`). Later steps keep the orientation. |
| `flip 3` under a species line | Mirror the molecule containing atom 3. |
| `---` between steps | Start a new row at the next step. |

Example: `lp 1 -> 2 below wide`.

## A multi-step example

```mechanism
mechanism "E2"

CC[O-:1] + [H:2][CH2:3][CH:4](C)[Br:5]
  lp 1 -> 2
  bond 2-3 -> 3-4
  bond 4-5 -> 5
--> "EtOH, heat"
CC[O:1][H:2] + [CH2:3]=[CH:4]C + [Br-:5]
```

## Errors

The checker reports problems by line, for example `line 4: error: arrow points at atom 7, which does not exist in this step`. Fix the line it names and check again. Drawing problems are warnings that come with fixes the checker has already tried, for example `line 4: warning: arrow runs through the C (atom 1) label; try "face 1 right" under the species line`. Typical causes: a missing map number, a map number that changed between steps, an arrow out of a bond that does not exist, or invalid SMILES (unbalanced brackets, wrong valence).
