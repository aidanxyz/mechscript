// Reference mechanisms (the five from the spec, plus a radical example).
export const examples = [
  {
    name: 'SN2',
    code: `mechanism "SN2: hydroxide + bromomethane"

[OH-:1] + [CH3:2][Br:3]
  lp 1 -> 2
  bond 2-3 -> 3
--> "backside attack"
[OH:1][CH3:2] + [Br-:3]
`,
  },
  {
    name: 'E2',
    code: `mechanism "E2: ethoxide + 2-bromopropane"

CC[O-:1] + [H:2][CH2:3][CH:4](C)[Br:5]
  lp 1 -> 2
  bond 2-3 -> 3-4
  bond 4-5 -> 5
--> "EtOH, heat"
CC[O:1][H:2] + [CH2:3]=[CH:4]C + [Br-:5]
`,
  },
  {
    name: 'Ester hydrolysis (acid)',
    code: `mechanism "Acid-catalysed ester hydrolysis"

C[C:1](=[O:2])[O:3]C + [H:4][OH2+:5]
  lp 2 -> 4
  bond 4-5 -> 5
--> "protonation"
C[C:1](=[O+:2][H:4])[O:3]C + [OH2:6]
  lp 6 -> 1
  bond 1-2 -> 2
--> "water adds"
C[C:1]([O:2][H:4])([O+:6]([H:7])[H:8])[O:3]C + [OH2:9]
  lp 9 -> 7
  bond 7-6 -> 6
--> "proton transfer"
C[C:1]([O:2][H:4])([O:6][H:8])[O:3]C + [H:7][OH2+:9]
  lp 3 -> 7
  bond 7-9 -> 9
--> "protonation"
C[C:1]([O:2][H:4])([O:6][H:8])[O+:3]([H:7])C + [OH2:9]
  lp 2 -> 2-1
  bond 1-3 -> 3
--> "loss of methanol"
C[C:1](=[O+:2][H:4])[O:6][H:8] + C[O:3][H:7] + [OH2:9]
  lp 9 -> 4
  bond 4-2 -> 2
--> "deprotonation"
C[C:1](=[O:2])[O:6][H:8] + C[O:3][H:7] + [H:4][OH2+:9]
`,
  },
  {
    name: 'Aldol addition',
    code: `mechanism "Aldol addition (base-catalysed)"

[OH-:1] + [H:2][CH2:3][CH:4]=[O:5]
  lp 1 -> 2
  bond 2-3 -> 3-4
  bond 4-5 -> 5
--> "enolate forms"
[H:2][OH:1] + [CH2:3]=[CH:4][O-:5] + C[CH:6]=[O:7]
  lp 5 -> 5-4
  bond 4-3 -> 3-6
  bond 6-7 -> 7
--> "C–C bond forms"
[H:2][OH:1] + [O:5]=[CH:4][CH2:3][CH:6](C)[O-:7]
  lp 7 -> 2
  bond 2-1 -> 1
--> "protonation"
[OH-:1] + [O:5]=[CH:4][CH2:3][CH:6](C)[O:7][H:2]
`,
  },
  {
    name: 'Aromatic bromination',
    code: `mechanism "Electrophilic aromatic bromination"

[C:4]1([H:10])=[CH:5][CH:6]=[CH:7][CH:8]=[CH:9]1 + [Br:1][Br+:2][Fe-:3](Br)(Br)Br
  bond 4-5 -> 4-1
  bond 1-2 -> 2
--> "slow"
[C:4]1([H:10])([Br:1])[CH+:5][CH:6]=[CH:7][CH:8]=[CH:9]1 + [Br:2][Fe-:3](Br)(Br)Br
  lp 2 -> 10
  bond 10-4 -> 4-5
--> "fast"
[C:4]1([Br:1])=[CH:5][CH:6]=[CH:7][CH:8]=[CH:9]1 + [H:10][Br:2] + Br[Fe:3](Br)Br
`,
  },
  {
    name: 'Radical bromination',
    code: `mechanism "Radical bromination of ethane"

[Br:1][Br:2]
  bond 1-2 ~> 1
  bond 1-2 ~> 2
--> "light"
C[CH2:3][H:4] + [Br:1] + [Br:2]
  rad 1 ~> 1-4
  bond 3-4 ~> 1-4
  bond 3-4 ~> 3
--> "H abstraction"
C[CH2:3] + [H:4][Br:1] + [Br:2]
`,
  },
  {
    name: 'Resonance',
    code: `mechanism "Resonance in the acetate ion"

C[C:1](=[O:2])[O-:3]
  lp 3 -> 3-1
  bond 1-2 -> 2
<->
C[C:1]([O-:2])=[O:3]
`,
  },
];
