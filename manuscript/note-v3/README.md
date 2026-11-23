# OEIS A240513 ladder formulas and related array recurrences

This package proves OEIS A240513, revision 11, FORMULA entries 1–3, and
extends Bala’s (2014) treatment of the related norm sequence to binary ladders.
It also proves the recurrences for widths three and four in A240514–A240515
and A240519, with exact minimal orders 20 and 48. The note gives corrected
Chebyshev identities, their transform, and further norm identities.

## Contents

- `note.tex`: the main mathematical source, with its bibliography.
- `note.pdf`: the typeset note, supplied alongside the source archive.
- `metadata.json`: stable identity, creators, license and mathematical coverage.
- `coverage.json`: exact posing locator, proof scopes and sibling coverage.
- `verify.py`: standard-library Python checks in exact arithmetic.
- `width-certificates.json`: the wider recurrence coefficients, initial values and Hankel witnesses.
- `checks.json`: deterministic output of the default checker command.
- `MANIFEST.json`: hashes and sizes of the scientific source files in the archive.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Or use a standard LaTeX installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is inside `note.tex`; no external bibliography file is required.

## Check

Python 3.10 or newer is sufficient; no third-party modules are needed:

```sh
python3 verify.py --output checks.json
```

The expected status is `pass`. The checker covers all 16 ladder states,
literal arrays through height 8, norm identities through index 210, all
rectangles with height and width from 1 to 4, and the exhaustive annihilator
identities on 64 and 256 wider states. The modular Hankel determinants are
306890401 and 347540197 modulo 1000000007. The written transfer and Hankel
arguments explain why these finite identities imply recurrences at every
height and their exact minimal eventual orders.

Optional comparisons with the published OEIS tables are:

```sh
python3 verify.py --bfile b240513.txt --width-three-bfile b240514.txt \
  --width-four-bfile b240515.txt --output external-checks.json
```

Obtain those tables from `https://oeis.org/A240513/b240513.txt`,
`https://oeis.org/A240514/b240514.txt`, and
`https://oeis.org/A240515/b240515.txt`. They are external comparison data,
not dependencies of the default checker or the LaTeX build.

## Mathematical coverage

The pinned declaration `HardinA240513Ladder.result` proves the literal
two-column count's recurrence for every height at least 6 and the formal
generating-function identity. Its declaration `finite_certificate` states
the sixteen-state vector identity:

https://github.com/the-omega-institute/trureturing/blob/70badecc7fca07aa0b70e1202720702f7d1f9670/D5/S3/Combinatorics/TraceAvoidingArrays/HardinA240513Ladder.lean

The quaternion relation, Chebyshev and auxiliary identities, transposition,
Fibonacci formula and fixed-width rationality have written proofs. The
width-three/four recurrences and minimality combine written proofs with
exhaustive finite integer/modular identities. Only the two-column recurrence
and generating function are covered by the linked Lean declaration.

The coverage table distinguishes asserted empirical formulas, corrected
source statements, previously proved norm-sequence results, and additional
consequences. The title credits Bala (2014); the text also credits
Chen–Shapiro (2007), Somos (2016) and Mathar (2024) at their exact scopes.

## Creators and License

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

CC BY 4.0 (`cc-by-4.0`).
