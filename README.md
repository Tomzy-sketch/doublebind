# DoubleBind starter (Days 1 and 2 of the 7-day plan)

Research prototype that predicts how bacteria evolve resistance to Drug A and which second
drug (Drug B) catches them. This starter holds the Day 1 data downloader, the Day 2 mutation
grammar, pocket distances, ESM-2 scoring and baselines, the energy-to-susceptibility maths from
the PRD, and tests for all of it. The full plan is in `docs/PRD_v2.md` (Section 17 is the day-by-day workflow).

## Day 1 (done): setup and data

1. **Python 3.10 or newer.** Check with `python --version`.
2. **Environment:** Windows PowerShell `python -m venv .venv` then `.venv\Scripts\Activate.ps1`
   (if blocked: `Set-ExecutionPolicy -Scope Process RemoteSigned`); macOS or Linux `python3 -m venv .venv && source .venv/bin/activate`
3. `pip install -e ".[dev]"`, then `python -m pytest`
4. `python -m doublebind data`, then `python -m doublebind verify` (every line should say PASS)

Day 1 is done when `verify` shows PASS for the manifest, the accessions and every coding sequence.

## Day 2: pocket, ESM-2 scores and baselines

First, after copying the update files over your project, run `pip install -e ".[dev]"` again (it adds numpy and Biopython).

1. **Look inside the structure** and see how it lines up with E. coli GyrA:
   `python -m doublebind inspect 2XCT --target gyra`
   You should see the protein chains with an identity and a coverage, and one candidate ligand (the quinolone).
2. **Measure how far each residue is from the drug:**
   `python -m doublebind pocket gyra`
   The last lines show `known S83L` and `known D87N` with their distances. If the wrong molecule was chosen as
   the drug, rerun with `--ligand NAME` (the names are printed by step 1). `python -m doublebind pocket tem1`
   only works if the TEM-1 structure holds a bound ligand; skip it if it says "no ligand found".
3. **Score every substitution with ESM-2.** Install the model libraries once (large download, skip if you
   will use Kaggle): `pip install -e ".[esm]"`. Smoke test first (small model, about a minute):
   `python -m doublebind esm tem1 --model facebook/esm2_t12_35M_UR50D`
   Then the real runs. TEM-1 takes a few minutes on a laptop CPU, GyrA 20 to 60 minutes:
   `python -m doublebind esm tem1` and `python -m doublebind esm gyra`
   Quick GyrA check while you wait: add `--mode wt` (one pass, seconds, slightly less accurate).
   For the 650M model use Kaggle (recipe below).
4. **Compare the baselines:** `python -m doublebind baselines gyra` and `python -m doublebind baselines tem1`
5. **Build the Day 2 table:** `python -m doublebind table gyra` and `python -m doublebind table tem1`

Day 2 is done when `data/processed` holds `*_variant_table.csv` and `*_baselines.json` for GyrA and TEM-1.

### How to read the baselines

- Positions are ranked 1 (best) to N. What counts is whether S83 and D87 land in the top 10.
- **random** is the expected rank (N+1)/2. **proximity only** ranks residues by distance to the drug.
  **ESM tolerated** ranks positions by their most plausible reachable substitution. **ESM constrained**
  is the opposite orientation, kept as a diagnostic. **proximity then ESM** is the first real candidate:
  inside the pocket radius, rank by ESM plausibility.
- With two known positions these numbers are anecdotes, not statistics. They tell you whether the pipeline
  behaves sensibly. Do not tune the cutoff or the scorer to make them look better.
- Distances are measured on a homolog structure (2XCT is believed to be S. aureus), mapped onto E. coli GyrA by
  alignment. Check the identity and coverage that `inspect` and `pocket` print.

### Kaggle recipe for the 650M model

1. kaggle.com > Create > New Notebook. Settings: Accelerator GPU, Internet On (needs phone verification).
2. Upload your project folder as a dataset, but zip it **without** `.venv`, `data/raw/ncbi`, `data/raw/chembl`
   and `data/raw/pdb` (they are large and not needed for scoring). `data/raw/cds` must be included.
3. First cell: `!cp -r /kaggle/input/YOUR-DATASET/doublebind /kaggle/working/ && pip install -q -e /kaggle/working/doublebind`
4. Next cells: `!cd /kaggle/working/doublebind && python -m doublebind esm gyra --model facebook/esm2_t33_650M_UR50D`
   and the same for `tem1`.
5. Download the new `data/processed/*_esm_*.csv` files from the Output panel into your local `data/processed` folder,
   then rerun `baselines` and `table` locally. The baselines pick the biggest model they find.

## Put it on GitHub

```
git init
git add .
git commit -m "Day 1 and 2"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/doublebind.git
git push -u origin main
```

## What needs your eyes (honest list)

- The ESM scoring wrapper and the real-structure code have never run on a real model or structure from this
  build sandbox (no access to Hugging Face or RCSB). They are tested with fake logits and synthetic structures.
  If a step fails, copy the error back to Claude.
- Accessions, protein lengths and PDB IDs in `src/doublebind/data/targets.py` are confirmed for GyrA, GyrB and TEM-1 by your Day 1 run.
- `data/known_mutations.csv` has labels from memory. Add a citation to every row before you rely on them.
- Sign conventions in `src/doublebind/signs.py` are marked `confirmed=False` until you read each tool's documentation.
- CARD data is free for academic use only (`--card` downloads it); read its license before any commercial use.

## Folder map

```
src/doublebind/
  constants.py          physical constants
  signs.py              energy sign conventions and direction checks
  csi/csi.py            binding energy -> IC50 fold change -> CSI, and the confidence rule
  csi/categories.py     collateral sensitivity / resistance-sparing / cross-resistance / inconclusive
  data/targets.py       GyrA, GyrB, TEM-1 definitions
  data/fetch.py         downloads and verification
  data/codons.py        genetic code, single-nucleotide variants, numbering offsets
  data/manifest.py      checksums for every raw file
  structure/proximity.py  chains, ligands, alignment to the target, residue-to-drug distances
  scoring/esm_llr.py    ESM-2 log-likelihood ratios (numpy core + Hugging Face wrapper)
  baselines.py          Day 2 baselines
  tables.py             the merged variant table
  cli.py                python -m doublebind <data|verify|variants|numbers|inspect|pocket|esm|baselines|table>
tests/                  offline tests (synthetic genes, structures and logits)
data/raw/               downloads land here (not committed)
data/processed/         generated tables (not committed)
data/known_mutations.csv
docs/PRD_v2.md          the full PRD
```

## Next

Day 2 adds ESM-2 plausibility scores for every reachable variant (run on a Kaggle GPU) and the three baselines.
Send the output of `verify` and `variants` and the next piece gets built on top of it.
