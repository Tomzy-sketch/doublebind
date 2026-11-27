# DoubleBind starter (Days 1 to 3 of the 7-day plan)

Research prototype that predicts how bacteria evolve resistance to Drug A and which second
drug (Drug B) catches them. This starter holds the data downloader, the mutation grammar, pocket distances,
ESM-2 scoring, ThermoMPNN stability import, baselines, the energy-to-susceptibility maths from the PRD, and
tests for all of it. The full plan is in `docs/PRD_v2.md` (Section 17 is the day-by-day workflow).

Run every command from the folder that contains `pyproject.toml`.

## Day 1 (done): setup and data

1. `python -m venv .venv`, then activate it (Windows PowerShell: `.venv\Scripts\Activate.ps1`; macOS or Linux: `source .venv/bin/activate`)
2. `pip install -e ".[dev]"`, then `python -m pytest`
3. `python -m doublebind data`, then `python -m doublebind verify` (every line should say PASS)

## Day 2 (done): pocket, ESM-2 scores and baselines

`inspect 2XCT --target gyra`, `pocket gyra`, `esm tem1`, `esm gyra`, `baselines gyra`, `baselines tem1`,
`table gyra`, `table tem1` (all with `python -m doublebind` in front). Use Kaggle for the 650M model:
upload the project folder as a dataset (without `.venv` and `data/raw/ncbi`), then
`pip install -q -e /kaggle/working/doublebind` and `python -m doublebind esm gyra --model facebook/esm2_t33_650M_UR50D`.

### Day 2 results (ESM-2 650M, pocket cutoff 8 A)

Rank of each known resistance position, 1 = best. The top 10 is what counts. Two to four known positions per
target make these anecdotes, not statistics.

| Method | GyrA S83 | GyrA D87 | TEM-1 A42 | TEM-1 E104 | TEM-1 M182 | TEM-1 G238 |
| --- | --- | --- | --- | --- | --- | --- |
| random (expected) | 437 | 437 | 143 | 143 | 143 | 143 |
| closest to the drug | 1 | 6 | 139 | 43 | 135 | 17 |
| ESM tolerated (max LLR) | 458 | 256 | 64 | 77 | 1 | 261 |
| inside the pocket, then ESM | 3 of 9 | 1 of 9 | n/a | n/a | n/a | 22 of 31 |

Variant level (ESM LLR; positive = more plausible than the wild type): GyrA S83L -6.15 and D87N -4.19;
TEM-1 A42G +0.70, E104K -1.95, M182T +5.30 (the single most plausible reachable variant of 1681) and G238S -5.94.

What this says:

- For a drug-binding target, distance to the drug finds the resistance positions, and ESM plausibility does not.
- Both known GyrA changes look implausible to ESM-2. Resistance mutations hit conserved residues and cost fitness,
  so ESM belongs in the model as a cost prior, not as the selector of resistance mutations.
- For TEM-1 (an enzyme), distance misses three of four known changes. ESM separates the stabilizing M182T
  from the destabilizing G238S, which is the trade-off Day 3 measures properly with ThermoMPNN.

### Scope, frozen at the end of Day 2

- Targets: GyrA (drug-binding resistance) and TEM-1 (enzyme, used to test the simulator on measured fitness).
- GyrA resistance candidates come from the pocket: positions within the cutoff of the drug.
- ESM-2 and ThermoMPNN are cost priors. Resistance benefit comes from the binding-change term (Day 4).
- Not this week: reinforcement learning, custom networks, metabolic models, other pathogens.

## Day 3: stability scores (ThermoMPNN)

ThermoMPNN predicts, for every possible substitution, how much it changes folding stability
(kcal/mol, positive = destabilizing, the same convention as DoubleBind). It needs a structure:

- TEM-1: the 1BTL crystal structure you already downloaded.
- GyrA: your only structure is 2XCT, a 54% fragment of an S. aureus homolog, so use the full-length
  AlphaFold DB model of E. coli GyrA (free, CC BY 4.0).

1. Update the install: `pip install -e ".[dev]"`
2. Download the AlphaFold models: `python -m doublebind data --alphafold --no-chembl`, then `python -m doublebind verify`.
   The AlphaFold lines should say PASS (the sequence matches UniProt).
3. Write the ThermoMPNN input files: `python -m doublebind stability-prep tem1` and `python -m doublebind stability-prep gyra`.
   Files land in `data/processed/stability/`. Each file has one chain, renumbered from 1, with a `.map.json`
   beside it that maps every position back to the original residue.
4. Run ThermoMPNN on each `.pdb` file (Kaggle recipe below).
5. Import the results (adjust the paths):
   `python -m doublebind stability-import tem1 --csv ThermoMPNN_inference_tem1_1BTL_A.csv`
   `python -m doublebind stability-import gyra --csv ThermoMPNN_inference_gyra_AF-P0AES4-F1_A.csv`
6. Rebuild the tables: `python -m doublebind table tem1` and `python -m doublebind table gyra`.

Day 3 is done when, for both targets, `stability-import` shows coverage of at least 95% inside the folded
region, the sign canary says PASS, and the literature-direction check for TEM-1 G238S (destabilizing) and
M182T (stabilizing) is reported. Report the direction check whether or not it agrees.

### Kaggle recipe for ThermoMPNN

1. New notebook, GPU on, Internet on. Upload the `.pdb` files from `data/processed/stability/` as a dataset.
2. First cell:
   `!git clone --depth 1 https://github.com/Kuhlman-Lab/ThermoMPNN /kaggle/working/ThermoMPNN`
   `!pip install -q "biopython==1.84" omegaconf pytorch-lightning torchmetrics wandb`
3. For each file:
   `!cd /kaggle/working/ThermoMPNN/analysis && PYTHONPATH=/kaggle/working/ThermoMPNN TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1 python custom_inference.py --pdb /kaggle/input/YOUR-DATASET/tem1_1BTL_A.pdb --chain A --out_dir /kaggle/working/out`
4. Download `/kaggle/working/out/ThermoMPNN_inference_*.csv` into your project and run step 5 above.

Known pitfalls:

- `PYTHONPATH=...` is needed because ThermoMPNN has its own `datasets.py`, and Kaggle also ships a library called `datasets`.
- `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1` avoids a `weights_only` or `Unsupported global` error when new PyTorch loads its checkpoint.
- `biopython==1.84` is pinned because ThermoMPNN imports `Bio.pairwise2`, which newer Biopython releases may remove.
- Python 3.14 is too new for most of this stack. Use Kaggle, or Python 3.10 to 3.12 locally.

## Put it on GitHub

```
git add .
git commit -m "Day 3: stability scores"
git push origin main
```

## What needs your eyes (honest list)

- The ThermoMPNN run itself was never executed from the build sandbox (no PyTorch there). The importer is tested
  against its real output format and a synthetic structure; its own commands are the likeliest place for an error.
- The AlphaFold download asks the AlphaFold DB API for the file address and falls back to common names. If it
  fails, download the PDB file from https://alphafold.ebi.ac.uk/entry/P0AES4 and save it as
  `data/raw/alphafold/AF-P0AES4-F1.pdb`.
- `data/known_mutations.csv` and `data/known_directions.csv` hold labels from memory. Add a citation to every row.
- Sign conventions in `src/doublebind/signs.py`: ThermoMPNN is confirmed (positive = destabilizing); the raw
  Megascale file is not yet.
- CARD data is free for academic use only (`--card` downloads it); read its license before any commercial use.

## Folder map

```
src/doublebind/
  constants.py            physical constants
  signs.py                energy sign conventions and direction checks
  csi/csi.py              binding energy -> IC50 fold change -> CSI, and the confidence rule
  csi/categories.py       collateral sensitivity / resistance-sparing / cross-resistance / inconclusive
  data/targets.py         GyrA, GyrB, TEM-1 definitions (structures, stability source)
  data/fetch.py           downloads (UniProt, NCBI, RCSB, ChEMBL, AlphaFold) and verification
  data/codons.py          genetic code, single-nucleotide variants, numbering offsets
  data/manifest.py        checksums for every raw file
  structure/proximity.py  chains, ligands, alignment to the target, residue-to-drug distances
  structure/prep.py       clean ThermoMPNN input file plus a map back to the original residues
  scoring/esm_llr.py      ESM-2 log-likelihood ratios
  stability.py            ThermoMPNN import, coverage, sign canary, literature directions
  baselines.py            Day 2 baselines
  tables.py               the merged variant table
  cli.py                  python -m doublebind <data|verify|variants|numbers|inspect|pocket|esm|baselines|table|stability-prep|stability-import>
tests/                    offline tests (synthetic genes, structures and logits)
data/raw/                 downloads land here (not committed)
data/processed/           generated tables (not committed)
data/known_mutations.csv
data/known_directions.csv
docs/PRD_v2.md           the full PRD
```

## Next

Day 4 builds the binding proxy: dock the approved drugs from ChEMBL into the pocket, re-score each pocket
variant, and compare with a protein-ligand mutation dataset. Send the output of `stability-import` for both
targets and the next piece gets built on top of it.
