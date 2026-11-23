# DoubleBind starter (Days 1 and 2 of the 7-day plan)

Research prototype that predicts how bacteria evolve resistance to Drug A and which second
drug (Drug B) catches them. This starter holds the Day 1 data downloader, the Day 2 mutation
grammar, the energy-to-susceptibility maths from the PRD, and tests for all of it.
The full plan is in `docs/PRD_v2.md` (Section 17 is the day-by-day workflow).

## Do this first (about 45 minutes)

1. **Python 3.10 or newer.** Check with `python --version`. On Windows, tick "Add Python to PATH" when installing.
2. **Open this folder in a terminal** (VS Code: Terminal > New Terminal) and create an environment:
   - Windows PowerShell: `python -m venv .venv` then `.venv\Scripts\Activate.ps1`
     (if PowerShell blocks it, run `Set-ExecutionPolicy -Scope Process RemoteSigned` first)
   - macOS or Linux: `python3 -m venv .venv && source .venv/bin/activate`
3. **Install:** `pip install -e ".[dev]"`
4. **Run the tests:** `python -m pytest` (expect every test to pass, none touch the internet)
5. **Check the PRD numbers:** `python -m doublebind numbers`
6. **Download the data:** `python -m doublebind data` (needs internet, a few MB, a few minutes)
7. **Check what arrived:** `python -m doublebind verify` and read every WARN or FAIL line
8. **Day 2 preview:** `python -m doublebind variants tem1` then `python -m doublebind variants gyra`

If step 6 is blocked on your network, run steps 3 to 7 in a Kaggle notebook with Internet turned on
(upload this folder as a dataset, or `git clone` your GitHub copy), then download `data/raw`.
To go faster you can add `--no-chembl` to step 6 (ChEMBL is only needed on Day 4).
NCBI is polite to scripts that identify themselves: set the `NCBI_EMAIL` environment variable if you can.

## Day 1 is done when

- `python -m doublebind verify` shows the manifest check as PASS (every raw file matches its checksum).
- The three UniProt accessions (P0AES4 GyrA, P0AES6 GyrB, P62593 TEM-1) print the protein you expect.
- Each target's coding sequence says PASS (translates exactly to the protein) or a WARN whose differences you understand.
- The tests pass and the repo is on GitHub.

## Put it on GitHub

Create an empty public repository called `doublebind` on github.com (public repos get free CI minutes), then:

```
git init
git add .
git commit -m "Day 1 starter"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/doublebind.git
git push -u origin main
```

## What needs your eyes (honest list)

- The downloader was written without access to UniProt, NCBI, RCSB or ChEMBL, so it has never run
  against the live servers. If a step fails, copy the printed URL and error back to Claude.
- Accessions, protein lengths and PDB IDs in `src/doublebind/data/targets.py` came from memory.
  `verify` prints what the servers actually say; fix the file if they disagree.
- `2XCT` is believed to be an S. aureus gyrase structure, a homolog of E. coli GyrA with different residue numbers.
  `verify` prints each structure's title so you can check.
- `data/known_mutations.csv` has labels from memory (GyrA S83L, D87N; TEM-1 A42G, E104K, M182T, G238S).
  Add a citation to every row before you rely on them. TEM-1 literature uses Ambler numbering, and the
  `variants` command finds the offset to UniProt numbering for you.
- Sign conventions in `src/doublebind/signs.py` are marked `confirmed=False` until you read each tool's documentation.
- CARD data is free for academic use only (`--card` downloads it); read its license before any commercial use.

## Folder map

```
src/doublebind/
  constants.py        physical constants
  signs.py            energy sign conventions and direction checks
  csi/csi.py          binding energy -> IC50 fold change -> CSI, and the confidence rule
  csi/categories.py   collateral sensitivity / resistance-sparing / cross-resistance / inconclusive
  data/targets.py     GyrA, GyrB, TEM-1 definitions
  data/fetch.py       downloads and verification
  data/codons.py      genetic code, single-nucleotide variants, numbering offsets
  data/manifest.py    checksums for every raw file
  cli.py              python -m doublebind <data|verify|variants|numbers>
tests/                offline tests
data/raw/             downloads land here (not committed)
data/processed/       generated tables (not committed)
data/known_mutations.csv
docs/PRD_v2.md        the full PRD
```

## Next

Day 2 adds ESM-2 plausibility scores for every reachable variant (run on a Kaggle GPU) and the three baselines.

