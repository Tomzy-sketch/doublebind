# DoubleBind: Product Requirements Document v2.0

Evolutionary trap and collateral-sensitivity engine for antimicrobial resistance (formerly Project EvoTrap). Draft dated 2026-11-20.

## 1. Document control

v2.0 keeps the v1.0 idea but cuts week-one scope to one pathogen and two protein targets, and fixes the technical, data-licensing and biosecurity gaps listed below.

| Field | Value |
| --- | --- |
| Product | DoubleBind (formerly Project EvoTrap) |
| Version | 2.0, supersedes v1.0 "MVP to Commercial API Specification" |
| Status | Draft for a 7-day MVP sprint (Day 1 = 2026-11-20) |
| Lead role | Principal ML Engineer / Technical Founder |
| Domain | Computational biology, AI drug discovery, antimicrobial resistance (AMR) |
| Week-1 deliverable | Working prototype: evidence-report notebook, web demo, API stub |
| Long-term deliverable | Cloud ML inference API plus interactive evolution sandbox |

**What changed from v1.0**

1. Renamed EvoTrap to DoubleBind (Section 2).
2. Added a 7-day MVP tier; the 24-week plan becomes Phases 1 to 4 (Sections 4, 14, 17).
3. The simulator starts as a population-genetics model; reinforcement learning becomes an optional adversarial explorer (Section 7).
4. Mutation supply is separated from selection: exact codon-level reachability replaces BLOSUM/PAM masks (Section 7).
5. The CSI > 2.0 rule is now uncertainty-aware: a 4× susceptibility change equals only about 0.85 kcal/mol, below the 1.2 kcal/mol error target (Sections 6, 9, 19).
6. Leakage policy rebuilt: clustering at 70% identity cannot separate strains of one species (Section 8).
7. Benchmarks corrected: SKEMPI 2.0 covers protein–protein binding, not drug–protein binding (Sections 8, 9).
8. Collateral-sensitivity scope clarified: structure models capture same-target effects only (Section 6).
9. Licensing review added, including the AlphaFold 3 validation gate (Section 8).
10. Responsible-use and biosecurity controls added (Section 12).
11. Equations restored: the reward and CSI formulas were blank in the v1.0 text received, so confirm the Section 19 definitions match your intent.
12. Added Part B: a 7-day build workflow and a zero-cost toolchain (Sections 17 and 18).

## 2. Executive summary and name

DoubleBind predicts how a bacterium will escape Drug A and which second drug catches it, so combination therapies can be designed before resistance spreads.

Resistance evolves faster than new antibiotics are approved, so discovery that reacts to resistance arrives late. DoubleBind simulates the likely resistance mutations under Drug A, estimates what each mutation costs the protein, and checks whether the mutant is more vulnerable to an existing Drug B (collateral sensitivity). The output is a ranked, explained list of drug pairs for a wet lab to test. It is not a treatment recommendation.

**Pitch:** Trap resistance before it evolves.

**Why this name.** The pathogen is in a double bind: the mutation that saves it from Drug A can expose it to Drug B. "Bind" also names the molecular binding energy the engine predicts. It is two plain words that are easy to say, spell and remember, and it fits a pharma slide deck.

| Option | Why it works | Watch-out |
| --- | --- | --- |
| DoubleBind (recommended) | Names the pathogen's dilemma and the binding energy the model predicts | "Double bind" is also a psychology term; check domain and trademark |
| Pincer | Two drugs closing on one pathogen | "Pincer ligand" is a common chemistry term, which muddies search |
| Achilles | Resistance creates a weak spot | A widely used name across industries |

A web search on 2026-11-20 found no antibiotic or AI drug-discovery product called DoubleBind. That is not a trademark check: before you publish, search the WIPO Global Brand Database, check domain names and the GitHub organization name. If you choose another name, find and replace "DoubleBind" in this document.

## 3. Problem and opportunity

Bacterial resistance directly caused an estimated 1.14 million deaths in 2021, and the Lancet GRAM study forecasts 1.91 million direct deaths a year by 2050 ([GBD 2021 AMR Collaborators, 2024](https://experts.umn.edu/en/publications/global-burden-of-bacterial-antimicrobial-resistance-19902021-a-sy/)).

Counting deaths where resistance played a role, the same study estimates 4.71 million in 2021 and 8.22 million in 2050. Among Gram-negative bacteria, carbapenem resistance grew faster than resistance to any other antibiotic class.

**The gap.** Antibiotic development is reactive: a drug is approved, bacteria evolve around it, and the pipeline starts again. Collateral sensitivity offers a different lever, because resistance to one drug can make bacteria more susceptible to another. Experimental maps of these relationships exist for E. coli, but they are slow to build and can vary between strains (Section 20). DoubleBind bets that structure-aware simulation can prioritize which pairs a lab should test first.

**Adjacent work.** AI antibiotic-discovery efforts such as the [MIT Antibiotics-AI Project](https://collinslab.mit.edu/antibiotics-ai-project) design new molecules. DoubleBind pulls a different lever: steering resistance evolution for drugs that already exist.

## 4. Goals, non-goals and scope tiers

The 7-day MVP proves one claim on one organism: DoubleBind can rank known resistance mutations and report whether a second drug looks like a trap, with honest uncertainty. Everything else waits for later tiers.

**Goals**

- G1. Rank single-nucleotide resistance mutations by likelihood of fixation under Drug A, scored against documented mutations.
- G2. Estimate the fitness cost of each path (protein stability and function) so costly dead ends are down-weighted.
- G3. For drugs that bind the same protein, estimate whether the resistant variant binds Drug B better than wild type, and rank pairs with uncertainty.
- G4. Produce an evidence report per pair that a wet-lab biologist can read in five minutes.

**Non-goals, all tiers**

- Clinical decision support, dosing or treatment advice. Outputs are research hypotheses only.
- De novo drug design. The engine ranks existing approved compounds.
- Predicting transmissibility, virulence or host range.
- Replacing wet-lab work. No pair is called validated without in vitro MIC data.

**Scope tiers**

| Tier | When | Scope | Exit test |
| --- | --- | --- | --- |
| MVP | Days 1–7 | E. coli K-12. GyrA for mutation ranking and a CSI demo; TEM-1 to validate the simulator. Pretrained models only. | Known-mutation recall reported, leakage tests pass, one evidence report |
| Beta | Weeks 2–12 | Adds parC, gyrB and CTX-M; custom surrogate on clade-split data; PPO explorer | Surrogate beats the pretrained baseline on held-out proteins |
| v1.0 | Weeks 13–24 | ESKAPE pathogens; metabolic-model penalty; API beta | Two design partners start in vitro tests |

**Pair outcome categories.** Every (Drug A, Drug B) pair is labeled collateral sensitivity, resistance-sparing, cross-resistance or inconclusive (definitions in Section 19). Same-target pairs often land in the last three, so an empty "trap" list is a valid scientific result, not a failure.

## 5. Users, use cases and acceptance stories

Four user types drive the requirements, and each acceptance story can be tested on the MVP unless it is marked Beta.

| User | Pain point | Acceptance story |
| --- | --- | --- |
| Computational biologist (pharma R&D) | Wet-lab attrition when screening drug pairs; no way to simulate multi-step evolution at scale | Given E. coli GyrA and ciprofloxacin, I get ranked Drug B candidates, each with a CSI interval and a category label (MVP) |
| Translational microbiologist | No view of the fitness trade-offs behind resistance mutations | Given a mutation such as S83L, I see the predicted stability change, the binding change per drug and a confidence label (MVP) |
| IP and portfolio strategist | Generic antibiotics lack patent protection | Given a ranked pair, I can export an evidence report with provenance that states the result is computational and unvalidated (MVP) |
| Wet-lab design partner | Needs protocol-ready candidates to test | Given a pair, I receive the mutation to engineer, the two drugs and a concentration range for an in vitro assay (Beta) |

**Primary use case (MVP).** An analyst selects a target and Drug A, runs the simulation, reviews the top resistance paths, and opens the evidence card for any Drug B that the engine labels collateral sensitivity or resistance-sparing.

## 6. Scientific foundations and assumptions

DoubleBind rests on six assumptions, and each has a known failure mode that the product must show to the user instead of hiding.

| # | Assumption | Known failure mode | Product response |
| --- | --- | --- | --- |
| A1 | Resistance arises from single-nucleotide missense changes in the drug target | Efflux, permeability, regulatory and gene-dosage routes are common | Label results "target-mediated only"; add efflux and metabolic terms in v1.0 |
| A2 | A binding-energy change converts to a fold change in IC50 (Section 19) | MIC also depends on uptake and efflux; for enzymes such as β-lactamases, catalytic efficiency matters more than binding | Use as a ranking signal; train an MIC-based proxy for enzymes in Beta |
| A3 | Collateral sensitivity comes from same-target structural change | Many documented cases come from physiology (membrane potential, efflux), not from Drug B's own target | Tier 1 (same-target, structural) is in scope; Tier 2 (physiology, metabolic models) is deferred; other pairs return "not assessable" |
| A4 | Evolution is an adaptive walk driven by selection coefficients | Clonal interference, recombination, plasmids, mutator strains and drug gradients | Start with the strong-selection weak-mutation model; add a population simulator in Beta |
| A5 | Mutation effects combine through a protein-stability threshold | Specific epistasis, such as TEM-1 M182T compensating for destabilizing changes | Use a stability-threshold fitness model and test it on the TEM-1 landscape |
| A6 | Predicted structures respond to single substitutions | Structure predictors often under-respond to point mutations | Start from experimental wild-type structures; treat predicted mutant structures as secondary evidence |

**Mutation supply is not selection.** v1.0 used BLOSUM/PAM-style masks to stop impossible jumps, but those matrices describe which substitutions evolution accepts, not which mutations the DNA can produce. v2.0 separates the two: the gene sequence defines exactly which amino-acid changes one nucleotide can reach (supply), and fitness decides which of them spread (selection).

**Honest expectations.** Same-target structural effects will find few true collateral-sensitivity pairs for a target like GyrA, because most gyrase mutations lower susceptibility to related drugs. Expect more cross-resistance and resistance-sparing calls than traps. Fluoroquinolone binding also depends on a metal-ion bridge that docking scores handle poorly, so report position-level recall (S83, D87) separately from substitution-level recall (S83L, D87N).

## 7. Functional requirements

Six modules (A to F) turn a target protein and Drug A into ranked, explained drug pairs; every requirement has a tier and a test.

### Module A: Data and mutation grammar

| ID | Requirement | Tier | Acceptance test |
| --- | --- | --- | --- |
| FR-A1 | Fetch the wild-type protein, its gene DNA and an experimental structure for a target | MVP | Gene translation matches the UniProt sequence exactly, or the run stops with an error |
| FR-A2 | Enumerate every single-nucleotide change per codon, keep missense variants, and weight each by a configurable transition/transversion ratio | MVP | Variant counts match the genetic-code table in unit tests |
| FR-A3 | Score each reachable variant with the ESM-2 log-likelihood ratio (zero-shot plausibility) and cache it | MVP | Table covers 100% of reachable variants for GyrA and TEM-1 |
| FR-A4 | Fine-tune ESM-2 on bacterial sequences | Beta | Must beat the zero-shot score on held-out clades, otherwise keep zero-shot |
| FR-A5 | Add insertions, deletions and efflux-upregulation actions | v1.0 | Each action has its own mutation-supply model and tests |

### Module B: Energetics surrogate

| ID | Requirement | Tier | Acceptance test |
| --- | --- | --- | --- |
| FR-B1 | Predict ΔΔG_fold for each variant with a pretrained model (ThermoMPNN) on the experimental structure | MVP | Sign-convention test passes; TEM-1 G238S and M182T directions are compared with the literature and reported |
| FR-B2 | Estimate a ΔΔG_bind proxy by re-scoring each drug in the wild-type and mutant pocket; optionally compare Boltz-2 affinity on a small subset | MVP | Spearman correlation against a Platinum/MdrDB subset is reported; output is labeled "proxy, unvalidated for mutations" |
| FR-B3 | Train an SE(3)-equivariant graph network on stability (Mega-scale) and binding (MdrDB, Platinum) data with deep-ensemble uncertainty | Beta | Beats FR-B1 and FR-B2 on held-out protein clusters; calibration error reported |
| FR-B4 | Re-score the top 0.1% terminal states asynchronously with open-licensed co-folding and docking (Boltz-2, AutoDock Vina) | Beta | Disagreements with the surrogate are flagged in the report |
| FR-B5 | Serve single-variant scores from cache at 15 ms or less on an A100 | v1.0 | Load test passes |

### Module C: Evolution simulator

| ID | Requirement | Tier | Acceptance test |
| --- | --- | --- | --- |
| FR-C1 | Strong-selection weak-mutation (SSWM) simulator over the reachable-variant graph, up to 5 steps, with drug concentration as input | MVP | Reproduces the qualitative path accessibility of the TEM-1 five-mutation landscape; rank correlation reported |
| FR-C2 | Fitness model: drug inhibition × fraction folded × metabolic cost (Section 19) | MVP | Unit tests: higher IC50 raises fitness when drug is present; fraction folded stays between 0 and 1 |
| FR-C3 | Gymnasium environment plus a PPO explorer that searches for high-risk paths the SSWM walk misses, keeping the +2.5 kcal/mol stability bound and uncertainty termination from v1.0 | Beta | Finds the top-5 SSWM paths; every extra path carries an uncertainty flag |
| FR-C4 | Population simulator (Wright–Fisher or Gillespie) with clonal interference and drug gradients | Beta | Compared with SSWM on the same landscape; differences documented |
| FR-C5 | Metabolic-model penalty from flux balance analysis | v1.0 | Penalized mutants stay viable in the genome-scale model |

### Module D: Collateral sensitivity engine

| ID | Requirement | Tier | Acceptance test |
| --- | --- | --- | --- |
| FR-D1 | Build the Drug B pool from approved compounds with annotated activity against the same target (ChEMBL); all other compounds return "not assessable" | MVP | Pool of 50 compounds or fewer, each with a source ID |
| FR-D2 | Compute CSI per drug and terminal state; aggregate the median, the 10th percentile and the share of states with CSI ≥ 2.0 | MVP | Matches hand calculation on 3 test cases |
| FR-D3 | Call "collateral sensitivity" only if P(CSI ≥ 2.0) ≥ 0.8 under surrogate uncertainty | MVP | Calls respect the rule; "inconclusive" is a valid output |
| FR-D4 | Label every pair as collateral sensitivity, resistance-sparing, cross-resistance or inconclusive (Section 19) | MVP | Every pair in the report carries exactly one label |
| FR-D5 | Trap score = path-probability-weighted share of terminal states that are confident traps; Beta adds an escape-risk factor for double-resistant mutants within 2 steps | MVP, Beta | Score reproducible from the seed |

### Module E: Explainability

| ID | Requirement | Tier | Acceptance test |
| --- | --- | --- | --- |
| FR-E1 | Evidence card per pair: top paths, mutated positions, ΔΔG values, CSI interval, label, data and model versions, warnings | MVP | A reviewer reads the card in 5 minutes (timed with one person) |
| FR-E2 | 3D view of the wild-type and mutant pocket with changed contacts highlighted (distance-based hydrogen-bond and clash heuristics) | MVP | Contact diff matches manual inspection on 2 examples |
| FR-E3 | Provenance bundle (input hashes, versions, seeds) exportable for IP filings | v1.0 | Re-running the bundle reproduces the report |

### Module F: API and interface

| ID | Requirement | Tier | Acceptance test |
| --- | --- | --- | --- |
| FR-F1 | Notebook, static results viewer and a FastAPI stub (score one mutation, run a simulation, list pairs) | MVP | A stranger can open the demo and re-run the notebook |
| FR-F2 | Async jobs with webhook callbacks, authentication and rate limits | Beta | Job status survives a restart |
| FR-F3 | GraphQL endpoint | v1.0 | Build only if a design partner asks |

## 8. Data, licensing and leakage control

Two licensing facts decide what you can ship: CARD data needs a paid license for commercial organizations, and AlphaFold 3 weights and outputs are for non-commercial use only. Both are fine for a student prototype, and both block a commercial API as v1.0 specified it.

### Data sources

| Source | Use | Terms and status |
| --- | --- | --- |
| UniProt, RCSB PDB | Wild-type sequences and experimental structures | Open (from memory; confirm on each site) |
| NCBI Pathogen Detection, BV-BRC | Genomes and MIC phenotypes | Public; check each record (from memory) |
| [CARD](https://card.mcmaster.ca/about) | Known resistance mutations | Free for academic, government and non-profit use; commercial organizations need a written license and fee; the ontology alone is CC BY 4.0 (verified 2026-11-20) |
| ChEMBL | Approved drugs and same-target activity | CC BY-SA 3.0, share-alike (from memory) |
| DrugBank | Not used in the MVP | Free access is non-commercial (from memory); use ChEMBL instead |
| Mega-scale stability data (Tsuboyama 2023) | Train and test ΔΔG_fold | Public release (from memory; confirm terms) |
| [MdrDB](https://www.nature.com/articles/s42004-023-00920-7) | Protein–ligand ΔΔG upon mutation: 100,537 samples, 240 proteins, 2,503 mutations, 440 drugs | Built partly from kinase and cancer-cell-line sets, so check bacterial coverage first; license to be confirmed |
| Platinum | Curated ΔΔG_bind upon mutation | Check terms (from memory) |
| SKEMPI 2.0 | Protein–protein binding only | Removed from the drug-binding benchmark |

### Models and tools

| Tool | Decision | Why |
| --- | --- | --- |
| [AlphaFold 3](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md) | Excluded from the validation gate | Weights and outputs are limited to non-commercial organizations, and access is at Google's discretion (verified 2026-11-20) |
| [Boltz-2](https://jclinic.mit.edu/boltz-2-towards-accurate-and-efficient-binding-affinity-prediction/) | Used for validation and affinity checks | MIT license for code and weights, commercial use allowed (verified 2026-11-20) |
| [ThermoMPNN](https://www.tamarind.bio/tools/thermompnn) | Used for ΔΔG_fold in the MVP | Code under MIT (read on a third-party page; confirm in the repository) |
| ESM-2 | Used for plausibility scores | MIT-licensed weights (from memory; confirm on Hugging Face) |
| AutoDock Vina | Used for docking and re-scoring | Open source, Apache 2.0 (from memory) |

Boltz-2 is trained mainly to rank different ligands against one protein. Its response to point mutations in the protein is not established, so treat its output as an unvalidated signal until you benchmark it.

### Leakage-control policy

v1.0 required clustering at 70% sequence identity. That cannot separate strains of one species, because their proteins are nearly identical and every strain would land in one cluster. v2.0 splits at several levels:

1. **Protein level (energetics surrogate).** Cluster proteins with MMseqs2 and keep whole clusters on one side: 30% identity or lower for the strict unseen-protein test, 70% for the related-family test. A UniProt entry or PDB complex never appears on both sides, and all mutations at one position stay together.
2. **Strain level (phenotype validation).** Split by phylogroup or sequence type from core-genome distance, not by protein identity.
3. **Drug level.** For binding data, also hold out whole drug scaffolds in one test split.
4. **Automated checks.** CI tests assert zero cluster overlap and compare every model with shuffled-label and position-only baselines.
5. **Sign and unit tests at ingestion.** Mega-scale and some tools may use the opposite ΔΔG sign (positive = stabilizing). Check each source's definition, convert once, then test with known cases: a proline inside a helix must read as destabilizing.
6. **Versioning.** Every dataset gets a DVC version and a dataset card (source, date, license, filters).

## 9. Evaluation plan and KPIs

The MVP reports honest baselines and does not gate on v1.0 numbers; the Beta and v1.0 targets below are proposals to recalibrate after the Day 4 baselines.

| Metric | MVP (week 1) | Beta | v1.0 | Benchmark and method |
| --- | --- | --- | --- | --- |
| Resistance-mutation recall | Report position-level recall for GyrA S83 and D87 in the top 10 positions, plus substitution-level recall | Recall@5 of 60% or more on held-out strains | Recall@5 of 75% or more | Documented mutations from CARD and the literature (E. coli GyrA, ParC, TEM, CTX-M); baselines: random, proximity-only, ESM-only |
| ΔΔG_fold error | Report RMSE and Spearman on a held-out protein subset | RMSE of 1.2 kcal/mol or less on protein-cluster holdout | 1.2 kcal/mol or less | Mega-scale held-out domains, S669 |
| ΔΔG_bind error | Report Spearman on a Platinum/MdrDB subset | RMSE of 1.5 kcal/mol or less | 1.2 kcal/mol or less (stretch) | Platinum and MdrDB split by protein and drug scaffold; SKEMPI 2.0 is no longer used |
| Uncertainty calibration | Report coverage of 80% intervals | Coverage between 70% and 90% | Coverage between 75% and 85% | Held-out clusters |
| Path-ordering fidelity | Report rank correlation with the TEM-1 landscape | 0.6 or more | 0.7 or more | Five-mutation TEM-1 data (Weinreich 2006) |
| Collateral-sensitivity precision@3 | Report only, plus the share of published pairs the structural module can assess | 50% or more | 65% or more | Published E. coli collateral-sensitivity networks, used for validation only |
| Simulation latency | 100 ms per step or less from cache on a notebook GPU | 50 ms | 15 ms on an A100 | Timing harness |
| Job turnaround | 10,000 trajectories in under 10 minutes | Same on a free GPU | 100,000-generation rollout plus a 500-compound screen in under 10 minutes | End to end |

**Evaluation rules**

- Every metric ships with three baselines (random, position-only, ESM-only) and a shuffled-label control.
- Test sets are touched once; tuning uses validation splits only.
- Failure analysis lists the five worst mispredictions per metric.
- A missed target is reported as missed. Negative and inconclusive results stay in the report.
- Every number regenerates from a seed and a data version.

## 10. Architecture and technology stack

The heavy models run once and write a cached score table; the simulator and ranking engine then read that table, which is why the MVP fits on free notebooks.

*(Diagram: architecture overview. Public data, A mutations, B energetics, then a score cache that C simulator, D CSI engine, E evidence and F API read from; a Beta validation gate re-checks the top 0.1% of terminal states. See the PRD document for the picture.)*

Data flows left to right in the top row and drops into the cache; the second row reads from it, and the dashed gate is a Beta add-on that re-checks the top 0.1% of terminal states.

**Technology by tier**

| Layer | MVP (free) | Production (v1.0) |
| --- | --- | --- |
| Orchestration | Makefile and notebooks | Airflow or Dagster |
| Compute | Kaggle or Colab GPU, laptop CPU | Kubernetes on AWS EKS or GCP GKE, with Ray RLlib workers |
| Storage | DVC with a Google Drive remote, Parquet files | Object storage with versioned feature graphs |
| Models | Pretrained ESM-2, ThermoMPNN and Vina; Boltz-2 on a subset | Custom equivariant-network ensemble plus the Boltz-2 validation gate |
| Serving | Static viewer plus a FastAPI stub | FastAPI REST endpoints, async jobs with webhooks |
| Viewer | py3Dmol or 3Dmol.js | Mol* dashboard |
| Tracking | MLflow on a laptop or CSV logs | Central tracking and audit logs |

## 11. Non-functional requirements

Week one runs at zero cost on notebooks and static hosting; production infrastructure waits until a design partner commits.

| Area | MVP (week 1) | v1.0 |
| --- | --- | --- |
| Cost | $0, with no payment card needed for the tools in Section 18 | Paid cloud once a partner funds it |
| Compute | Kaggle or Colab GPU for batch jobs; laptop CPU for simulation | Containerized services on Kubernetes (EKS or GKE) with Ray RLlib workers |
| Reproducibility | Pinned environment file, fixed seeds, DVC data versions; one command regenerates all results | Signed container images; a run manifest stored with each report |
| Availability | Static site always on; the API stub may sleep | Targets set with the first partner |
| Security | Public data only; no secrets in the repository | Tenant-isolated VPCs, TLS 1.3 in transit, AES-256 at rest, so partner compounds never enter shared model weights |
| Explainability | Evidence card per pair (Section 7, Module E) | Structural diff with hydrogen-bond and steric-clash notes for every recommended pair |
| Compliance | License table kept current (Section 8) | Commercial licenses signed before the API beta |
| Observability | Metrics and seeds logged per run (MLflow or a CSV file) | Central logs, tracing and cost dashboards |

Performance targets live in Section 9.

## 12. Responsible use and biosecurity

DoubleBind predicts how pathogens evolve resistance, which is dual-use knowledge, so scope and access controls are product requirements and not afterthoughts.

- **Scope.** Antimicrobial-resistance research only. No feature may increase transmissibility, virulence, host range or immune escape.
- **Organisms.** The MVP covers E. coli K-12; later tiers add ESKAPE pathogens. Select agents and other regulated pathogens are excluded from inputs.
- **Outputs.** Resistance paths are always paired with the countermeasure analysis (Drug B). Public demos show only curated, already-published targets and mutations; novel-mutation predictions go to vetted partners.
- **Access.** The partner API requires identity verification, an institutional affiliation, rate limits and audit logs.
- **Claims.** Every report states "research hypothesis, not clinical guidance" and lists its unvalidated assumptions.
- **Review.** Before any public release, ask an institutional biosafety or research-ethics contact to review the dual-use risk, and publish model and data cards with the release.
- **Wet-lab handoff.** Partners test in vitro under their own biosafety approval. DoubleBind does not supply culture or strain-engineering protocols.

## 13. Risk register

The biggest risks are scientific, not engineering: the thresholds are tighter than the models' error, and much collateral sensitivity comes from mechanisms the structural model cannot see.

| Risk | Likelihood | Impact | Mitigation |
| --- | --- | --- | --- |
| Threshold smaller than model error: a 4× change is about 0.85 kcal/mol, below the 1.2 kcal/mol error target | High | High | Probability-based calls (Section 19); report "inconclusive" honestly |
| Mechanism coverage: much collateral sensitivity comes from physiology, not the drug target | High | Medium | Tier 1 label on every result; efflux and metabolic modules in later tiers |
| Sim-to-real gap: predicted mutants fail in living cells because of epistasis and network effects | High | High | Metabolic-model penalty (v1.0); label outputs as hypotheses; partner in vitro validation |
| Surrogate exploitation: the explorer hunts for model errors (reward hacking) | Medium | High | Hard stability bound of +2.5 kcal/mol, ensemble uncertainty termination, validation-gate re-scoring |
| Sparse collateral-sensitivity ground truth | High | Medium | Frame as a biophysical zero-shot task; use published matrices only for validation, never for training |
| Data leakage inflates metrics | Medium | High | Multi-level split policy and CI leakage tests (Section 8) |
| Licensing blocks commercialization (CARD, AlphaFold 3, others) | High | High | License table in Section 8; prefer MIT-licensed models; obtain a CARD license or curate from the literature before an API beta |
| One-week scope creep | High | High | Cut-lines and scope freeze at the end of Day 2 (Section 17) |
| Free-tier limits stop jobs (Kaggle weekly quota, Colab unpublished quota) | Medium | Medium | Cache results, checkpoint every run, keep jobs under the 12-hour session limit, run CPU stages locally |
| Dual-use misuse | Low | High | Controls in Section 12 |
| Overclaiming | Medium | High | Evidence cards show uncertainty and baselines; nothing is "validated" without MIC data |

## 14. Roadmap and exit criteria

The 24-week plan stays, but it now starts with a one-week Phase 0, and every phase closes with a testable gate instead of a date.

*(Diagram: five phases, each closed by a gate. Gate 0: recall reported and leakage tests pass. Gate 1: ΔΔG_fold RMSE at most 1.2 on held-out clusters. Gate 2: reproduces gyrA and TEM-1 paths. Gate 3: ranked pairs with evidence cards. Gate 4: two partners start in vitro tests.)*

Phase 0 is this week's build. A failed gate holds the next phase back, and the week ranges are planning estimates.

| Phase | Core deliverables |
| --- | --- |
| 0 (week 1) | Variant score table for GyrA and TEM-1, SSWM simulator, one evidence card, static viewer, honest metrics report |
| 1 (weeks 2 to 6) | Data pipeline with clade-aware splits, custom equivariant surrogate for ΔΔG_fold and ΔΔG_bind, calibration |
| 2 (weeks 7 to 12) | Gymnasium environment, PPO explorer, population simulator, ciprofloxacin and cefotaxime resistance trajectories in E. coli |
| 3 (weeks 13 to 18) | Collateral engine across the Drug B pool, evidence reports with structural diffs, Mol* dashboard, FastAPI backend |
| 4 (weeks 19 to 24) | ESKAPE pathogens, metabolic-model penalty, API beta with tenant isolation, design-partner onboarding |

## 15. IP and commercial considerations

Computational predictions alone rarely support a patent; commercial value comes from partner-validated data. This is general information, not legal advice, so talk to a patent attorney before filing anything.

- **v1.0 over-claimed.** It said the platform enables composition-of-matter patents on combinations of off-patent drugs. Mixtures of known compounds are generally hard to claim as compositions, and method-of-use or regimen claims usually need experimental evidence of an unexpected effect.
- **Build the evidence trail now.** Dated lab notebooks, in vitro MIC data from partners and the provenance bundle (FR-E3) are what a filing needs.
- **Licenses gate revenue.** A commercial product cannot ship CARD data, AlphaFold 3 outputs or DrugBank data without licenses (Section 8).
- **Early revenue paths.** Academic collaborations, grants and paid pilots with biotech design partners come before any API subscription.

## 16. Open decisions

Seven decisions change the plan; each row shows the default this document assumes, so edit the row if your answer differs.

| Decision | Default assumed here | Why it matters |
| --- | --- | --- |
| Start and finish | Day 1 is Friday 2026-11-20; v0.1 ships on Thursday 2026-11-26 | Sets the daily plan in Section 17 |
| Team | One builder, with one optional helper | A second person can take the energetics days (3 and 4) while you build the simulator |
| Purpose of the week | A working prototype with an honest report | A competition, a course project and a lab pitch need different demos |
| Final name | DoubleBind | Find and replace if you change it; run trademark and domain checks |
| First targets | E. coli GyrA and TEM-1 | Smaller scope; ESKAPE pathogens come later |
| Negative results | Publish them | Honest "inconclusive" calls are expected for same-target pairs |
| Commercial intent | After partner validation | Decides when CARD and other licenses must be bought |

## 17. Part B: 7-day build workflow

Build in this order: data, mutation grammar, stability, binding proxy, simulator, collateral engine, then packaging. Each day ends with a saved result and a test, so a failed day never silently breaks the next one.

*(Diagram: seven days, one saved result each, with checkpoints at the end of Day 2 (scope freeze), Day 4 (binding go/no-go) and Day 7 (release v0.1).)*

Scope freezes at the end of Day 2, the binding proxy gets a go/no-go at the end of Day 4, and v0.1 ships on Day 7. Days 3 and 4 are the riskiest because of tool installs and docking, so start them early and do not polish Days 1 and 2.

### Day 1, Fri Nov 20: set up and collect data

1. Create a GitHub repository called doublebind and a Kaggle account; verify your phone so Kaggle unlocks GPUs.
2. Copy the layout from Section 18, pin package versions in environment.yml, and add a Makefile with one target per day.
3. Download the wild-type data: UniProt records for E. coli GyrA, GyrB and TEM-1 (confirm the accessions P0AES4, P0AES6 and P62593 on UniProt), the gene DNA from NCBI, and experimental structures from RCSB PDB (TEM-1 and a gyrase–DNA–quinolone complex; if no E. coli complex exists, use the closest homolog and say so in the report).
4. Download ChEMBL's approved drugs with their activity against gyrase and TEM-1, and CARD's gyrA and TEM mutation lists (academic use only for now).
5. Start DVC with a Google Drive remote and commit the raw data with checksums.
6. Write the sign-convention test before any model code.

**Done when:** `make data` rebuilds every raw file with matching checksums and the three accessions are confirmed.

### Day 2, Sat Nov 21: mutation grammar and plausibility scores

1. In codons.py, read the gene DNA, list the 9 single-nucleotide changes per codon, translate them and keep the missense variants. Make the transition/transversion weight a parameter.
2. Test on a small hand-made gene: counts must match the genetic-code table.
3. On Kaggle, run ESM-2 (start with the 150M model, then 650M) to get the log-likelihood ratio for every reachable variant, and save it as Parquet.
4. Mark residues within about 8 Å of the bound drug in the structure (Biopython).
5. Compute three GyrA baselines: random, ligand-proximity only, ESM-only. Record whether S83 and D87 appear in the top 10 positions.

**Done when:** the variant table exists for GyrA and TEM-1, baseline numbers are saved in a notebook, and you freeze scope.

### Day 3, Sun Nov 22: stability scores

1. Install ThermoMPNN on Kaggle by following its README, then check its sign on three known mutations.
2. Score every reachable variant on the experimental structure and add a ΔΔG_fold column to the table.
3. Compare with the literature for TEM-1, where G238S is reported as destabilizing and M182T as stabilizing. Record whether your numbers agree, even if they do not.
4. Save the table and the evaluation notebook, and commit.

**Done when:** ΔΔG_fold exists for at least 95% of reachable variants and the sign test passes.

### Day 4, Mon Nov 23: binding proxy

1. Build the Drug B pool: approved drugs with ChEMBL activity against the same target (50 or fewer).
2. Prepare ligands with RDKit and Meeko, dock each into the wild-type pocket with AutoDock Vina and keep the best pose.
3. For each variant near the pocket, build the mutant side chain with PDBFixer, run a short OpenMM minimization, and re-score the same pose. The difference from wild type is your ΔΔG_bind proxy.
4. Run Boltz-2 on 10 to 20 variants as a cross-check.
5. Compare against a Platinum or MdrDB subset and record the Spearman correlation.

**Go/no-go:** if the correlation is near zero or docking fails for most drugs, drop the proxy, mark CSI as "not assessed", and ship mutation ranking, stability and the simulator.

### Day 5, Tue Nov 24: evolution simulator

1. Write fitness.py from the Section 19 formula, with unit tests for the monotonic cases.
2. Write sswm.py: nodes are genotypes with up to 5 mutations, edges are single changes, transition probabilities follow Section 19, and a fixed seed drives 10,000 trajectories.
3. Test the path logic on measured data first: rebuild the TEM-1 five-mutation landscape from the published MIC values and compare simulated path frequencies with the paper's pattern. Then swap in your predicted fitness and see what degrades.
4. Run GyrA under low, medium and high ciprofloxacin concentrations and save the top paths.

**Done when:** the path-ordering score is reported and a 10,000-trajectory run finishes in under 10 minutes.

### Day 6, Wed Nov 25: collateral engine and evidence card

1. For each terminal mutant from Day 5, compute CSI for every drug in the pool, then the median, the 10th percentile and the share above 2.0.
2. Apply the decision rule and label every pair.
3. Write evidence_card.py: top three paths, ΔΔG values, CSI interval, label, warnings and versions.
4. Add a 3D view of the wild-type and mutant pocket with changed contacts highlighted (py3Dmol).
5. Generate cards for GyrA with ciprofloxacin as Drug A.

**Done when:** one command produces the evidence cards from the cached table.

### Day 7, Thu Nov 26: package and ship

1. Build a static viewer (HTML, Plotly, 3Dmol.js) that reads your Parquet and JSON outputs, and publish it on GitHub Pages.
2. Add the FastAPI stub (score one mutation, run a short simulation, list pairs) and run it inside a notebook for the demo.
3. Write the README, MODEL_CARD, DATA_CARD and a one-page limitations list copied from Sections 6 and 13.
4. Record a 3-minute demo: the problem, one evidence card, and the honest limitations.
5. Tag v0.1.

**Done when:** a stranger can open the demo link, read one evidence card and re-run the notebook.

### If you fall behind

Cut in this order: (1) the Boltz-2 cross-check, (2) the binding proxy, with CSI marked not assessed, (3) the web viewer, replaced by a notebook report. Never cut the sign and leakage tests, the baselines or the limitations page.

### Daily rhythm

- **Start:** pull, run the tests, and write one line with today's goal in the README log.
- **During:** commit after each step, and save every Kaggle output to a dataset or Drive before the session ends.
- **End:** update the metrics table, push, and write three lines: what worked, what failed, what is next.

## 18. Part B: Free toolchain and repository layout

Everything in the MVP runs on free tools that need no payment card, but two free tiers changed recently, so read the notes below before you plan deployment.

### Compute and hosting (checked 2026-11-20)

| Service | What the free tier gives you | Catch |
| --- | --- | --- |
| [Kaggle Notebooks](https://aimultiple.com/free-cloud-gpu) | Two T4 GPUs (16 GB each), typically 30 GPU hours a week, 12-hour sessions, 20 GB saved storage | Phone verification needed; the P100 option was retired on 2026-09-15; queues at busy times |
| Google Colab (free) | Sessions up to 12 hours; Google assigns the GPU | No published reset schedule; access varies with demand; runtime files vanish, so save to Drive |
| Lightning AI | Five starter credits; no card needed to sign up | Phone verification needed; credits are promotional and do not accumulate |
| [Hugging Face Spaces](https://huggingface.co/docs/hub/en/spaces-overview) | Static Spaces are free for everyone | Gradio and Docker Spaces now need a paid plan, except up to two Gradio Spaces on ZeroGPU for eligible free accounts |
| Modal | $30 of monthly credits | Requires a payment method, so skip it |
| Streamlit Community Cloud, GitHub Pages | Free app and static-site hosting | Shared resources; confirm current limits (from memory) |

**A deployment plan that stays free.** Run every heavy step on Kaggle or Colab, save the results as Parquet and JSON, and publish a static viewer (GitHub Pages or a Hugging Face Static Space) that reads those files. Use a notebook with a tunnel only when you record the live demo.

### Software stack

| Need | Tool | Notes |
| --- | --- | --- |
| Code and teamwork | Git, GitHub (repositories, Actions, Projects), VS Code | Free personal accounts |
| Protein language model | ESM-2 through Hugging Face Transformers | Start with the 35M and 150M models, then try 650M |
| Stability | ThermoMPNN | MIT code |
| Docking and affinity | AutoDock Vina with Meeko and RDKit; Boltz-2 on a small subset | Open source |
| Structures | Biopython, PDBFixer with OpenMM, py3Dmol | Open source |
| Sequence clustering | MMseqs2 | GPL-3.0; call it as a command-line tool |
| Deep learning | PyTorch, PyTorch Geometric, e3nn (Beta) | Open source |
| Reinforcement learning (Beta) | Gymnasium, Stable-Baselines3 | Open source |
| Data | Pandas or Polars, PyArrow (Parquet), DVC with a Google Drive remote | Free tiers |
| Experiment tracking | MLflow on your machine, or a Weights & Biases personal plan | Free |
| Pipelines | Makefile in the MVP; Prefect or Dagster in Beta; Airflow only at v1.0 | Open source |
| API and viewer | FastAPI with Uvicorn, Streamlit, 3Dmol.js; Mol* in Beta | Open source |
| Writing and diagrams | MkDocs Material, Mermaid, Zotero, Overleaf free tier | Free |

### Repository layout

```text
doublebind/
  data/            raw/ and processed/ (DVC-tracked)
  src/doublebind/
    data/          fetch.py, codons.py
    scoring/       esm_llr.py, stability.py, binding.py
    sim/           fitness.py, sswm.py, env.py (Beta)
    csi/           csi.py, categories.py, trap_score.py
    report/        evidence_card.py, contacts.py
    api/           main.py
  app/             static viewer or Streamlit app
  notebooks/       01_data.ipynb to 06_report.ipynb
  tests/           test_codons.py, test_signs.py, test_leakage.py, test_csi.py
  Makefile  pyproject.toml  environment.yml
  README.md  MODEL_CARD.md  DATA_CARD.md
```

## 19. Appendix: formulas, conventions, glossary and API sketch

All equations use kcal/mol and 37 °C (310 K), where RT = 0.616 kcal/mol; at 25 °C, RT = 0.593 kcal/mol.

### Sign conventions

- ΔΔG_fold = ΔG_fold(mutant) − ΔG_fold(wild type). Positive means destabilizing.
- ΔΔG_bind = ΔG_bind(mutant) − ΔG_bind(wild type). Positive means weaker binding, which means more resistance.
- A positive CSI means the mutant is more susceptible to the drug.
- Tools differ. ThermoMPNN's documented output is negative for stabilizing changes, which matches this PRD, and Boltz-2 reports a log(IC50) in µM. Check each tool's convention and unit-test it.

### Binding energy to IC50

```formula
IC50(S) / IC50(S0)  ≈  exp( ΔΔG_bind / RT )
```

This treats IC50 as proportional to the dissociation constant. Use it as a ranking signal, not an absolute prediction (assumption A2).

### Collateral Sensitivity Index

```formula
CSI_i(S) = log2 [ IC50,i(S0) / IC50,i(S) ] = − ΔΔG_bind,i / ( RT · ln 2 )
```

RT ln 2 is 0.427 kcal/mol, so a CSI of 2.0 (a 4× gain in susceptibility) needs a ΔΔG_bind of −0.85 kcal/mol or lower. The surrogate returns a mean μ and a standard deviation σ for ΔΔG_bind, and the call rule is:

```formula
P( CSI_i ≥ 2 ) = Φ( ( −2 · RT · ln 2 − μ ) / σ )  ≥  0.8
```

Worked example: with σ = 1.2 kcal/mol, the predicted mean must be about −1.9 kcal/mol or lower, a 20-fold predicted gain, before the engine calls a trap.

### Pair categories

Labels are assigned on terminal states that are already resistant to Drug A.

| Label | Rule |
| --- | --- |
| Collateral sensitivity | CSI_B of 2.0 or more, with P(CSI_B ≥ 2.0) of 0.8 or more |
| Resistance-sparing | CSI_B between −1 and 1: less than a 2× change either way |
| Cross-resistance | CSI_B of −1 or less |
| Inconclusive | Everything else, including CSI_B from 1 to 2 and low-confidence calls |

### Fitness model

```formula
W(S) = 1 / ( 1 + ( D / IC50(S) )^h )  ×  1 / ( 1 + exp( ΔG_fold(S) / RT ) )  ×  exp( −c_met(S) )
```

The first factor is drug inhibition at concentration D with Hill coefficient h (default 1). The second is the fraction of protein that is folded: ΔG_fold(S) is the wild-type value plus the sum of ΔΔG_fold over the mutations (an additive approximation), and ΔG_fold is negative for a stable protein. Use −5 kcal/mol as a placeholder wild-type value until you find a measured one, and run a sensitivity check around it. The third factor is a metabolic cost, constant in the MVP.

### Evolution step

```formula
s_ij = W(S_j) / W(S_i) − 1
π_ij = ( 1 − exp(−2 s_ij) ) / ( 1 − exp(−2 N_e s_ij) )
P(i → j) = μ_ij · π_ij / Σ_k ( μ_ik · π_ik )
```

μ_ij is the mutation supply (per-site rate times a transition/transversion weight), N_e is the effective population size, and π_ij is the chance that a new mutant takes over the population; use 1/N_e when s_ij is 0. This strong-selection weak-mutation approximation ignores clonal interference, which is why Beta adds a population simulator.

### Explorer reward (Beta)

```formula
R_t = log W(S_t+1) − log W(S_t) − λ · U(S_t+1)
```

U is the ensemble uncertainty. End an episode when the cumulative ΔΔG_fold passes +2.5 kcal/mol or U passes a set threshold.

### Trap score

```formula
T(D_A, D_B) = Σ over S_T of  w(S_T) · 1[ confident collateral sensitivity to D_B at S_T ]
```

w(S_T) is the probability, from the simulator, of reaching terminal state S_T under Drug A.

### Glossary

| Term | Plain meaning |
| --- | --- |
| AMR | Antimicrobial resistance: bacteria surviving drugs that used to kill them |
| MIC and IC50 | MIC is the lowest drug concentration that stops visible growth; IC50 is the concentration that halves an effect. Lower means more susceptible |
| ΔΔG | The energy change a mutation causes, in kcal/mol; the model's unit for stability and binding |
| Collateral sensitivity | Resistance to Drug A that makes the same bacteria more vulnerable to Drug B |
| CSI | Collateral Sensitivity Index: log base 2 of how much more susceptible the mutant is to a drug |
| SSWM | Strong selection, weak mutation: a simple model where one mutation at a time takes over |
| QRDR | Quinolone resistance-determining region, where fluoroquinolone resistance mutations cluster in gyrase genes |
| ESKAPE | Six dangerous hospital pathogens: Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae, Acinetobacter baumannii, Pseudomonas aeruginosa and Enterobacter species |
| Epistasis | When the effect of one mutation depends on which other mutations are present |
| LLR | Log-likelihood ratio from a protein language model: how plausible an amino-acid change looks |
| Clade split | Dividing data by evolutionary group so test cases are truly unseen |

### API sketch (illustrative values)

```json
POST /v1/score-mutation
{ "target": "P0AES4", "mutation": "S83L", "ligand_smiles": "..." }
-> { "ddg_fold": 0.4, "ddg_bind": 1.1, "sigma": 0.9, "model_version": "..." }

POST /v1/simulate
{ "target": "P0AES4", "drug_a": "ciprofloxacin", "n_trajectories": 10000, "seed": 7 }
-> { "job_id": "..." }

GET /v1/pairs/{job_id}
-> [ { "drug_b": "...", "label": "inconclusive", "csi_median": 0.6, "csi_p10": -0.8, "trap_score": 0.0 } ]
```

## 20. References

The first list holds the pages behind the numbers and license terms in this document, read on 2026-11-20 (several through search extracts, so re-read each license page before you sign or ship anything). The second list is background literature written from memory, so confirm each entry before you cite it.

**Sources read**

- [GBD 2021 Antimicrobial Resistance Collaborators, The Lancet 404(10459), 2024](https://experts.umn.edu/en/publications/global-burden-of-bacterial-antimicrobial-resistance-19902021-a-sy/)
- [CARD terms of use](https://card.mcmaster.ca/about) and [CARD licensing page](https://research.mcmaster.ca/?p=9267)
- [AlphaFold 3 model parameters terms of use](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)
- [Boltz-2 announcement, MIT Jameel Clinic](https://jclinic.mit.edu/boltz-2-towards-accurate-and-efficient-binding-affinity-prediction/)
- [ThermoMPNN summary and license](https://www.tamarind.bio/tools/thermompnn) and the [ThermoMPNN preprint](https://www.biorxiv.org/content/10.1101/2023.07.27.550881.full.pdf)
- [MdrDB, Communications Chemistry 2023](https://www.nature.com/articles/s42004-023-00920-7)
- [Hugging Face Spaces overview](https://huggingface.co/docs/hub/en/spaces-overview)
- [Free cloud GPUs: models, limits, access requirements, AIMultiple, updated 2026-11-06](https://aimultiple.com/free-cloud-gpu)
- [MIT Antibiotics-AI Project](https://collinslab.mit.edu/antibiotics-ai-project)

**Background literature (from memory)**

- Szybalski W, Bryson V. Cross resistance of Escherichia coli to fifteen antibiotics. J Bacteriol, 1952.
- Imamovic L, Sommer MOA. Use of collateral sensitivity networks to design drug cycling protocols. Sci Transl Med, 2013.
- Lázár V et al. Bacterial evolution of antibiotic hypersensitivity. Mol Syst Biol, 2013.
- Podnecky NL et al. Conserved collateral antibiotic susceptibility networks in diverse clinical strains of Escherichia coli. Nat Commun, 2018.
- Mira PM et al. Rational design of antibiotic treatment plans (TEM-1 alleles). PLoS ONE, 2015.
- Weinreich DM et al. Darwinian evolution can follow only very few mutational paths to fitter proteins. Science, 2006.
- Wang X, Minasov G, Shoichet BK. Evolution of an antibiotic resistance enzyme constrained by stability and activity trade-offs. J Mol Biol, 2002.
- Bloom JD et al. Protein stability promotes evolvability. PNAS, 2006.
- DePristo MA, Weinreich DM, Hartl DL. Missense meanderings in sequence space. Nat Rev Genet, 2005.
- McCandlish DM, Stoltzfus A. Modeling evolution using the probability of fixation. Q Rev Biol, 2014.
- Tsuboyama K et al. Mega-scale experimental analysis of protein folding stability in biology and design. Nature, 2023.
- Pires DEV, Ascher DB. Platinum: a database of experimentally measured effects of mutations on structurally defined protein-ligand complexes. Nucleic Acids Res, 2015.
- Lin Z et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science, 2023.
- Eberhardt J et al. AutoDock Vina 1.2.0. J Chem Inf Model, 2021.
