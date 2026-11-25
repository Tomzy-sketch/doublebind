"""ESM-2 plausibility scores (FR-A3): log-likelihood ratio of every amino-acid substitution.

LLR(pos, mut) = log P(mut | context) - log P(wt | context). Positive means the language model
finds the change at least as plausible as the wild type. It measures tolerance, not resistance.

Two modes (Meier et al. 2021 compare them):
  masked : hide one position at a time and read its distribution (L forward passes, most accurate)
  wt     : one forward pass on the unmasked sequence (fast, slightly less accurate)

The core below uses numpy only and takes any `logits_fn`, so it is testable without model weights.
`load_hf_esm` wraps a Hugging Face ESM-2 checkpoint (needs torch and transformers).
"""
from __future__ import annotations

import math

import numpy as np

AA20 = "ACDEFGHIKLMNPQRSTVWY"
MAX_RESIDUES = 1022  # ESM-2 position limit


def position_logits(logits_fn, token_ids, aa_ids, mode="masked", mask_id=None, batch_size=8,
                    positions=None, progress=None) -> dict:
    """Return {position: array of 20 logits (ordered like AA20)}.

    token_ids : [<cls>, residue tokens..., <eos>], so residue number p sits at index p
    logits_fn : callable (int array [B, T]) -> float array [B, T, V]
    """
    base = np.asarray(token_ids, dtype=np.int64)
    n_res = len(base) - 2
    if n_res > MAX_RESIDUES:
        raise ValueError(f"{n_res} residues is longer than the ESM-2 limit of {MAX_RESIDUES}")
    ids20 = np.array([aa_ids[a] for a in AA20])
    todo = list(positions) if positions is not None else list(range(1, n_res + 1))
    rows: dict[int, np.ndarray] = {}
    if mode == "wt":
        logits = logits_fn(base[None, :])
        for pos in todo:
            rows[pos] = np.asarray(logits[0, pos, ids20], dtype=np.float64)
        return rows
    if mode != "masked":
        raise ValueError("mode must be 'masked' or 'wt'")
    if mask_id is None:
        raise ValueError("masked mode needs mask_id")
    for start in range(0, len(todo), batch_size):
        chunk = todo[start:start + batch_size]
        batch = np.tile(base, (len(chunk), 1))
        for b, pos in enumerate(chunk):
            batch[b, pos] = mask_id
        logits = logits_fn(batch)
        for b, pos in enumerate(chunk):
            rows[pos] = np.asarray(logits[b, pos, ids20], dtype=np.float64)
        if progress:
            progress(min(start + batch_size, len(todo)), len(todo))
    return rows


def _logsumexp(x: np.ndarray) -> float:
    m = float(np.max(x))
    return m + math.log(float(np.sum(np.exp(x - m))))


def llr_table(protein: str, rows: dict):
    """Turn per-position logits into (substitution rows, per-position stats).

    substitution rows : (pos, wt_aa, mut_aa, llr) for all 19 alternatives
    per-position stats: (pos, wt_aa, wt_logp, entropy) over the 20 standard amino acids
    """
    subs, stats = [], []
    for pos in sorted(rows):
        wt = protein[pos - 1]
        if wt not in AA20:
            continue  # unknown residue, nothing to compare against
        logit = rows[pos]
        logp = logit - _logsumexp(logit)
        wt_i = AA20.index(wt)
        entropy = float(-(np.exp(logp) * logp).sum())
        stats.append((pos, wt, float(logp[wt_i]), entropy))
        for j, aa in enumerate(AA20):
            if aa != wt:
                subs.append((pos, wt, aa, float(logit[j] - logit[wt_i])))
    return subs, stats


def load_hf_esm(model_name: str, device: str = "auto"):
    """Load an ESM-2 checkpoint. Returns (logits_fn, tokenizer, mask_id, device_name)."""
    try:
        import torch
        from transformers import AutoTokenizer, EsmForMaskedLM
    except ImportError as e:  # pragma: no cover - depends on the user's environment
        raise ImportError('ESM scoring needs torch and transformers: pip install -e ".[esm]"') from e

    dev = device if device != "auto" else ("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name).eval().to(dev)

    def forward(batch: np.ndarray) -> np.ndarray:
        try:
            with torch.inference_mode():
                out = model(input_ids=torch.as_tensor(batch, device=dev))
            return out.logits.float().cpu().numpy()
        except RuntimeError as e:  # out of memory: split the batch in two and retry
            if "out of memory" in str(e).lower() and len(batch) > 1:
                if dev == "cuda":
                    torch.cuda.empty_cache()
                half = len(batch) // 2
                return np.concatenate([forward(batch[:half]), forward(batch[half:])], axis=0)
            raise

    return forward, tokenizer, tokenizer.mask_token_id, dev


def score_protein(protein: str, model_name: str, mode: str = "masked", device: str = "auto",
                  batch_size: int | None = None, progress=None):
    """Full pipeline for one protein. Returns (substitution rows, per-position stats, device name)."""
    logits_fn, tokenizer, mask_id, dev = load_hf_esm(model_name, device)
    token_ids = tokenizer(protein)["input_ids"]
    aa_ids = dict(zip(AA20, tokenizer.convert_tokens_to_ids(list(AA20))))
    if batch_size is None:
        batch_size = 16 if dev == "cuda" else 4
    rows = position_logits(logits_fn, token_ids, aa_ids, mode=mode, mask_id=mask_id,
                           batch_size=batch_size, progress=progress)
    subs, stats = llr_table(protein, rows)
    return subs, stats, dev
