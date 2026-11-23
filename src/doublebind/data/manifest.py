"""Checksums for every raw file, so a rebuild can be proven identical (Day 1 'done when')."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Manifest:
    """JSON file stored inside the raw-data folder; paths are relative to that folder."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.root = self.path.parent
        self.entries: dict[str, dict] = {}
        if self.path.exists():
            self.entries = json.loads(self.path.read_text(encoding="utf-8"))

    def record(self, file: Path, source: str) -> None:
        file = Path(file)
        rel = file.relative_to(self.root).as_posix()
        self.entries[rel] = {
            "sha256": sha256_file(file),
            "bytes": file.stat().st_size,
            "source": source,
            "recorded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        self.save()

    def save(self) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.entries, indent=2, sort_keys=True), encoding="utf-8")

    def verify(self) -> list[tuple[str, str]]:
        """Return (relative path, 'ok' | 'missing' | 'changed') for every recorded file."""
        results = []
        for rel, info in sorted(self.entries.items()):
            f = self.root / rel
            if not f.exists():
                results.append((rel, "missing"))
            elif sha256_file(f) != info["sha256"]:
                results.append((rel, "changed"))
            else:
                results.append((rel, "ok"))
        return results
