#!/usr/bin/env python3
"""Rebuild the public CSV datasets from size-limited GitHub archive parts."""
from __future__ import annotations
import argparse
import base64
import hashlib
from pathlib import Path
import tarfile
import io

EXPECTED_SHA256 = "5361826554c55de3f0f05b725f6c86eaa0d1c701a985728b28366fc12d6fb6bc"
PARTS = 23

def restore(root: Path, dest: Path):
    source = root / "data" / "archive_parts"
    files = [source / f"part_{i:03}.b64" for i in range(1, PARTS+1)]
    missing = [f.name for f in files if not f.exists()]
    if missing:
        raise FileNotFoundError(f"Missing dataset parts: {missing}")
    text = "".join(f.read_text(encoding="ascii").strip() for f in files)
    data = base64.b64decode(text, validate=True)
    sha = hashlib.sha256(data).hexdigest()
    if sha != EXPECTED_SHA256:
        raise ValueError(f"Checksum mismatch: {sha}")
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
        for member in tf.getmembers():
            if not member.isfile() or "/" in member.name or ".." in Path(member.name).parts:
                raise ValueError(f"Unsafe archive member: {member.name}")
        tf.extractall(dest, filter="data")
    print(f"Restored {len(list(dest.glob('*.csv')))} CSVs in {dest}")
    print("SHA256:", sha)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    restore(a.root, a.out or (a.root / "data" / "csv"))
