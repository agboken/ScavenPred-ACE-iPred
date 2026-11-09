#!/usr/bin/env python3
"""Restore publicly archived ScavenPred/ACE-iPred research datasets.

No proprietary SAPIEN-Pep material is contained in these archives.
"""
from __future__ import annotations
import argparse
import base64
import gzip
import hashlib
import io
import tarfile
from pathlib import Path

EXPECTED_SHA256 = '5361826554c55de3f0f05b725f6c86eaa0d1c701a985728b28366fc12d6fb6bc'


def restore(directory: Path, destination: Path) -> None:
    bundle = directory / 'ScavenPred_ACEiPred_dataset_bundle.tar.gz.b64'
    test = directory / 'scavenpred_test_256_features.csv.gz.b64'
    data = base64.b64decode(bundle.read_text(encoding='utf8').strip(), validate=True)
    found = hashlib.sha256(data).hexdigest()
    if found != EXPECTED_SHA256:
        raise ValueError(f'Dataset checksum mismatch: {found}')
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        for member in archive.getmembers():
            if not member.isfile() or '/' in member.name or '\\' in member.name:
                raise ValueError(f'Unexpected archive entry: {member.name}')
        archive.extractall(destination, filter='data')
    raw = gzip.decompress(base64.b64decode(test.read_text(encoding='utf8').strip(), validate=True))
    (destination / 'scavenpred_test_88.csv').write_bytes(raw)
    print(f'Restored {len(list(destination.glob("*.csv")))} CSV datasets to {destination}')
    print(f'Verified compressed archive SHA256: {found}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1] / 'data')
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'restored')
    args = parser.parse_args()
    restore(args.source, args.out)
