#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math, os, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
EXPECTED = ROOT / 'outputs_expected'
PHI = (1.0 + math.sqrt(5.0)) / 2.0
A1 = 2.0 * math.pi * PHI * math.log10(math.e) / 3.0
B_MASS = -4.7790

def read_csv_dict(path: Path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def write_json(obj, path: Path | None = None):
    text = json.dumps(obj, indent=2, ensure_ascii=False)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding='utf-8')
    print(text)

def rms(vals):
    vals = list(vals)
    return math.sqrt(sum(v*v for v in vals)/len(vals))

def parse_word(s: str):
    s = s.strip().replace('[','').replace(']','')
    if not s:
        return []
    return [int(x.strip()) for x in s.split(',') if x.strip()]

def inverse_braid_word(word):
    return [-x for x in reversed(word)]

def load_expected(name: str):
    return json.loads((EXPECTED / name).read_text(encoding='utf-8'))
