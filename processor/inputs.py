"""
Finding the recordings under INPUT_DIR: the ieeg-prep derivative (recordings.csv) or a plain
BIDS-iEEG dataset (*_ieeg.edf/.vhdr/.fif/.bdf/.set), with their _channels.tsv and _events.tsv.
"""
from __future__ import annotations

import csv
import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("processor-ieeg-bandpower")

SIGNAL_EXTS = (".edf", ".bdf", ".vhdr", ".fif", ".set")
EEG_TYPES = {"SEEG", "ECOG", "EEG", "DBS", ""}
NA = {"", "n/a", "na", "nan", "none"}


def walk_files(root: Path) -> list[Path]:
    """Every file under root, following symlinked folders and files (INPUT_DIR entries are symlinks)."""
    found: list[Path] = []
    seen: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        for name in sorted(filenames):
            if not name.startswith("."):
                p = Path(dirpath) / name
                if p.is_file():
                    found.append(p)
    return found


def read_table(path: Path, delimiter: str | None = None) -> list[dict[str, str]]:
    """Rows of a TSV/CSV as dicts (BOM stripped, keys and values stripped, n/a -> '')."""
    delim = delimiter or ("\t" if path.suffix.lower() == ".tsv" else ",")
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = []
        for raw in csv.DictReader(fh, delimiter=delim):
            row = {}
            for k, v in raw.items():
                if k is None:
                    continue
                v = (v or "").strip()
                row[k.strip()] = "" if v.lower() in NA else v
            rows.append(row)
        return rows


def fnum(v: str | None) -> float | None:
    try:
        x = float(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return x if x == x else None


def recording_id_from(path: Path) -> str:
    """BIDS stem of a signal file: name without _ieeg.<ext> and without the desc-preproc entity."""
    name = path.name
    for ext in SIGNAL_EXTS:
        if name.lower().endswith("_ieeg" + ext):
            name = name[: -len("_ieeg" + ext)]
            break
    else:
        name = path.stem
    return re.sub(r"_desc-preproc$", "", name)


def sidecar(signal_path: Path, suffix: str) -> Path | None:
    """<prefix>_<suffix> next to <prefix>_ieeg.<ext> (suffix e.g. 'channels.tsv')."""
    name = signal_path.name
    m = re.match(r"^(.*)_ieeg\.[A-Za-z0-9]+$", name)
    if not m:
        return None
    p = signal_path.with_name(f"{m.group(1)}_{suffix}")
    return p if p.is_file() else None


@dataclass
class Channel:
    name: str
    type: str = ""
    status: str = "good"
    low_cutoff: float | None = None
    high_cutoff: float | None = None


@dataclass
class Recording:
    recording_id: str
    path: Path
    subject: str = ""
    kind: str = "unknown"
    onset_s: float | None = None
    baseline: tuple[float, float] | None = None
    source: str = ""            # recordings.csv or bids
    channels: list[Channel] | None = None
    seizures: list[tuple[float, float | None]] = field(default_factory=list)
    artifacts: list[tuple[float, float]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def subject_of(recording_id: str) -> str:
    m = re.match(r"^sub-([A-Za-z0-9]+)", recording_id)
    return m.group(1) if m else ""


def _is_seizure(trial_type: str) -> bool:
    """Seizure events; also stimulation (`stim`: a stimulation-induced seizure in Discover 530),
    which must not count as baseline either."""
    t = trial_type.strip().lower()
    return t in ("seizure", "ictal", "sz", "stim", "stimulation") or t.startswith("seizure")


def _is_artifact(trial_type: str) -> bool:
    t = trial_type.strip().lower()
    return "artifact" in t or "artefact" in t


def load_sidecars(rec: Recording) -> None:
    """Channels (_channels.tsv) and events (_events.tsv) next to the signal file, if present."""
    ch = sidecar(rec.path, "channels.tsv")
    if ch:
        rec.channels = [Channel(r.get("name", ""), r.get("type", "").upper(), (r.get("status") or "good").lower(),
                                fnum(r.get("low_cutoff")), fnum(r.get("high_cutoff")))
                        for r in read_table(ch) if r.get("name")]
    ev = sidecar(rec.path, "events.tsv")
    if ev:
        for r in read_table(ev):
            on, dur, tt = fnum(r.get("onset")), fnum(r.get("duration")), r.get("trial_type", "")
            if on is None:
                continue
            if _is_seizure(tt):
                rec.seizures.append((on, dur if dur and dur > 0 else None))
            elif tt.strip().lower() == "baseline" and rec.baseline is None and dur:
                rec.baseline = (on, on + dur)
            elif _is_artifact(tt) and dur:
                rec.artifacts.append((on, on + dur))
    if rec.onset_s is not None and not any(abs(s - rec.onset_s) < 1e-6 for s, _ in rec.seizures):
        rec.seizures.insert(0, (rec.onset_s, None))
    if rec.onset_s is None and rec.seizures:
        rec.onset_s = min(s for s, _ in rec.seizures)


def find_recordings(input_dir: Path) -> list[Recording]:
    """Recordings listed in any recordings.csv, plus signal files not listed (plain BIDS)."""
    files = walk_files(input_dir)
    recs: dict[str, Recording] = {}
    listed: set[str] = set()
    for idx in [p for p in files if p.name == "recordings.csv"]:
        for r in read_table(idx, ","):
            rid, rel = r.get("recording_id", ""), r.get("path", "")
            if not rid or not rel:
                continue
            path = (idx.parent / rel)
            if not path.is_file():
                log.warning("%s: %s listed in %s but not found", rid, rel, idx)
                continue
            listed.add(os.path.realpath(path))
            if rid in recs:
                log.warning("%s listed twice; using the first (%s)", rid, recs[rid].path)
                continue
            bs, be = fnum(r.get("baseline_start_s")), fnum(r.get("baseline_end_s"))
            recs[rid] = Recording(rid, path, subject=r.get("subject") or subject_of(rid),
                                  kind=(r.get("kind") or "unknown").lower(), onset_s=fnum(r.get("onset_s")),
                                  baseline=(bs, be) if bs is not None and be is not None and be > bs else None,
                                  source="recordings.csv")
    for p in files:
        if not any(p.name.lower().endswith("_ieeg" + e) for e in SIGNAL_EXTS):
            continue
        if os.path.realpath(p) in listed:
            continue
        rid = recording_id_from(p)
        if rid in recs:
            log.warning("%s: %s duplicates a recording already found (%s); skipped", rid, p, recs[rid].path)
            continue
        recs[rid] = Recording(rid, p, subject=subject_of(rid), source="bids")
    for rec in recs.values():
        load_sidecars(rec)
        if rec.kind == "unknown" and rec.seizures:
            rec.kind = "ictal"
    return [recs[k] for k in sorted(recs)]
