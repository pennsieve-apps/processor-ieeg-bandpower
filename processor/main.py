"""
processor-ieeg-bandpower: Welch spectra and band power per good channel of every recording.

  INPUT_DIR   the ieeg-prep BIDS-iEEG derivative (recordings.csv + EDF + sidecars), or any
              BIDS-iEEG dataset (*_ieeg.edf/.vhdr/.fif/.bdf/.set with _channels.tsv/_events.tsv)
  OUTPUT_DIR[/LABEL]/bandpower.csv    recording_id, channel, band, abs_power, rel_power, norm_median
  OUTPUT_DIR[/LABEL]/coherence.csv    recording_id, channel_a, channel_b, band, coherence (COHERENCE=true)
  OUTPUT_DIR[/LABEL]/provenance.json

One recording at a time, streamed in chunks; a recording that fails is logged and skipped.
"""
from __future__ import annotations

import csv
import json
import logging
import os
import platform
import sys
import time
import traceback
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

import numpy as np

from processor import inputs as I
from processor import spectral as S
from processor import window as W

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-ieeg-bandpower")
log.setLevel(logging.INFO)
logging.getLogger("mne").setLevel(logging.ERROR)

APP_NAME = "processor-ieeg-bandpower"
REPO_ROOT = Path(__file__).resolve().parent.parent
TRUE = ("1", "true", "yes", "on")

DEFAULTS = {
    "WINDOW": "baseline", "MARGIN_S": "30", "RESAMPLE_HZ": "200", "WELCH_WINDOW_S": "2", "WELCH_OVERLAP": "0.5",
    "BANDPASS": "0.5-80", "BANDS": "conventions", "COHERENCE": "false", "CHUNK_S": "60", "LABEL": "",
}


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def tool_versions() -> dict:
    out = {"python": platform.python_version()}
    for pkg in ("numpy", "scipy", "mne"):
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    out["spectral"] = S.SPECTRAL_VERSION
    return out


def read_params(env=os.environ) -> tuple[dict, S.SpectralParams, dict]:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    if raw["WINDOW"] not in W.MODES:
        raise ValueError(f"WINDOW must be one of {', '.join(W.MODES)}, got {raw['WINDOW']!r}")
    params = S.SpectralParams(resample_hz=float(raw["RESAMPLE_HZ"]), window_s=float(raw["WELCH_WINDOW_S"]),
                              overlap=float(raw["WELCH_OVERLAP"]), bandpass=S.parse_bandpass(raw["BANDPASS"]),
                              bands=S.parse_bands(raw["BANDS"]), chunk_s=float(raw["CHUNK_S"]))
    other = {"window": raw["WINDOW"], "margin_s": float(raw["MARGIN_S"]),
             "coherence": str(raw["COHERENCE"]).strip().lower() in TRUE,
             "label": raw["LABEL"].strip()}
    if other["label"] and not all(c.isalnum() or c in "-_" for c in other["label"]):
        raise ValueError("LABEL: letters, digits, - and _ only")
    return raw, params, other


def match_channels(raw_names: list[str], rec: I.Recording) -> tuple[list[int], list[str], list[str]]:
    """(raw indices, output names, notes) of the good channels."""
    notes: list[str] = []
    norm = {n.strip().lower(): i for i, n in enumerate(raw_names)}
    if rec.channels is None:
        skip = ("ecg", "ekg", "emg", "eog", "resp", "pleth", "spo2", "status", "trigger", "annotation", "dc")
        idx = [i for i, n in enumerate(raw_names) if not n.strip().lower().startswith(skip)]
        notes.append("no _channels.tsv: every signal channel used")
        return idx, [raw_names[i] for i in idx], notes
    idx, names, missing = [], [], []
    for ch in rec.channels:
        if ch.status != "good" or ch.type not in I.EEG_TYPES:
            continue
        i = raw_names.index(ch.name) if ch.name in raw_names else norm.get(ch.name.strip().lower())
        if i is None:
            missing.append(ch.name)
            continue
        idx.append(i)
        names.append(ch.name)
    if missing:
        notes.append(f"{len(missing)} good channel(s) of _channels.tsv not in the signal file: {', '.join(missing[:5])}")
    return idx, names, notes


def already_filtered(rec: I.Recording, names: list[str], band: tuple[float, float] | None) -> bool:
    """True when _channels.tsv says every used channel is already band-limited within BANDPASS
    (low_cutoff >= low and high_cutoff <= high), e.g. the HUP atlas clips (0.5-80 Hz)."""
    if band is None or not rec.channels:
        return False
    by = {c.name: c for c in rec.channels}
    for n in names:
        c = by.get(n)
        if c is None or c.low_cutoff is None or c.high_cutoff is None:
            return False
        if c.low_cutoff < band[0] - 1e-6 or c.high_cutoff > band[1] + 1e-6 or c.high_cutoff <= 0:
            return False
    return True


def process(rec: I.Recording, params: S.SpectralParams, other: dict, bp_writer, coh_writer) -> dict:
    import mne

    t0 = time.monotonic()
    raw = mne.io.read_raw(str(rec.path), preload=False, verbose="error")
    sfreq = float(raw.info["sfreq"])
    duration = raw.n_times / sfreq
    if sfreq < 2 * max(hi for _, _, hi in params.bands + (S.BROADBAND,)):
        raise ValueError(f"sampling rate {sfreq:g} Hz is too low for bands up to 80 Hz")
    idx, names, notes = match_channels(raw.ch_names, rec)
    if not idx:
        raise ValueError("no good channel")
    intervals_s, wnote = W.select(duration, rec.seizures, rec.baseline, rec.artifacts, other["window"], other["margin_s"])
    intervals = [(int(np.ceil(a * sfreq)), min(raw.n_times, int(np.floor(b * sfreq)))) for a, b in intervals_s]

    def read(a: int, b: int) -> np.ndarray:
        return raw.get_data(picks=idx, start=a, stop=b, units="uV", verbose="error")

    prefiltered = already_filtered(rec, names, params.bandpass)
    if prefiltered:
        notes.append(f"_channels.tsv cutoffs within {params.bandpass[0]:g}-{params.bandpass[1]:g} Hz: BANDPASS not applied again")
    acc, info = S.stream(read, intervals, sfreq, len(idx), params, prefiltered=prefiltered, cross=other["coherence"])
    if acc.n_windows == 0:
        raise ValueError(f"no complete {params.window_s:g} s window in the selected interval ({wnote})")
    freqs, psd = acc.psd()
    bt = S.band_table(freqs, psd, params.bands)
    for c, name in enumerate(names):
        for b, band in enumerate(bt.bands):
            bp_writer.writerow([rec.recording_id, name, band, f"{bt.abs_power[c, b]:.6g}", f"{bt.rel_power[c, b]:.6g}",
                                f"{bt.norm_median[c, b]:.6g}"])
    if other["coherence"]:
        cnames, cb = S.coherence_bands(freqs, acc.coherence(), params.bands)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                for b, band in enumerate(cnames):
                    coh_writer.writerow([rec.recording_id, names[i], names[j], band, f"{cb[i, j, b]:.6g}"])
    secs = time.monotonic() - t0
    n_good = len(idx)
    log.info("%s: %d good of %d channels, %g Hz, %s, %.1f s used (%d windows)%s, %.1fs", rec.recording_id, n_good,
             len(raw.ch_names), sfreq, wnote, info["seconds_used"], info["welch_windows"],
             ", already band-limited" if prefiltered else "", secs)
    return {"status": "ok", "n_channels": len(raw.ch_names), "n_good": n_good, "sfreq": sfreq,
            "duration_s": round(duration, 3), "window": wnote,
            "intervals_s": [[round(a, 3), round(b, 3)] for a, b in intervals_s], "bandpass_applied": not prefiltered,
            **info,
            "notes": notes + rec.notes, "seconds": round(secs, 2)}


def main() -> int:
    t_start = time.monotonic()
    started = datetime.now(timezone.utc)
    input_dir = Path(os.environ.get("INPUT_DIR", "/data/input"))
    output_dir = Path(os.environ.get("OUTPUT_DIR", "/data/output"))
    try:
        raw_params, params, other = read_params()
    except ValueError as e:
        log.error("parameters: %s", e)
        return 1
    out = output_dir / other["label"] if other["label"] else output_dir
    out.mkdir(parents=True, exist_ok=True)
    log.info("%s %s: INPUT_DIR=%s OUTPUT_DIR=%s window=%s bands=%s", APP_NAME, app_version(), input_dir, out,
             other["window"], ",".join(f"{n} {lo:g}-{hi:g}" for n, lo, hi in params.bands))
    recs = I.find_recordings(input_dir) if input_dir.is_dir() else []
    if not recs:
        log.error("no recording found under %s (recordings.csv or *_ieeg.edf/.vhdr/.fif/.bdf/.set)", input_dir)
        return 1

    def rel(p: Path) -> str:
        try:
            return str(p.relative_to(input_dir))
        except ValueError:
            return str(p)

    per_rec: dict[str, dict] = {}
    n_ok = 0
    bp_path, coh_path = out / "bandpower.csv", out / "coherence.csv"
    with open(bp_path, "w", newline="") as bp_fh, \
            (open(coh_path, "w", newline="") if other["coherence"] else open(os.devnull, "w")) as coh_fh:
        bp = csv.writer(bp_fh)
        bp.writerow(["recording_id", "channel", "band", "abs_power", "rel_power", "norm_median"])
        coh = csv.writer(coh_fh)
        if other["coherence"]:
            coh.writerow(["recording_id", "channel_a", "channel_b", "band", "coherence"])
        for rec in recs:
            try:
                per_rec[rec.recording_id] = {"path": rel(rec.path), **process(rec, params, other, bp, coh)}
                n_ok += 1
            except Exception as e:  # noqa: BLE001 - a recording that fails is logged and skipped
                log.error("%s: failed: %s", rec.recording_id, e)
                log.debug(traceback.format_exc())
                per_rec[rec.recording_id] = {"path": rel(rec.path), "status": "failed", "error": str(e)}
            bp_fh.flush()

    provenance = {
        "app": APP_NAME, "appVersion": app_version(), "tools": tool_versions(),
        "startedAt": started.isoformat(timespec="seconds"),
        "finishedAt": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t_start, 2),
        "parameters": raw_params, "spectral": params.as_dict(),
        "inputsFound": sorted(rel(r.path) for r in recs),
        "outputs": ["bandpower.csv"] + (["coherence.csv"] if other["coherence"] else []),
        "recordings": per_rec,
    }
    with open(out / "provenance.json", "w") as fh:
        json.dump(provenance, fh, indent=2)
    log.info("%d of %d recording(s) done, %.1fs; wrote %s", n_ok, len(recs), time.monotonic() - t_start, out)
    if n_ok == 0:
        log.error("every recording failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
