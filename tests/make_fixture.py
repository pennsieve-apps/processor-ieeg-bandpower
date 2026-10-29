"""
A tiny synthetic ieeg-prep derivative (conventions.md layout) for the tests and `make run`:

  <out>/linked/ieeg/            the derivative, as the `ieeg` input port would deliver it
    dataset_description.json, recordings.csv
    sub-P01/ses-implant1/ieeg/  an ictal clip (512 Hz, 100 s, 6 bipolar channels, one bad):
                                seizure at 60 s for 20 s (20 Hz burst on LA02-LA03, LB02-LB03),
                                baseline 5-50 s, an artifact at 20-22 s with a 40 Hz burst on
                                every channel, a 10 Hz rhythm on LA01-LA02, 2 Hz on LB01-LB02;
                                and an interictal clip (40 s, no seizure)
    sub-HUPatlas001/ses-atlas/ieeg/  an atlas-like clip (200 Hz, 60 s, 4 channels) whose
                                electrodes.tsv carries atlas_patient = 1
  <out>/plain/                  the same ictal clip as a plain BIDS-iEEG dataset (no recordings.csv,
                                `ictal` event, no baseline) for the no-derivative path

Usage: python -m tests.make_fixture <out-dir>
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

RNG_SEED = 42


def write_edf(path: Path, data_uv: np.ndarray, sfreq: int, labels: list[str]) -> None:
    """Minimal EDF (16-bit, 1 s records, per-channel physical range from the data, units uV)."""
    n_ch, n = data_uv.shape
    n_rec = n // sfreq
    data_uv = data_uv[:, : n_rec * sfreq]
    pmax = np.maximum(np.abs(data_uv).max(axis=1) * 1.1, 1.0)
    pmin = -pmax
    dmin, dmax = -32768, 32767
    digital = np.round((data_uv - pmin[:, None]) / (pmax - pmin)[:, None] * (dmax - dmin) + dmin).astype("<i2")

    def f(s, w):
        return str(s)[:w].ljust(w).encode("ascii")

    hdr = b"".join([f("0", 8), f("X X X X", 80), f("Startdate 01-JAN-2020 X X X", 80), f("01.01.20", 8),
                    f("00.00.00", 8), f(256 * (n_ch + 1), 8), f("", 44), f(n_rec, 8), f(1, 8), f(n_ch, 4)])
    fields = [[f(lab, 16) for lab in labels], [f("", 80)] * n_ch, [f("uV", 8)] * n_ch,
              [f(f"{v:.3f}", 8) for v in pmin], [f(f"{v:.3f}", 8) for v in pmax], [f(dmin, 8)] * n_ch,
              [f(dmax, 8)] * n_ch, [f("", 80)] * n_ch, [f(sfreq, 8)] * n_ch, [f("", 32)] * n_ch]
    hdr += b"".join(b"".join(col) for col in fields)
    body = digital.reshape(n_ch, n_rec, sfreq).transpose(1, 0, 2).tobytes()
    path.write_bytes(hdr + body)


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows([["n/a" if v is None or v == "" else v for v in r] for r in rows])


def background(rng, n_ch: int, n: int, sfreq: int) -> np.ndarray:
    """1/f-like noise, ~20 uV rms."""
    white = rng.standard_normal((n_ch, n))
    spec = np.fft.rfft(white, axis=1)
    f = np.fft.rfftfreq(n, 1 / sfreq)
    spec /= np.maximum(f, 0.5) ** 0.5
    x = np.fft.irfft(spec, n=n, axis=1)
    return 20 * x / x.std(axis=1, keepdims=True)


CH_COLS = ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "status",
           "status_description", "contacts", "soz", "resected", "onset", "spread", "x", "y", "z", "space", "roi"]
ICTAL_CHANNELS = [
    # name, status, soz, resected, onset, xyz (MNI, midpoint of the pair)
    ("LA01-LA02", "good", "false", "true", "false", (-24.0, -8.0, -20.0)),
    ("LA02-LA03", "good", "true", "true", "true", (-28.0, -8.0, -20.0)),
    ("LB01-LB02", "good", "false", "false", "false", (-30.0, -25.0, -10.0)),
    ("LB02-LB03", "good", "true", "false", "true", (-34.0, -25.0, -10.0)),
    ("LC01-LC02", "bad", "false", "false", "false", (-40.0, 10.0, 30.0)),
    ("RA01-RA02", "good", "false", "false", "false", (24.0, -8.0, -20.0)),
]


def ictal_signal(rng, sfreq=512, dur=100, onset=60.0, sz_dur=20.0) -> np.ndarray:
    n = sfreq * dur
    t = np.arange(n) / sfreq
    x = background(rng, len(ICTAL_CHANNELS), n, sfreq)
    x[0] += 60 * np.sin(2 * np.pi * 10 * t)
    x[2] += 80 * np.sin(2 * np.pi * 2 * t)
    sz = (t >= onset) & (t < onset + sz_dur)
    for c in (1, 3):
        x[c, sz] += 150 * np.sin(2 * np.pi * 20 * t[sz])
    art = (t >= 20) & (t < 22)
    x[:, art] += 300 * np.sin(2 * np.pi * 40 * t[art])
    for c in (1, 3):  # interictal spikes
        for s in rng.uniform(5, 55, 8):
            i = int(s * sfreq)
            x[c, i: i + 20] += 200 * np.hanning(20)
    return x


def make_derivative(root: Path, rng) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "dataset_description.json").write_text(json.dumps({
        "Name": "synthetic ieeg-prep fixture", "BIDSVersion": "1.9.0", "DatasetType": "derivative",
        "GeneratedBy": [{"Name": "ieeg-prep", "Version": "fixture"}]}, indent=2))
    rows = []

    # sub-P01 ictal clip
    d = root / "sub-P01" / "ses-implant1" / "ieeg"
    d.mkdir(parents=True, exist_ok=True)
    stem = "sub-P01_ses-implant1_task-ictal01_run-01"
    x = ictal_signal(rng)
    names = [c[0] for c in ICTAL_CHANNELS]
    write_edf(d / f"{stem}_desc-preproc_ieeg.edf", x, 512, names)
    write_tsv(d / f"{stem}_desc-preproc_channels.tsv", CH_COLS, [
        [n, "SEEG", "µV", 0.5, 200, 512, st, "flat" if st == "bad" else "", n.replace("-", "|"), soz, res, on, "",
         *xyz, "MNI152NLin2009aSym", ""] for n, st, soz, res, on, xyz in ICTAL_CHANNELS])
    write_tsv(d / f"{stem}_desc-preproc_events.tsv", ["onset", "duration", "trial_type", "source"], [
        [60.0, 20.0, "seizure", "fixture"], [5.0, 45.0, "baseline", "fixture"], [20.0, 2.0, "artifact", "fixture"]])
    (d / f"{stem}_desc-preproc_ieeg.json").write_text(json.dumps({
        "SamplingFrequency": 512, "iEEGReference": "bipolar", "PowerLineFrequency": 60}, indent=2))
    rows.append([stem, "P01", "implant1", "ictal01", "01", f"sub-P01/ses-implant1/ieeg/{stem}_desc-preproc_ieeg.edf",
                 512, 6, 5, 100, "bipolar", "ictal", 60, 5, 50, "bids-edf"])

    # sub-P01 interictal clip
    stem2 = "sub-P01_ses-implant1_task-interictal01_run-01"
    x2 = background(rng, len(ICTAL_CHANNELS), 512 * 40, 512)
    write_edf(d / f"{stem2}_desc-preproc_ieeg.edf", x2, 512, names)
    write_tsv(d / f"{stem2}_desc-preproc_channels.tsv", CH_COLS, [
        [n, "SEEG", "µV", 0.5, 200, 512, st, "", n.replace("-", "|"), soz, res, "", "", *xyz, "MNI152NLin2009aSym", ""]
        for n, st, soz, res, on, xyz in ICTAL_CHANNELS])
    write_tsv(d / f"{stem2}_desc-preproc_events.tsv", ["onset", "duration", "trial_type", "source"], [])
    rows.append([stem2, "P01", "implant1", "interictal01", "01",
                 f"sub-P01/ses-implant1/ieeg/{stem2}_desc-preproc_ieeg.edf", 512, 6, 5, 40, "bipolar", "interictal",
                 "", "", "", "bids-edf"])
    write_tsv(d / "sub-P01_ses-implant1_desc-preproc_electrodes.tsv", ["name", "x", "y", "z", "size", "atlas_patient"],
              [[c, -22.0 - 2 * i, -8.0, -20.0, "n/a", ""] for i, c in enumerate(["LA01", "LA02", "LA03"])])

    # atlas-like subject
    d3 = root / "sub-HUPatlas001" / "ses-atlas" / "ieeg"
    d3.mkdir(parents=True, exist_ok=True)
    stem3 = "sub-HUPatlas001_ses-atlas_task-interictal_run-01"
    atlas_ch = [("LA1-LA2", (-25.0, -8.0, -19.0), "true"), ("LA2-LA3", (-29.0, -8.0, -19.0), "false"),
                ("RB1-RB2", (30.0, -20.0, -12.0), "false"), ("RB2-RB3", (34.0, -20.0, -12.0), "false")]
    x3 = background(rng, 4, 200 * 60, 200)
    write_edf(d3 / f"{stem3}_desc-preproc_ieeg.edf", x3, 200, [c[0] for c in atlas_ch])
    write_tsv(d3 / f"{stem3}_desc-preproc_channels.tsv", CH_COLS, [
        [n, "SEEG", "µV", 0.5, 80, 200, "good", "", n.replace("-", "|"), soz, soz, "", "", *xyz, "MNI152", ""]
        for n, xyz, soz in atlas_ch])
    write_tsv(d3 / f"{stem3}_desc-preproc_events.tsv", ["onset", "duration", "trial_type", "source"], [])
    write_tsv(d3 / "sub-HUPatlas001_ses-atlas_desc-preproc_electrodes.tsv",
              ["name", "x", "y", "z", "size", "atlas_patient"],
              [["LA1", -23, -8, -19, "n/a", 1], ["LA2", -27, -8, -19, "n/a", 1], ["LA3", -31, -8, -19, "n/a", 1]])
    rows.append([stem3, "HUPatlas001", "atlas", "interictal", "01",
                 f"sub-HUPatlas001/ses-atlas/ieeg/{stem3}_desc-preproc_ieeg.edf", 200, 4, 4, 60, "bipolar",
                 "interictal", "", "", "", "mat"])

    with open(root / "recordings.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["recording_id", "subject", "session", "task", "run", "path", "sfreq", "n_channels", "n_good",
                    "duration_s", "montage", "kind", "onset_s", "baseline_start_s", "baseline_end_s", "source_format"])
        w.writerows(rows)


def make_plain(root: Path, rng) -> None:
    d = root / "sub-P02" / "ses-implant1" / "ieeg"
    d.mkdir(parents=True, exist_ok=True)
    stem = "sub-P02_ses-implant1_task-ictal01_run-01"
    x = ictal_signal(rng)
    write_edf(d / f"{stem}_ieeg.edf", x, 512, [c[0] for c in ICTAL_CHANNELS])
    write_tsv(d / f"{stem}_channels.tsv", ["name", "type", "units", "sampling_frequency", "status"],
              [[c[0], "SEEG", "uV", 512, c[1]] for c in ICTAL_CHANNELS])
    write_tsv(d / f"{stem}_events.tsv", ["onset", "duration", "trial_type"], [[60.0, 20.0, "ictal"]])


def main(out: str) -> None:
    rng = np.random.default_rng(RNG_SEED)
    out_p = Path(out)
    make_derivative(out_p / "linked" / "ieeg", rng)
    make_plain(out_p / "plain", rng)
    print(f"fixture written to {out_p}/linked/ieeg and {out_p}/plain")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/fixture")
