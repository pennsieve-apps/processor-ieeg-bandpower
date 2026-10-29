"""
Which part of a recording the spectra are computed over (WINDOW):

- baseline (default): for a recording with a known seizure, the baseline interval
  (recordings.csv baseline_start_s/_end_s, or a `baseline` event); when there is none, the
  stretch from the start to MARGIN_S before the first seizure onset. A recording without a
  known seizure is used whole.
- all: the whole recording.
- interictal-only: the whole recording minus every seizure and MARGIN_S around it (from
  MARGIN_S before onset to MARGIN_S after its end; a seizure of unknown duration is cut to
  the end of the recording).

Artifact events are always cut out.
"""
from __future__ import annotations

MODES = ("baseline", "all", "interictal-only")


def subtract(intervals: list[tuple[float, float]], cut: tuple[float, float]) -> list[tuple[float, float]]:
    out = []
    c0, c1 = cut
    for a, b in intervals:
        if c1 <= a or c0 >= b:
            out.append((a, b))
            continue
        if c0 > a:
            out.append((a, c0))
        if c1 < b:
            out.append((c1, b))
    return out


def select(duration: float, seizures: list[tuple[float, float | None]], baseline: tuple[float, float] | None,
           artifacts: list[tuple[float, float]], mode: str, margin: float) -> tuple[list[tuple[float, float]], str]:
    """(intervals in seconds, a note saying what was chosen)."""
    if mode not in MODES:
        raise ValueError(f"WINDOW must be one of {', '.join(MODES)}")
    whole = [(0.0, duration)]
    if mode == "all":
        iv, note = whole, "all"
    elif mode == "interictal-only":
        iv = whole
        for on, dur in seizures:
            end = on + dur if dur else duration
            iv = subtract(iv, (on - margin, end + margin))
        note = f"interictal-only ({len(seizures)} seizure(s) cut with {margin:g} s margin)" if seizures else \
            "interictal-only (no seizure known: whole recording)"
    else:  # baseline
        if not seizures:
            iv, note = whole, "baseline: no seizure known, whole recording used"
        elif baseline:
            iv, note = [(max(0.0, baseline[0]), min(duration, baseline[1]))], "baseline interval"
        else:
            first = min(on for on, _ in seizures)
            iv, note = [(0.0, max(0.0, first - margin))], f"baseline: none given, 0 s to {margin:g} s before onset"
    for art in artifacts:
        iv = subtract(iv, art)
    iv = [(a, b) for a, b in iv if b > a]
    return iv, note
