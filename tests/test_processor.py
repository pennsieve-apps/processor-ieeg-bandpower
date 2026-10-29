import csv
import json
import os
from pathlib import Path

import numpy as np
import pytest
from scipy import signal

from processor import main as M
from processor import spectral as S
from processor import window as W
from tests import make_fixture as F

ICTAL = "sub-P01_ses-implant1_task-ictal01_run-01"
INTERICTAL = "sub-P01_ses-implant1_task-interictal01_run-01"
ATLAS = "sub-HUPatlas001_ses-atlas_task-interictal_run-01"


@pytest.fixture(scope="session")
def fixture_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("fixture")
    F.main(str(d))
    return d


def run(monkeypatch, input_dir, output_dir, **env):
    for k in M.DEFAULTS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for k, v in env.items():
        monkeypatch.setenv(k, str(v))
    return M.main()


def table(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def value(rows, rid, ch, band, col):
    for r in rows:
        if (r["recording_id"], r["channel"], r["band"]) == (rid, ch, band):
            return float(r[col])
    raise KeyError((rid, ch, band))


# --- spectral --------------------------------------------------------------------------

def test_accumulator_matches_scipy_welch():
    rng = np.random.default_rng(0)
    x = rng.standard_normal((3, 200 * 30))
    p = S.SpectralParams()
    acc = S.Accumulator(3, 200.0, p)
    acc.add(x[:, :1234])          # in uneven pieces: windows run across them
    acc.add(x[:, 1234:4000])
    acc.add(x[:, 4000:])
    f, psd = acc.psd()
    f2, p2 = signal.welch(x, fs=200, window=signal.windows.hamming(400, sym=True), nperseg=400, noverlap=200,
                          detrend="constant")
    assert np.allclose(f, f2) and np.allclose(psd, p2, rtol=1e-10)
    assert acc.n_windows == 29


def test_stream_chunks_equal_one_piece():
    rng = np.random.default_rng(1)
    x = rng.standard_normal((2, 1024 * 200))
    p_small = S.SpectralParams(chunk_s=13)
    p_big = S.SpectralParams(chunk_s=1000)
    a1, _ = S.stream(lambda a, b: x[:, a:b], [(0, x.shape[1])], 1024, 2, p_small)
    a2, _ = S.stream(lambda a, b: x[:, a:b], [(0, x.shape[1])], 1024, 2, p_big)
    assert a1.n_windows == a2.n_windows
    assert np.allclose(a1.psd()[1], a2.psd()[1], rtol=1e-3)


def test_coherence_matches_scipy_and_paper_windows():
    """1 s Hamming windows, no overlap (mscohere(x, y, 200, 0, 1:80, 200) at 200 Hz) = scipy.signal.coherence."""
    rng = np.random.default_rng(3)
    common = rng.standard_normal(200 * 60)
    x = np.vstack([common + rng.standard_normal(len(common)) * k for k in (0.5, 1.0, 3.0)])
    p = S.SpectralParams()
    acc, info = S.stream(lambda a, b: x[:, a:b], [(0, x.shape[1])], 200, 3, p, prefiltered=True, cross=True)
    assert info["coherence_windows"] == 60
    c = acc.coh.coherence()
    f2, c2 = signal.coherence(x[0], x[2], fs=200, window=signal.windows.hamming(200, sym=True), nperseg=200,
                              noverlap=0, detrend="constant")
    assert np.allclose(c[0, 2], c2[: c.shape[-1]], rtol=1e-8)
    names, cb = S.coherence_bands(acc.coh.freqs[: acc.coh.nf], c, p.bands)
    alpha = (f2 >= 8) & (f2 < 13)
    assert cb[0, 2, 2] == pytest.approx(c2[alpha].mean())
    strength = S.coherence_strength(cb)
    assert strength[0, 0] == pytest.approx((cb[0, 1, 0] + cb[0, 2, 0]) / 2)
    assert strength[0, 0] > strength[2, 0]          # the noisiest channel is the least coherent


def test_feature_family():
    assert S.feature_family("normpow_delta") == "spectral" and S.feature_family("coh_gamma") == "coherence"
    assert S.feature_family("spike_rate") == "other" and S.feature_family("x") == "other"


def test_band_table_definitions():
    f = np.arange(0, 100.5, 0.5)
    psd = np.ones((1, len(f)))
    t = S.band_table(f, psd, S.BAND_PRESETS["conventions"])
    assert t.bands == ["delta", "theta", "alpha", "beta", "gamma", "broadband"]
    # flat spectrum: absolute power = width, relative power = width / 79
    assert np.allclose(t.abs_power[0], [3, 4, 5, 17, 50, 79])
    assert np.allclose(t.rel_power[0], np.array([3, 4, 5, 17, 50, 79]) / 79)
    # normalised to sum 1 over 0.5-80 Hz (160 bins), so every median is 1/160
    assert np.allclose(t.norm_median[0], 1 / 160)
    # a peak at 10 Hz moves alpha's share up and leaves the other medians
    psd[0, f == 10] = 100
    t2 = S.band_table(f, psd, S.BAND_PRESETS["conventions"])
    assert t2.rel_power[0, 2] > 0.4 > t.rel_power[0, 2]
    assert t2.norm_median[0, 0] < 1 / 160


def test_sine_lands_in_its_band():
    fs = 512
    t = np.arange(fs * 60) / fs
    x = np.vstack([np.sin(2 * np.pi * hz * t) * 50 + np.random.default_rng(2).standard_normal(len(t))
                   for hz in (2.5, 6, 10, 20, 45)])
    p = S.SpectralParams()
    acc, info = S.stream(lambda a, b: x[:, a:b], [(0, x.shape[1])], fs, 5, p)
    assert info["fs"] == 200
    tab = S.band_table(*acc.psd(), p.bands)
    assert list(np.argmax(tab.rel_power[:, :5], axis=1)) == [0, 1, 2, 3, 4]


def test_parse_bands_and_bandpass():
    assert S.parse_bands("bernabei2022")[2] == ("alpha", 8.0, 12.0)
    assert S.parse_bands('{"slow": [1, 4], "fast": [30, 70]}') == (("slow", 1.0, 4.0), ("fast", 30.0, 70.0))
    with pytest.raises(ValueError):
        S.parse_bands('{"broadband": [1, 80]}')
    with pytest.raises(ValueError):
        S.parse_bands("theta-ish")
    assert S.parse_bandpass("none") is None and S.parse_bandpass("0.5-80") == (0.5, 80.0)
    with pytest.raises(ValueError):
        S.SpectralParams(resample_hz=100)


def test_resample_ratio():
    assert S.resample_ratio(1024, 200) == (25, 128)
    assert S.resample_ratio(200, 200) == (1, 1)
    assert S.resample_ratio(500, 200) == (2, 5)


# --- window ----------------------------------------------------------------------------

def test_window_modes():
    sz = [(60.0, 20.0)]
    iv, _ = W.select(100, sz, (5, 50), [(20, 22)], "baseline", 30)
    assert iv == [(5, 20), (22, 50)]
    iv, note = W.select(100, sz, None, [], "baseline", 30)
    assert iv == [(0.0, 30.0)] and "before onset" in note
    iv, _ = W.select(100, [], None, [], "baseline", 30)
    assert iv == [(0.0, 100)]
    iv, _ = W.select(200, sz, None, [], "interictal-only", 30)
    assert iv == [(0.0, 30.0), (110.0, 200)]
    iv, _ = W.select(200, [(60.0, None)], None, [], "interictal-only", 10)
    assert iv == [(0.0, 50.0)]
    iv, _ = W.select(100, sz, (5, 50), [], "all", 30)
    assert iv == [(0.0, 100)]


# --- end to end ------------------------------------------------------------------------

def test_run_on_derivative(fixture_dir, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(monkeypatch, fixture_dir / "linked", out) == 0          # COHERENCE defaults to true
    rows = table(out / "bandpower.csv")
    assert list(rows[0]) == ["recording_id", "channel", "band", "abs_power", "rel_power", "norm_median"]
    rids = {r["recording_id"] for r in rows}
    assert rids == {ICTAL, INTERICTAL, ATLAS}
    chans = {r["channel"] for r in rows if r["recording_id"] == ICTAL}
    assert "LC01-LC02" not in chans and len(chans) == 5        # bad channel left out
    assert {r["band"] for r in rows} == {"delta", "theta", "alpha", "beta", "gamma", "broadband"}
    for rid, ch in [(ICTAL, "LA01-LA02"), (ATLAS, "LA1-LA2")]:
        shares = sum(value(rows, rid, ch, b, "rel_power") for b in ("delta", "theta", "alpha", "beta", "gamma"))
        assert shares == pytest.approx(1.0, abs=1e-3)
        assert value(rows, rid, ch, "broadband", "rel_power") == pytest.approx(1.0)
    # 10 Hz rhythm and 2 Hz rhythm show up in their bands
    assert value(rows, ICTAL, "LA01-LA02", "alpha", "rel_power") > 0.5
    assert value(rows, ICTAL, "LB01-LB02", "delta", "rel_power") > 0.5
    prov = json.loads((out / "provenance.json").read_text())
    rec = prov["recordings"][ICTAL]
    # baseline 5-50 s minus the artifact 20-22 s
    assert rec["intervals_s"] == [[5.0, 20.0], [22.0, 50.0]] and rec["status"] == "ok"
    assert rec["bandpass_applied"] is True
    assert prov["recordings"][ATLAS]["bandpass_applied"] is False   # cutoffs 0.5-80 in _channels.tsv
    assert prov["spectral"]["bands"]["alpha"] == [8.0, 13.0]
    coh = table(out / "coherence.csv")
    assert len([r for r in coh if r["recording_id"] == ICTAL]) == 10 * 5   # 5 channels -> 10 pairs x 5 bands
    assert all(0 <= float(r["coherence"]) <= 1 for r in coh)
    feats = table(out / "features.csv")
    assert list(feats[0]) == ["recording_id", "channel", "feature", "value"]
    names = {r["feature"] for r in feats}
    bands = ("delta", "theta", "alpha", "beta", "gamma")
    assert names == {f"{p}_{b}" for p in ("relpow", "normpow", "logpow", "coh") for b in bands}
    fv = {(r["recording_id"], r["channel"], r["feature"]): float(r["value"]) for r in feats}
    assert fv[(ICTAL, "LA01-LA02", "relpow_alpha")] == pytest.approx(value(rows, ICTAL, "LA01-LA02", "alpha", "rel_power"), rel=1e-4)
    assert fv[(ICTAL, "LA01-LA02", "normpow_delta")] == pytest.approx(value(rows, ICTAL, "LA01-LA02", "delta", "norm_median"), rel=1e-4)
    assert fv[(ICTAL, "LA01-LA02", "logpow_beta")] == pytest.approx(np.log10(value(rows, ICTAL, "LA01-LA02", "beta", "abs_power")), abs=1e-4)
    # coh_<band> = mean of the channel's pair coherences in coherence.csv
    pairs = [float(r["coherence"]) for r in coh if r["recording_id"] == ICTAL and r["band"] == "theta"
             and "LA01-LA02" in (r["channel_a"], r["channel_b"])]
    assert len(pairs) == 4 and fv[(ICTAL, "LA01-LA02", "coh_theta")] == pytest.approx(np.mean(pairs), rel=1e-4)
    assert prov["recordings"][ICTAL]["coherence_windows"] > 0 and prov["spectral"]["coherence"]["window_s"] == 1.0


def test_window_changes_seizure_band(fixture_dir, tmp_path, monkeypatch):
    """The 20 Hz seizure burst on LA02-LA03 is outside the baseline and inside WINDOW=all."""
    assert run(monkeypatch, fixture_dir / "linked", tmp_path / "b") == 0
    assert run(monkeypatch, fixture_dir / "linked", tmp_path / "a", WINDOW="all") == 0
    base = value(table(tmp_path / "b" / "bandpower.csv"), ICTAL, "LA02-LA03", "beta", "rel_power")
    allw = value(table(tmp_path / "a" / "bandpower.csv"), ICTAL, "LA02-LA03", "beta", "rel_power")
    assert allw > 0.8 > 0.5 > base


def test_plain_bids_and_label(fixture_dir, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(monkeypatch, fixture_dir / "plain", out, LABEL="pre", COHERENCE="false") == 0
    prov = json.loads((out / "pre" / "provenance.json").read_text())
    rec = prov["recordings"]["sub-P02_ses-implant1_task-ictal01_run-01"]
    assert rec["intervals_s"] == [[0.0, 30.0]]          # `ictal` event at 60 s, no baseline, 30 s margin
    assert not (out / "pre" / "coherence.csv").exists()
    assert not any(r["feature"].startswith("coh_") for r in table(out / "pre" / "features.csv"))


def test_failed_recording_is_skipped(fixture_dir, tmp_path, monkeypatch):
    import shutil
    src = tmp_path / "in"
    shutil.copytree(fixture_dir / "plain", src)
    broken = src / "sub-P03" / "ieeg"
    broken.mkdir(parents=True)
    (broken / "sub-P03_task-x_ieeg.edf").write_bytes(b"not an edf")
    out = tmp_path / "out"
    assert run(monkeypatch, src, out) == 0
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["recordings"]["sub-P03_task-x"]["status"] == "failed"
    # only broken input: the run fails
    shutil.rmtree(src / "sub-P02")
    assert run(monkeypatch, src, tmp_path / "out2") == 1


def test_no_input_fails(tmp_path, monkeypatch):
    (tmp_path / "empty").mkdir()
    assert run(monkeypatch, tmp_path / "empty", tmp_path / "out") == 1


def test_bad_parameter_fails(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir / "plain", tmp_path / "o", WINDOW="ictal") == 1
    assert run(monkeypatch, fixture_dir / "plain", tmp_path / "o", RESAMPLE_HZ="100") == 1


def test_symlinked_input(fixture_dir, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    inp.mkdir()
    os.symlink(fixture_dir / "linked" / "ieeg", inp / "ieeg")
    assert run(monkeypatch, inp, tmp_path / "out") == 0
    assert len({r["recording_id"] for r in table(tmp_path / "out" / "bandpower.csv")}) == 3


def test_lambda_handler(fixture_dir, tmp_path, monkeypatch):
    from processor import handler
    monkeypatch.setenv("COHERENCE", "false")          # left over from a previous invocation
    res = handler.handler({"inputDir": str(fixture_dir / "plain"), "outputDir": str(tmp_path / "out"),
                           "WINDOW": "all"}, None)
    assert res["status"] == "ok"
    assert os.environ.get("COHERENCE") is None and (tmp_path / "out" / "coherence.csv").exists()
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(tmp_path / "nothing"), "outputDir": str(tmp_path / "o2")}, None)


def test_spectral_module_is_vendored_identically():
    other = Path(__file__).resolve().parents[2] / "processor-ieeg-normative" / "processor" / "spectral.py"
    if not other.exists():
        pytest.skip("processor-ieeg-normative not checked out next to this repo")
    mine = Path(__file__).resolve().parents[1] / "processor" / "spectral.py"
    assert mine.read_bytes() == other.read_bytes(), "spectral.py differs from processor-ieeg-normative's copy"
