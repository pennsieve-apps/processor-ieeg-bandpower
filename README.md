# processor-ieeg-bandpower

A step of the **Interictal Biomarker Mapping** workflow
(`ieeg-prep → {ieeg-spikes, ieeg-bandpower → ieeg-normative, ieeg-hfo} → ieeg-report`).
It measures, for every good channel of every intracranial EEG recording (SEEG or ECoG), how
the signal's power is spread over the classic frequency bands.

- **For the epileptologist:** epileptogenic tissue often shows a different interictal spectrum
  from healthy cortex at the same place: more slow (delta) activity, sometimes less beta and
  gamma. This step turns each channel into five numbers per measure, so `ieeg-normative` can ask
  "is this channel's spectrum unusual for this brain region?" against an atlas of normal
  channels. It uses the pre-ictal baseline of a seizure recording by default, so seizure
  activity does not dominate the spectrum.
- **For the researcher:** a faithful reimplementation of the univariate features of Bernabei
  et al. 2022 (Brain): 0.5–80 Hz band-pass, resampling to 200 Hz, Welch PSD (2 s Hamming
  windows, 1 s overlap), PSD normalised to sum to one over 0.5–80 Hz and the median taken per
  band; plus absolute and relative band power and, optionally, band coherence. The spectral code
  is shared byte-for-byte with `ieeg-normative`, which computes the atlas reference with it, so
  test and reference values are computed the same way.

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `ieeg` | input | subfolder | The ieeg-prep BIDS-iEEG derivative (`recordings.csv`, `*_desc-preproc_ieeg.edf` with `_channels.tsv` and `_events.tsv`), or any BIDS-iEEG dataset (`*_ieeg.edf/.bdf/.vhdr/.fif/.set`, sidecars optional). Found anywhere under `INPUT_DIR`, symlinks followed. |
| `bandpower` | output | `OUTPUT_DIR` (or `OUTPUT_DIR/<LABEL>/`) | `bandpower.csv`, `coherence.csv` (when `COHERENCE=true`), `provenance.json`. |

## Parameters

| Name | Default | What it does |
|---|---|---|
| `WINDOW` | `baseline` | Which part of each recording is used. `baseline`: for a recording with a known seizure, its baseline interval (`recordings.csv` `baseline_start_s`/`baseline_end_s`, or a `baseline` event); without one, from the start to `MARGIN_S` before the first onset. A recording with no known seizure is used whole. `all`: the whole recording. `interictal-only`: the whole recording minus every seizure from `MARGIN_S` before onset to `MARGIN_S` after its end (a seizure of unknown duration is cut to the end). `artifact` events are always cut out. |
| `MARGIN_S` | `30` | Seconds kept clear around seizures (see `WINDOW`). |
| `RESAMPLE_HZ` | `200` | Every recording is resampled (polyphase, anti-aliased) to this rate first, so 512/1024/2048 Hz recordings and the 200 Hz atlas are computed on the same grid. ≥ 160. |
| `WELCH_WINDOW_S` | `2` | Welch window length (Hamming). 2 s gives 0.5 Hz bins. |
| `WELCH_OVERLAP` | `0.5` | Window overlap as a fraction (0.5 = 1 s for 2 s windows). |
| `BANDPASS` | `0.5-80` | First-order Butterworth band-pass, zero-phase, before resampling; `none` skips it. Not applied again when every used channel's `_channels.tsv` `low_cutoff`/`high_cutoff` already lie within it (the HUP atlas clips are 0.5–80 Hz). |
| `BANDS` | `conventions` | `conventions`: delta 1–4, theta 4–8, alpha 8–13, beta 13–30, gamma 30–80 Hz. `bernabei2022`: the edges in the paper's Methods, 0.5–4, 4–8, 8–12, 12–30, 30–80 Hz. Or JSON, e.g. `{"delta": [1, 4], "gamma": [30, 80]}`. A `broadband` 1–80 Hz row is always added. |
| `COHERENCE` | `false` | Also write `coherence.csv`: median magnitude-squared coherence per band for every channel pair (same Welch windows). |
| `CHUNK_S` | `60` | Seconds read at a time; memory per chunk ≈ channels × seconds × rate × 8 bytes. |
| `LABEL` | (empty) | When set, outputs go to `OUTPUT_DIR/<LABEL>/` (a workflow running this step twice, e.g. baseline and interictal-only). |

## What it computes

Per recording (one at a time):

1. **Channels:** the `good` channels of `_channels.tsv` of type SEEG/ECOG/EEG/DBS (all signal
   channels when there is no `_channels.tsv`, minus ECG/EMG/EOG/trigger-like names).
2. **Interval:** chosen by `WINDOW` (above), in seconds, then cut to samples.
3. **Preprocessing:** first-order Butterworth band-pass 0.5–80 Hz (`scipy.signal.sosfiltfilt`),
   then `resample_poly` to 200 Hz. The data are streamed in `CHUNK_S` chunks with 4 s of
   padding on each side that absorbs filter and resampler transients and is cut off again;
   chunked and one-piece results agree to < 0.1 % (tested).
4. **Welch PSD:** 2 s symmetric Hamming windows, 1 s step, mean removed per window, one-sided
   density (µV²/Hz), averaged over every window of the interval (windows run across chunk
   boundaries; separate intervals, e.g. either side of an artifact, are not joined). Identical
   to `scipy.signal.welch` (tested).
5. **Per band:**
   - `abs_power`: the PSD integrated over the band (trapezoid, edges included), µV².
   - `rel_power`: `abs_power` / broadband (1–80 Hz) power: the band's share. The five shares sum to 1.
   - `norm_median`: **Bernabei et al. 2022's feature.** The Methods say: PSD from 0.5 to 80 Hz in
     0.5 Hz steps (Welch, 2 s Hamming window, 1 s overlap), "normalized the spectral density to
     have a sum of 1", then the median value within each band. So the PSD bins from 0.5 to 80 Hz
     are divided by their sum, and the median over the band's bins (lo ≤ f < hi; the top band
     includes 80 Hz) is reported. This measures the *shape* of the spectrum, insensitive to overall
     amplitude, and it is what `ieeg-normative` z-scores by default.

### Choices, and where the paper and the conventions differ

- **Bands.** `conventions.md` lists delta 1–4, theta 4–8, alpha 8–13, beta 13–30, gamma 30–80
  Hz "(Bernabei 2022)". These are the bands of Taylor et al. 2022 (Brain) and of the legend of
  Bernabei's Fig. 4; the text of Bernabei's Methods gives delta 0.5–4, alpha 8–12 and beta 12–30
  Hz. The default follows the contract (band names are what downstream steps key on, and every
  step should use the same edges); `BANDS=bernabei2022` reproduces the Methods text. On the atlas
  the choice changes ieeg-normative's AUCs by ≤ 0.03 (see processor-ieeg-normative's README).
- **Normalisation range.** `norm_median` uses 0.5–80 Hz as in the paper; `rel_power` uses 1–80 Hz
  as the conventions and the task asked. They are different quantities and both are reported.
- **Filter.** The paper used a first-order Butterworth 0.5–80 Hz band-pass and a 60 Hz notch;
  the notch is ieeg-prep's job. The paper does not say whether the filter was zero-phase; this
  app uses `sosfiltfilt` (zero-phase, which squares the magnitude response). For data whose
  `_channels.tsv` already says 0.5–80 Hz (the atlas) the filter is not applied twice.
- **Window of the recording.** The atlas used artifact-free wake clips ≥ 2 h from seizures. The
  demo clips are 2 min of pre-ictal baseline at most; pre-ictal activity is not interictal
  activity, so treat comparisons against the atlas with that in mind (or use longer interictal
  recordings with `WINDOW=all`).

## Outputs

**`bandpower.csv`**: one row per good channel and band (5 bands + `broadband`).

| Column | |
|---|---|
| `recording_id` | BIDS entity chain of the source recording (`sub-HUP224_ses-postimplant_task-ictal339135_run-00`) |
| `channel` | `_channels.tsv` name (e.g. `LA01-LA02`) |
| `band` | `delta`, `theta`, `alpha`, `beta`, `gamma`, `broadband` |
| `abs_power` | band power, µV² (PSD integrated over the band) |
| `rel_power` | share of the 1–80 Hz power (broadband = 1) |
| `norm_median` | median of the sum-1-normalised 0.5–80 Hz PSD over the band's bins (Bernabei 2022) |

`norm_median` is an extra column beyond `conventions.md` (which lists `abs_power, rel_power`).

**`coherence.csv`** (`COHERENCE=true`): `recording_id, channel_a, channel_b, band, coherence`
(median magnitude-squared coherence over the band's bins, 0–1), each unordered pair once.

**`provenance.json`**: app and version, tool versions, parameters, the spectral settings
(`spectral`, read by ieeg-normative to compute its reference identically), inputs found, and per
recording: status, channel counts, sampling rate, the window chosen and why, the intervals used
(s), seconds used, Welch windows, whether the band-pass was applied, notes, seconds; or the error.

A recording that fails (unreadable, no good channel, too short, sampling rate < 160 Hz) is logged
and skipped; the run fails only when every recording failed or none was found.

## Runtime

`computeTypes: [standard, lambda]`, 1 vCPU / 2 GiB, 900 s timeout. One recording at a time,
streamed: memory depends on channel count × `CHUNK_S`, not on recording length.

Measured in Docker (linux/amd64 **emulated on an M-series Mac**, `--cpus 2 --memory 2g`):

| Input | Window | Wall time | cgroup `memory.peak` |
|---|---|---|---|
| 4 demo clips of HUP224 (230 channels, 1024 Hz, 4.6–5.7 min each; plain BIDS) | `all` | 14.7 s | 0.90 GB |
| same, `COHERENCE=true` (26,335 pairs × 5 bands per clip) | `all` | 71 s | 0.94 GB |
| 5 atlas patients (22–74 channels, 200 Hz, 60 s) | `baseline` | 0.9 s | 0.33 GB |
| ieeg-prep output of the demo (24 clips, 6 patients, 72–206 channels, 1024 Hz) | `baseline` (90 s each) | 16 s (12 s in the app) | 0.91 GB |

The peak includes the page cache of the EDF files read; the process itself stays well below.
A 1 h, 128-channel, 1024 Hz recording is ~11× the samples of a demo clip at half the channels:
~25 s emulated (~100 s with coherence). The image is 124 MB compressed (275 MB unpacked).
For studies with many hour-long recordings and `COHERENCE=true`, run it as `standard`
(the 900 s timeout is Lambda's limit).

## Running it locally

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python      # pytest: synthetic fixture (derivative + plain BIDS), Lambda handler
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/, into data/output/
make run-demo DEMO=../workflow-epilepsy-ieeg/data/work/prep/ieeg
```

`tests/make_fixture.py` writes an ieeg-prep-style derivative (an ictal clip at 512 Hz with a
bad channel, a 10 Hz and a 2 Hz rhythm, a 20 Hz seizure burst after a baseline, an artifact
event; an interictal clip; an atlas-like 200 Hz subject with 0.5–80 Hz cutoffs) and a plain
BIDS-iEEG dataset with an `ictal` event and no baseline.

As a Lambda function (`processor/entry.py` hands over to `awslambdaric` when
`AWS_LAMBDA_RUNTIME_API` is set; `processor/handler.py` takes `inputDir`, `outputDir` and the
parameters from the payload and runs the same `main()`; a failed run raises):

```sh
curl -XPOST localhost:9000/2015-03-31/functions/function/invocations \
  -d '{"inputDir": "/in", "outputDir": "/tmp/out", "WINDOW": "baseline", "COHERENCE": false}'
```

## Demo data result

Over the ieeg-prep derivative of the six demo patients (Discover 530; 24 clips, 2882 bipolar
channels, all good), default parameters: every clip has a baseline interval in
`recordings.csv` (0–90 s, ending ≥ 10 s before the expert onset), so 90 s (89 Welch windows) per
clip were used; 17,292 rows. Median over channels of `rel_power`: delta 0.414, theta 0.186,
alpha 0.148, beta 0.138, gamma 0.016; median broadband power 386 µV².

The clips carry the experts' consensus onset channels (`onset` in `_channels.tsv`, 119 of 2882
channel-clips). In the pre-ictal baseline those channels already have a larger delta share
(median 0.58 vs 0.41; AUC 0.71 for `rel_power`, 0.74 for `norm_median`), a smaller alpha and
beta share (AUC 0.29, 0.37) and more absolute power in every band (AUC 0.74–0.83). These are
pooled over patients and seizures of the same patients, so they are descriptive, not a
validated localiser.

## Limitations

- Band power describes the spectrum of the chosen interval only. A pre-ictal baseline (the
  default on seizure clips) is not interictal background; sleep stage, medication, state of
  vigilance and recent seizures all change spectra (Wang et al. 2023 show interictal abnormality
  varies over hours).
- Bipolar and referential montages give different spectra; compare like with like (the HUP atlas
  is bipolar).
- `norm_median` and `rel_power` describe spectral shape; they are blind to a channel that is
  uniformly louder or quieter. `abs_power` is not, but depends on electrode type, impedance and
  montage.
- Coherence between neighbouring bipolar channels that share a contact is inflated by the
  shared contact.
- Line noise is not removed here (ieeg-prep notches it); with `BANDPASS=none` and raw input,
  the gamma band includes 50/60 Hz.

## Licences

Apache-2.0 for this app's code (`LICENSE`). The image bundles MNE-Python, NumPy, SciPy (BSD-3-
Clause) and their dependencies; see `THIRD_PARTY.md`. Nothing GPL is linked in-process.

## Citations

Method (features, preprocessing, normalisation):

> Bernabei JM, Sinha N, Arnold TC, Conrad E, Ong I, Pattnaik AR, Stein JM, Shinohara RT, Lucas TH,
> Bassett DS, Davis KA, Litt B. **Normative intracranial EEG maps epileptogenic tissues in focal
> epilepsy.** Brain 145(6):1949–1961 (2022). doi:10.1093/brain/awab480, PMID 35640886.

Band edges of the conventions (delta 1–4, theta 4–8, alpha 8–13, beta 13–30, gamma 30–80 Hz):

> Taylor PN, Papasavvas CA, Owen TW, Schroeder GM, Hutchings FE, Chowdhury FA, Diehl B, Duncan JS,
> McEvoy AW, Miserocchi A, de Tisi J, Vos SB, Walker MC, Wang Y. **Normative brain mapping of
> interictal intracranial EEG to localize epileptogenic tissue.** Brain 145(3):939–949 (2022).
> doi:10.1093/brain/awab380, PMID 35075485.

Welch's method:

> Welch P. **The use of fast Fourier transform for the estimation of power spectra: A method based
> on time averaging over short, modified periodograms.** IEEE Trans Audio Electroacoust
> 15(2):70–73 (1967). doi:10.1109/TAU.1967.1161901.

Software and data format:

> Gramfort A, Luessi M, Larson E, Engemann DA, Strohmeier D, Brodbeck C, Goj R, Jas M, Brooks T,
> Parkkonen L, Hämäläinen M. **MEG and EEG data analysis with MNE-Python.** Front Neurosci 7:267
> (2013). doi:10.3389/fnins.2013.00267, PMID 24431986. (Crossref lists only the first author;
> `CITATION.cff` follows Crossref, the full list here is from PubMed.)

> Virtanen P, Gommers R, Oliphant TE, et al. (SciPy 1.0 Contributors). **SciPy 1.0: fundamental
> algorithms for scientific computing in Python.** Nat Methods 17(3):261–272 (2020).
> doi:10.1038/s41592-019-0686-2, PMID 32015543.

> Holdgraf C, Appelhoff S, Bickel S, et al. **iEEG-BIDS, extending the Brain Imaging Data Structure
> specification to human intracranial electrophysiology.** Sci Data 6:102 (2019).
> doi:10.1038/s41597-019-0105-7, PMID 31239438.

Every DOI was checked against Crossref (title, authors, year, volume, pages) and every PMID
against PubMed E-utilities; `CITATION.cff` lists all authors in Crossref's order.
