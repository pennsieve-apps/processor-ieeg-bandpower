# sample-input

Real input for this app is a BIDS-iEEG dataset, ideally the `ieeg` output of
**processor-ieeg-prep** (bipolar, notch-filtered, with `recordings.csv`, baseline intervals and
labels). The workflows' demo data are public on Pennsieve Discover and fetched by
`workflow-epilepsy-ieeg/scripts/fetch-sample-data.sh`:

> Ojemann W, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. Multi-Expert Seizure Annotation.
> Pennsieve Discover 530, doi:10.26275/dely-fzwx (CC BY 4.0). SEEG seizure clips, 1024 Hz,
> 120 s before the clinical onset.

```sh
cd ../workflow-epilepsy-ieeg && scripts/fetch-sample-data.sh data      # ~2.7 GB
# run ieeg-prep into data/work/prep/, then:
make run-demo DEMO=../workflow-epilepsy-ieeg/data/work/prep/ieeg      # -> data/output-demo/
```

The raw dataset also works directly (plain BIDS: `*_ieeg.edf` with `_channels.tsv` and an
`ictal`/`stim` event; without a baseline, the start of the clip up to `MARGIN_S` before the onset is used):

```sh
make run-demo DEMO=../workflow-epilepsy-ieeg/data/seizures/DATA
```

For tests and `make run`, `tests/make_fixture.py` writes a small synthetic fixture instead
(`make fixture` -> `data/fixture/linked/ieeg/` and `data/fixture/plain/`).
