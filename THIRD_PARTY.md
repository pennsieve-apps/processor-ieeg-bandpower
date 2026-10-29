# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components, installed from PyPI as binary wheels at the versions pinned in `requirements.txt`.
All are permissively licensed (MPL-2.0 for certifi and part of tqdm, file-level copyleft that
does not extend to this app); nothing GPL is linked in-process.

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| MNE-Python | 1.13.2 | Reads EDF/BDF, BrainVision, FIF and EEGLAB recordings (`mne.io.read_raw`), channel data in µV | BSD-3-Clause | https://mne.tools |
| NumPy | 2.5.4 | Arrays, FFT of the Welch windows | BSD-3-Clause (bundled parts 0BSD, MIT, Zlib, CC0-1.0) | https://numpy.org |
| SciPy | 1.18.1 | Butterworth band-pass (`sosfiltfilt`), polyphase resampling (`resample_poly`), Hamming window | BSD-3-Clause | https://scipy.org |
| awslambdaric | 4.2.3 | Lambda runtime interface client (only when run as a Lambda function) | Apache-2.0 | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| simplejson, snapshot-restore-py | 4.2.0, 1.0.0 | awslambdaric dependencies | MIT OR AFL-2.1; Apache-2.0 | PyPI |
| Matplotlib | 3.11.2 | MNE dependency (not used for output) | Matplotlib License (PSF-based, BSD-compatible) | https://matplotlib.org |
| contourpy, cycler, fonttools, kiwisolver, pyparsing, pillow | 1.4.0, 0.12.1, 4.68.0, 1.5.1, 3.3.3, 13.0.0 | Matplotlib dependencies | BSD-3-Clause; BSD-3-Clause; MIT; BSD-3-Clause; MIT; MIT-CMU | PyPI |
| decorator, lazy-loader, pooch, tqdm, Jinja2, MarkupSafe, packaging | 5.3.1, 0.6, 1.9.0, 4.70.1, 3.1.6, 3.0.4, 26.3 | MNE dependencies | BSD-2-Clause; BSD-3-Clause; BSD-3-Clause; MPL-2.0 AND MIT; BSD-3-Clause; BSD-3-Clause; Apache-2.0 OR BSD-2-Clause | PyPI |
| platformdirs, requests, certifi, charset-normalizer, idna, urllib3 | 4.14.0, 2.34.2, 2026.10.23, 3.5.2, 3.20, 2.8.0 | pooch dependencies (no network use at run time) | MIT; Apache-2.0; MPL-2.0; MIT; BSD-3-Clause; MIT | PyPI |
| python-dateutil, six, tzdata | 3.0.0, 1.17.0, 2026.6 | Matplotlib dependencies | Apache-2.0/BSD dual; MIT; Apache-2.0 | PyPI |
| python:3.12.15-slim-bookworm | — | Base image (Debian 12 + CPython 3.12, pip) | PSF-2.0 (Python); Debian packages under their own free licences; pip MIT | https://hub.docker.com/_/python |

The method follows Bernabei et al. 2022 (Brain, doi:10.1093/brain/awab480) and was
reimplemented from the paper; no code was taken from the authors' repository
(jbernabei/iEEG_atlas, which has no licence).

Demo data used for testing (not in the image or the repository): Discover 530 "Multi-Expert
Seizure Annotation" (doi:10.26275/dely-fzwx, CC BY 4.0) through ieeg-prep, and Discover 179
"HUP Interictal iEEG Atlas" (doi:10.26275/lzhy-ox20, ODC-By 1.0).
