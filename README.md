# EvoCult — cross-species ligand compatibility predictor

A web tool that takes a target species' ligand and receptor sequences and predicts,
per organoid medium (ICH / ICM / WEN), whether human-derived niche factors will support
**Never**, **Short-term**, or **Long-term** maintenance — and flags the specific
ligand–receptor pairs that break compatibility (the candidates for a rescue cocktail),
exactly as we do for zebrafish.

The repo has two halves:

| Part | Runs where | What it does |
|------|-----------|--------------|
| **Static prototype** (`index.html`, `assets/`) | GitHub Pages — no server | Computes the sequence-identity (`X1`) axis in the browser and scores each pair against a compatibility floor. Works today. |
| **Backend** (`backend/`) | Any Python host | Runs the full trained model — the importance-weighted, covariance-eliminated ordinal Frank–Hall classifier from the notebook — over all five feature axes. The production path. |

---

## Deploy the prototype on GitHub Pages

1. Push this folder to a GitHub repo.
2. **Settings → Pages → Build and deployment → Source: Deploy from a branch**, pick
   `main` / root.
3. Open `https://<user>.github.io/<repo>/`.

That's it — the prototype is fully static (no build step, no backend needed). The
`.nojekyll` file keeps GitHub Pages from touching the `assets/` folder.

### Using it
- **Load zebrafish example** — scores *Danio rerio* from the real precomputed `X1`
  values in the training set. EGF·EGFR and RSPO1·LGR5 come back below the floor: the
  known teleost bottlenecks.
- **Your own species** — paste or upload FASTA with the ligand and receptor orthologs.
  Headers are matched to the core niche factors by name (`>WNT3A`, `>FZD8`, `>EGF`, …).

---

## What the prototype computes (and what it doesn't)

The browser scorer uses only the **`X1` sequence-identity axis** — the cheapest, fastest
feature in the pipeline — aligning each ortholog to a reference and scoring each pair by
its limiting identity against an illustrative compatibility floor
(`Never < 0.55 ≤ Short-term < 0.75 ≤ Long-term`).

> The reference sequences in `assets/references.js` are **illustrative placeholders** so
> the demo runs out of the box. Replace them with the curated human UniProt interface
> panel, or move identity computation to the backend. The zebrafish example uses real
> values and is unaffected.

The full model additionally uses structure (`X2`), docking energetics (`X3`), biophysical
(`X4`) and evolutionary-context (`X5`) features. Those need the heavy pipeline
(AlphaFold, PRODIGY, PAML/HyPhy, expression data), so the prototype leaves them to the
backend.

---

## Run the backend (full model)

```bash
cd backend
pip install -r requirements.txt
python app.py            # http://127.0.0.1:8000
```

```bash
curl -X POST http://127.0.0.1:8000/api/score \
  -H "Content-Type: application/json" \
  -d '{"species":"Xenopus laevis","fasta":">WNT3A\nMAPL...\n>FZD8\nMEWG..."}'
```

To make the frontend call the backend instead of the in-browser proxy, replace the body
of `runFromFasta()` in `assets/app.js` with a `fetch("http://<host>/api/score", …)` and
render the returned `media` + `pairs`. (GitHub Pages can't host the Python process
itself; deploy the backend anywhere that runs Python and point the frontend at it.)

### Backend layout
- `predict.py` — pure-numpy forward pass; loads `evocult_weighted_model.json`. No sklearn
  or pickle at serve time. Verified in the notebook to match the sklearn pipeline exactly.
- `compute_features.py` — FASTA → feature row. `X1` is real; `X2`–`X5` are `TODO` hooks
  into the pipeline. The model imputes anything still missing.
- `app.py` — the `/api/score` and `/api/health` Flask endpoints.
- `evocult_weighted_model.json` — the trained model parameters (see below).

---

## Where the model comes from

`evocult_weighted_model.json` is exported by **Cell 14** of
`EvoCult_Feature_Evaluation.ipynb`. That cell:
1. takes the covariance-eliminated feature set (Cell 12),
2. weights each feature by its univariate `LOSO_QWK` importance (Cell 13),
3. fits one ordinal model per medium on all training species, and
4. writes the imputer medians, scaler, weights and per-boundary logistic coefficients.

Re-run the notebook to refresh the model; the file drops straight into `backend/`.

---

## Structure

```
evocult-platform/
├── index.html                    # static app (GitHub Pages)
├── .nojekyll
├── assets/
│   ├── style.css
│   ├── app.js                    # FASTA parse · alignment · scoring · readout
│   └── references.js             # niche factors, pairs, floors, zebrafish X1 values
└── backend/
    ├── app.py                    # Flask /api/score
    ├── predict.py                # numpy model forward pass
    ├── compute_features.py       # FASTA → features (X1 real, X2–X5 stubbed)
    ├── requirements.txt
    └── evocult_weighted_model.json
```
