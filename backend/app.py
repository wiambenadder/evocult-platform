"""EvoCult scoring API (production path for the web platform).

    POST /api/score
      body: {"species": "Xenopus laevis", "fasta": ">WNT3A\\nMAPL...\\n>FZD8\\n..."}
      returns: per-medium ordinal prediction + per-pair compatibility from the
               importance-weighted covariance-eliminated model trained in the notebook.

Run locally:
    pip install -r requirements.txt
    python app.py                      # serves http://127.0.0.1:8000
Then point the frontend's scorer at this endpoint instead of the in-browser proxy.
"""
from flask import Flask, request, jsonify
from flask_cors import CORS

from predict import EvoCultModel
from compute_features import features_from_fasta

app = Flask(__name__)
CORS(app)
MODEL = EvoCultModel()

# ligand · receptor pairs, mirrored from the frontend
PAIRS = [("WNT3A", "FZD8"), ("EGF", "EGFR"), ("RSPO1", "LGR5"),
         ("NOG", "BMPR1A"), ("NOG", "BMPR1B")]
FLOOR_NEVER, FLOOR_LONG = 0.55, 0.75


def _classify(v):
    if v is None:
        return None
    if v < FLOOR_NEVER:
        return "Never"
    return "Long-term" if v >= FLOOR_LONG else "Short-term"


@app.get("/api/health")
def health():
    return jsonify(status="ok", features=len(MODEL.features), media=MODEL.media)


@app.post("/api/score")
def score():
    data = request.get_json(force=True) or {}
    fasta = data.get("fasta", "")
    if not fasta.strip():
        return jsonify(error="No FASTA provided."), 400

    row = features_from_fasta(fasta)                       # X1 real, X2-X5 imputed
    media_pred = MODEL.predict(row)                        # trained weighted model

    # per-pair readout from the X1 identities we computed (compatibility floor)
    pairs = []
    for lig, rec in PAIRS:
        li, ri = row.get(f"X1_{lig}"), row.get(f"X1_{rec}")
        if li is None or ri is None:
            continue
        limiting = min(li, ri)
        pairs.append({"ligand": lig, "receptor": rec, "identity": round(limiting, 3),
                      "class": _classify(limiting)})

    return jsonify(species=data.get("species", ""), media=media_pred, pairs=pairs)


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8000, debug=True)
