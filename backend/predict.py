"""Dependency-light predictor for the EvoCult weighted ordinal model.

Loads evocult_weighted_model.json (exported by Cell 14 of the notebook) and runs
the forward pass in pure numpy — no sklearn or pickle version coupling. Verified in
the notebook to reproduce the sklearn pipeline's predictions exactly.
"""
import json
import os
import numpy as np

_MODEL_PATH = os.path.join(os.path.dirname(__file__), "evocult_weighted_model.json")


class EvoCultModel:
    def __init__(self, path=_MODEL_PATH):
        with open(path) as fh:
            self.b = json.load(fh)
        self.features = self.b["kept_features"]
        self.media = self.b["media"]
        self.class_order = self.b["class_order"]

    def _forward(self, mp, x_raw):
        x = np.asarray(x_raw, dtype=float)
        med = np.asarray(mp["impute_median"], dtype=float)
        x = np.where(np.isnan(x), med, x)                       # median imputation
        x = (x - np.asarray(mp["scaler_center"])) / np.asarray(mp["scaler_scale"])
        x = x * np.asarray(mp["weights"])                       # feature weighting
        cum = []
        for bnd in mp["boundaries"]:                            # cumulative logits
            if bnd["type"] == "const":
                cum.append(bnd["value"])
            else:
                z = float(np.dot(x, np.asarray(bnd["coef"])) + bnd["intercept"])
                cum.append(1.0 / (1.0 + np.exp(-z)))
        k = len(mp["classes"])
        P = np.zeros(k)
        prev = 1.0
        for j in range(k - 1):
            P[j] = prev - cum[j]
            prev = cum[j]
        P[-1] = prev
        P = np.clip(P, 1e-9, 1.0)
        return int(np.argmax(P)), P

    def predict(self, feature_row):
        """feature_row: dict {feature_name: value or None}. Returns per-medium dict."""
        x = [feature_row.get(f, np.nan) for f in self.features]
        x = [np.nan if v is None else v for v in x]
        out = {}
        for medium in self.media:
            idx, proba = self._forward(self.b["models"][medium], x)
            out[medium] = {
                "class": self.class_order[idx],
                "proba": {c: round(float(p), 4) for c, p in zip(self.class_order, proba)},
            }
        return out


if __name__ == "__main__":
    # smoke test: an all-missing row still predicts (falls back to per-medium medians)
    m = EvoCultModel()
    print(f"loaded {len(m.features)} features, media={m.media}")
    print(m.predict({}))
