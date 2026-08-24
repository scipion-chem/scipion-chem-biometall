class _PlattCalibratedPipeline:
    """RF pipeline wrapped with Platt logistic regression for probability calibration."""
    def __init__(self, base_pipeline, platt_lr):
        self.base_pipeline = base_pipeline
        self.platt_lr = platt_lr

    def predict_proba(self, X):
        raw = self.base_pipeline.predict_proba(X)[:, 1].reshape(-1, 1)
        return self.platt_lr.predict_proba(raw)

    def predict(self, X):
        return (self.predict_proba(X)[:, 1] >= 0.5).astype(int)
