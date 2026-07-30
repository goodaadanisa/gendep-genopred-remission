"""Dependency-light binary-prediction metrics used for independent aggregation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit
from scipy.stats import rankdata


@dataclass(frozen=True)
class Calibration:
    """Represent Calibration as an immutable workflow record."""
    intercept: float
    slope: float


def _arrays(outcome: np.ndarray, prediction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(outcome, dtype=float)
    p = np.asarray(prediction, dtype=float)
    if y.ndim != 1 or p.ndim != 1 or y.size != p.size:
        raise ValueError("Outcome and prediction must be aligned one-dimensional arrays")
    if y.size == 0 or not np.isfinite(y).all() or not np.isfinite(p).all():
        raise ValueError("Outcome and prediction must be finite and non-empty")
    if not set(np.unique(y)).issubset({0.0, 1.0}):
        raise ValueError("Outcome must be binary")
    return y, p


def auc(outcome: np.ndarray, prediction: np.ndarray) -> float:
    """Calculate ROC AUC using average ranks so tied predictions are handled correctly."""
    y, p = _arrays(outcome, prediction)
    positives = y == 1
    n_positive = int(positives.sum())
    n_negative = int((~positives).sum())
    if n_positive == 0 or n_negative == 0:
        raise ValueError("AUC requires both outcome classes")
    ranks = rankdata(p, method="average")
    return float(
        (ranks[positives].sum() - n_positive * (n_positive + 1) / 2)
        / (n_positive * n_negative)
    )


def brier(outcome: np.ndarray, prediction: np.ndarray) -> float:
    """Calculate the mean squared error of predicted event probabilities."""
    y, p = _arrays(outcome, prediction)
    return float(np.mean((p - y) ** 2))


def log_loss(outcome: np.ndarray, prediction: np.ndarray, epsilon: float = 1e-15) -> float:
    """Calculate binary logarithmic loss after bounding probabilities away from zero and one."""
    y, p = _arrays(outcome, prediction)
    bounded = np.clip(p, epsilon, 1 - epsilon)
    return float(-np.mean(y * np.log(bounded) + (1 - y) * np.log(1 - bounded)))


def _newton_logistic(
    design: np.ndarray,
    outcome: np.ndarray,
    *,
    offset: np.ndarray | None = None,
    max_iter: int = 100,
    tolerance: float = 1e-12,
) -> np.ndarray:
    beta = np.zeros(design.shape[1], dtype=float)
    base = np.zeros(outcome.size, dtype=float) if offset is None else np.asarray(offset, dtype=float)
    for _ in range(max_iter):
        eta = base + design @ beta
        probability = expit(eta)
        weights = np.clip(probability * (1 - probability), 1e-12, None)
        score = design.T @ (outcome - probability)
        information = (design.T * weights) @ design
        step = np.linalg.pinv(information) @ score
        beta_next = beta + step
        if np.max(np.abs(step)) <= tolerance:
            return beta_next
        beta = beta_next
    return beta


def calibration(
    outcome: np.ndarray,
    prediction: np.ndarray,
    epsilon: float = 1e-15,
) -> Calibration:
    """Estimate calibration-in-the-large and calibration slope from fixed predictions."""
    y, p = _arrays(outcome, prediction)
    bounded = np.clip(p, epsilon, 1 - epsilon)
    logits = np.log(bounded / (1 - bounded))
    intercept = _newton_logistic(np.ones((y.size, 1)), y, offset=logits)[0]
    slope_beta = _newton_logistic(np.column_stack([np.ones(y.size), logits]), y)
    return Calibration(intercept=float(intercept), slope=float(slope_beta[1]))


def metric_bundle(outcome: np.ndarray, prediction: np.ndarray, epsilon: float = 1e-15) -> dict[str, float]:
    """Return discrimination, overall-accuracy and calibration metrics for one prediction vector."""
    cal = calibration(outcome, prediction, epsilon)
    return {
        "AUC": auc(outcome, prediction),
        "Brier": brier(outcome, prediction),
        "log_loss": log_loss(outcome, prediction, epsilon),
        "calibration_intercept": cal.intercept,
        "calibration_slope": cal.slope,
    }
