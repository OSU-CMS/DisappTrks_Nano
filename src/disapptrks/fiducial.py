"""Numerical core of the electron and muon fiducial-map construction."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HotSpot:
    eta: float
    phi: float
    radius: float
    sigma: float


@dataclass(frozen=True)
class FiducialMapSummary:
    mean_inefficiency: float
    stddev_inefficiency: float
    inefficiency: np.ndarray
    significance: np.ndarray
    hot_spots: tuple[HotSpot, ...]


def summarize_fiducial_map(
    before: np.ndarray,
    after: np.ndarray,
    eta_edges: np.ndarray,
    phi_edges: np.ndarray,
    threshold: float = 2.0,
) -> FiducialMapSummary:
    """Reproduce the legacy bin-by-bin inefficiency hot-spot calculation."""
    before = np.asarray(before, dtype=float)
    after = np.asarray(after, dtype=float)
    eta_edges = np.asarray(eta_edges, dtype=float)
    phi_edges = np.asarray(phi_edges, dtype=float)

    if before.shape != after.shape:
        raise ValueError("before and after histograms must have identical shapes")
    if before.shape != (len(eta_edges) - 1, len(phi_edges) - 1):
        raise ValueError("histogram shape does not match the supplied axes")

    occupied = before > 0.0
    if not np.any(occupied):
        raise ValueError("before histogram has no occupied bins")

    mean = float(after[occupied].sum() / before[occupied].sum())
    inefficiency = np.zeros_like(after)
    inefficiency[occupied] = after[occupied] / before[occupied]

    n_occupied = int(np.count_nonzero(occupied))
    stddev = (
        float(
            np.sqrt(
                np.sum((inefficiency[occupied] - mean) ** 2)
                / (n_occupied - 1)
            )
        )
        if n_occupied > 1
        else 0.0
    )
    significance = np.zeros_like(after)
    if stddev > 0.0:
        significance[occupied] = (inefficiency[occupied] - mean) / stddev

    hot_spots = []
    for ix, iy in np.argwhere(
        occupied & (inefficiency > 0.0) & (significance > threshold)
    ):
        eta_width = eta_edges[ix + 1] - eta_edges[ix]
        phi_width = phi_edges[iy + 1] - phi_edges[iy]
        hot_spots.append(
            HotSpot(
                eta=float((eta_edges[ix] + eta_edges[ix + 1]) / 2.0),
                phi=float((phi_edges[iy] + phi_edges[iy + 1]) / 2.0),
                radius=float(np.hypot(eta_width / 2.0, phi_width / 2.0)),
                sigma=float(significance[ix, iy]),
            )
        )

    return FiducialMapSummary(
        mean_inefficiency=mean,
        stddev_inefficiency=stddev,
        inefficiency=inefficiency,
        significance=significance,
        hot_spots=tuple(hot_spots),
    )
