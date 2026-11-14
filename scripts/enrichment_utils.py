"""Shared enrichment statistics."""

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def bh_group(values: pd.Series) -> pd.Series:
    p = values.astype(float).to_numpy()
    order = np.argsort(np.nan_to_num(p, nan=1.0))
    ranked = np.nan_to_num(p[order], nan=1.0)
    m = len(ranked)
    q = np.empty_like(ranked)
    running = 1.0
    for i in range(m - 1, -1, -1):
        running = min(running, ranked[i] * m / (i + 1))
        q[i] = running
    out = np.empty_like(q)
    out[order] = q
    return pd.Series(out, index=values.index)


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    return float(hypergeom.sf(k - 1, N, K, n))


def fisher_greater(a: int, b: int, c: int, d: int) -> float:
    return hypergeom_sf(a, a + b + c + d, a + c, a + b)


def partial_conjunction(p_values: list[float]) -> float:
    """Bonferroni partial conjunction: enrichment in at least two species.

    Valid without assuming independence among the related species.
    """
    ps = sorted([float(p) for p in p_values] + [1.0] * (3 - len(p_values)))
    if len(ps) < 2:
        return 1.0
    return min(1.0, (len(ps) - 1) * ps[1])
