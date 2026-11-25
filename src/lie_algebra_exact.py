"""Exact-Q adjoint-irrep variance laws, SO(2n) Majorana Hodge-chiral selection rules,
U(1) Hamming-weight block formulas, and SU(2^n) unstructured barren-plateau bounds.

ZERO third-party imports, ZERO floating-point arithmetic: pure Python stdlib Fractions (Q).
Author: Lluis Eriksson
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Dict, List, Sequence, Tuple


def so2n_majorana_irrep_dim(n: int, k: int, chirality: str = "full") -> int:
    """Exact dimension of the SO(2n) adjoint irrep Lambda^k (or chiral Lambda^n_+/- for even n, k=n)."""
    if n < 2 or not (1 <= k <= 2 * n - 1):
        raise ValueError(f"Invalid SO(2n) parameters n={n}, k={k}")
    full_dim = math.comb(2 * n, k)
    if chirality == "full":
        return full_dim
    if chirality in ("plus", "minus"):
        if k != n or (n % 2 != 0) or n < 4:
            raise ValueError(
                f"Universal Hodge-chiral equipartition requires even n >= 4 and k == n "
                f"(so(4) = su(2)_+ + su(2)_- is not simple); got n={n}, k={k}"
            )
        return full_dim // 2
    raise ValueError(f"Unknown chirality {chirality!r}")


def so2n_casimir_ad_norm_sq(n: int, k: int, gate_energy: Fraction, chirality: str = "full") -> Fraction:
    """Exact rational ||ad_H|_{V_lambda}||_HS^2 for H = sum_{a<b} h_{ab} (i/2) gamma_a gamma_b
    with Frobenius gate energy E_H = sum_{a<b} h_{ab}^2:
      - Full Lambda^k:  2 * E_H * C(2n - 2, k - 1)
      - Chiral Lambda^n_+/- (even n, k = n): E_H * C(2n - 2, n - 1)  (half of full!)
    """
    if gate_energy < 0:
        raise ValueError(f"gate_energy must be >= 0, got {gate_energy}")
    full_cas = Fraction(2 * math.comb(2 * n - 2, k - 1), 1) * gate_energy
    if chirality == "full":
        return full_cas
    if chirality in ("plus", "minus"):
        if k != n or (n % 2 != 0) or n < 4:
            raise ValueError(
                f"Universal Hodge-chiral equipartition requires even n >= 4 and k == n "
                f"(so(4) splits into commuting ideals); got n={n}, k={k}"
            )
        return full_cas / 2
    raise ValueError(f"Unknown chirality {chirality!r}")


def so2n_sector_variance_exact(
    n: int,
    k: int,
    gate_energy: Fraction,
    rho_norm_sq: Fraction,
    M_norm_sq: Fraction,
    chirality: str = "full",
) -> Tuple[int, Fraction, Fraction]:
    """Compute exact (dim_V, casimir_ad_norm_sq, variance_contribution) in Q for an SO(2n) sector.

    Here rho_norm_sq = ||Pi_{V}(rho)||_HS^2 and M_norm_sq = ||Pi_{V}(M)||_HS^2 in Hilbert-Schmidt norm
    Tr(X^2) (so for a Pauli/Majorana operator M with M^2 = I, ||M||_HS^2 = 2^n; for a single basis
    projector on a pure product state in Lambda^k, ||Pi_k(rho)||_HS^2 <= 1).
    """
    if rho_norm_sq < 0 or M_norm_sq < 0:
        raise ValueError("Norms squared must be nonnegative")
    dim_v = so2n_majorana_irrep_dim(n, k, chirality)
    cas = so2n_casimir_ad_norm_sq(n, k, gate_energy, chirality)
    var = cas * rho_norm_sq * M_norm_sq / Fraction(dim_v * dim_v, 1)
    return dim_v, cas, var


def u1_sector_variance_exact(
    n: int,
    w: int,
    H_traceless_norm_sq: Fraction,
    rho_traceless_norm_sq: Fraction,
    M_traceless_norm_sq: Fraction,
) -> Tuple[int, Fraction, Fraction]:
    """Compute exact (dim_su_dw, casimir_ad_norm_sq, variance_contribution) in Q for a U(1)
    Hamming-weight sector w (Hilbert block dim d_w = C(n, w), adjoint irrep su(d_w) dim d_w^2 - 1):
      ||ad_H|_{su(d_w)}||_HS^2 = 2 * d_w * ||H_w^o||_HS^2,
      Var_w = 2 * d_w * ||H_w^o||_HS^2 * ||rho_w^o||_HS^2 * ||M_w^o||_HS^2 / (d_w^2 - 1)^2.
    """
    if n < 2 or not (0 <= w <= n):
        raise ValueError(f"Invalid U(1) parameters n={n}, w={w}")
    if H_traceless_norm_sq < 0 or rho_traceless_norm_sq < 0 or M_traceless_norm_sq < 0:
        raise ValueError("Traceless norms squared must be nonnegative")
    dw = math.comb(n, w)
    if dw < 2:
        if H_traceless_norm_sq != 0 or rho_traceless_norm_sq != 0 or M_traceless_norm_sq != 0:
            raise ValueError(
                f"1D U(1) block (n={n}, w={w}, d_w={dw}) has zero traceless subspace; "
                "all traceless norms must be 0"
            )
        return 0, Fraction(0, 1), Fraction(0, 1)
    dim_adj = dw * dw - 1
    cas = Fraction(2 * dw, 1) * H_traceless_norm_sq
    var = cas * rho_traceless_norm_sq * M_traceless_norm_sq / Fraction(dim_adj * dim_adj, 1)
    return dim_adj, cas, var


def su_unstructured_variance_exact(
    n: int,
    H_traceless_norm_sq: Fraction,
    rho_traceless_norm_sq: Fraction,
    M_traceless_norm_sq: Fraction,
) -> Tuple[int, Fraction, Fraction]:
    """Exact rational variance for unstructured SU(2^n) 2-design circuits (d = 2^n):
      Var = 2 * d * ||H^o||_HS^2 * ||rho^o||_HS^2 * ||M^o||_HS^2 / (d^2 - 1)^2.
    """
    if n < 1:
        raise ValueError(f"Invalid qubit count n={n}")
    if H_traceless_norm_sq < 0 or rho_traceless_norm_sq < 0 or M_traceless_norm_sq < 0:
        raise ValueError("Traceless norms squared must be nonnegative")
    d = 1 << n
    dim_adj = d * d - 1
    cas = Fraction(2 * d, 1) * H_traceless_norm_sq
    var = cas * rho_traceless_norm_sq * M_traceless_norm_sq / Fraction(dim_adj * dim_adj, 1)
    return dim_adj, cas, var


def exact_sample_complexity_lower_bound(tau: Fraction, fourth_moment_floor: Fraction) -> Fraction:
    """Theorem 3.4 exact Chebyshev lower bound N >= 4 * B^2 / tau^2 to certify Var >= tau."""
    if tau <= 0 or fourth_moment_floor <= 0:
        raise ValueError("tau and fourth_moment_floor must be strictly positive")
    return Fraction(4, 1) * fourth_moment_floor * fourth_moment_floor / (tau * tau)
