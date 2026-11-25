r"""Zero-floating-point-trust exact-Q certificate builder and verifier for Paper 09.

Schema: LIE-TRAIN-CERT-1.0
Author: Lluis Eriksson

Security & Soundness Guarantees:
  - ZERO third-party imports: pure Python stdlib (`fractions.Fraction`, `math.comb`, `re`).
  - Strict recursive type whitelisting (`type(x) is ...`): rejects floats, bools, Decimals,
    tuples, subclasses, and non-str dict keys.
  - Two-regex anchored ASCII canonical rational parser: rejects `-0`, `2/1`, `4/14`, `+1`,
    leading zeros, whitespace, trailing newlines, and Unicode digits.
  - Exact key-set enforcement at EVERY dictionary level.
  - Replays exact binomial dimensions, exact Casimir adjoint norms, exact physical norm
    upper bounds (||rho||_HS^2 <= 1, ||Pi_lambda(rho)||_HS^2 <= 1, non-duplicate sectors),
    exact total rational variance Var in Q, and exact threshold comparison across 3 mutually
    exclusive verdict modes (soundness formalized in Lean 4 Theorems 14-16):
      1. `CERTIFIED_TRAINABLE`: exact Var >= tau (with tau > 0).
      2. `CERTIFIED_BARREN_PLATEAU`: 0 < Var <= epsilon (with epsilon < tau).
      3. `CERTIFIED_BLIND_SPOT`: exact Var == 0 (via chiral selection or scalar-gate commutant).
"""

from __future__ import annotations

import math
import re
from fractions import Fraction
from typing import Any, Dict, List, Sequence, Tuple

from .lie_algebra_exact import (
    so2n_sector_variance_exact,
    su_unstructured_variance_exact,
    u1_sector_variance_exact,
)

SCHEMA_VERSION = "LIE-TRAIN-CERT-1.0"
_CERT_ID_RE = re.compile(r"\A[A-Za-z0-9._-]{1,120}\Z", re.ASCII)
_INT_RE = re.compile(r"\A(?:0|-?[1-9][0-9]*)\Z", re.ASCII)
_FRAC_RE = re.compile(r"\A(?P<num>-?[1-9][0-9]*)/(?P<den>[1-9][0-9]*)\Z", re.ASCII)

_TOP_KEYS = frozenset({
    "schema_version",
    "cert_id",
    "family",
    "n_qubits",
    "dla_dimension",
    "threshold_tau",
    "verdict",
    "exact_variance",
    "sectors",
    "physical_checks",
})

_SO2N_SECTOR_KEYS = frozenset({
    "k_degree",
    "chirality",
    "irrep_dim",
    "gate_energy",
    "casimir_ad_norm_sq",
    "rho_proj_norm_sq",
    "M_proj_norm_sq",
    "sector_variance",
})

_U1_SECTOR_KEYS = frozenset({
    "hamming_weight",
    "block_dim",
    "irrep_dim",
    "H_traceless_norm_sq",
    "casimir_ad_norm_sq",
    "rho_traceless_norm_sq",
    "M_traceless_norm_sq",
    "sector_variance",
})

_SU_SECTOR_KEYS = frozenset({
    "hilbert_dim",
    "irrep_dim",
    "H_traceless_norm_sq",
    "casimir_ad_norm_sq",
    "rho_traceless_norm_sq",
    "M_traceless_norm_sq",
    "sector_variance",
})

_PHYS_KEYS = frozenset({
    "total_rho_proj_norm_sq",
    "rho_purity_bound_ok",
    "all_sectors_distinct",
})

_VALID_FAMILIES = frozenset({"so2n_majorana", "u1_hamming", "su_unstructured"})
_VALID_VERDICTS = frozenset({"CERTIFIED_TRAINABLE", "CERTIFIED_BARREN_PLATEAU", "CERTIFIED_BLIND_SPOT"})


def _assert_strict_json_types(obj: Any, path: str = "$") -> None:
    if obj is None:
        return
    if type(obj) is bool:
        raise ValueError(f"Forbidden bool at {path}: {obj!r}")
    if type(obj) is int:
        return
    if type(obj) is str:
        return
    if type(obj) is list:
        for idx, item in enumerate(obj):
            _assert_strict_json_types(item, f"{path}[{idx}]")
        return
    if type(obj) is dict:
        for k, v in obj.items():
            if type(k) is not str:
                raise ValueError(f"Non-string dict key at {path}: {k!r}")
            _assert_strict_json_types(v, f"{path}.{k}")
        return
    raise ValueError(f"Forbidden type {type(obj).__name__} at {path}: {obj!r}")


def canonical_fraction_str(x: Fraction) -> str:
    if not isinstance(x, Fraction):
        raise TypeError(f"Expected Fraction, got {type(x).__name__}")
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def parse_canonical_fraction(s: Any, field_name: str = "rational") -> Fraction:
    if type(s) is not str:
        raise ValueError(f"{field_name} must be a canonical str, got {type(s).__name__}")
    if _INT_RE.match(s):
        return Fraction(int(s), 1)
    m = _FRAC_RE.match(s)
    if not m:
        raise ValueError(f"Non-canonical rational string for {field_name}: {s!r}")
    num = int(m.group("num"))
    den = int(m.group("den"))
    if den < 2:
        raise ValueError(f"Denominator must be >= 2 for {field_name}: {s!r}")
    if math.gcd(abs(num), den) != 1:
        raise ValueError(f"Unreduced rational string for {field_name}: {s!r}")
    return Fraction(num, den)


def expected_dla_dim(family: str, n: int) -> int:
    if family == "so2n_majorana":
        return n * (2 * n - 1)
    if family == "u1_hamming":
        # Direct sum of su(C(n, w)) over w=0..n (plus n-dim abelian center of block phases, tracked as semisimple sum)
        return sum(math.comb(n, w) ** 2 - 1 for w in range(n + 1) if math.comb(n, w) >= 2)
    if family == "su_unstructured":
        return (1 << n) ** 2 - 1
    raise ValueError(f"Unknown family: {family!r}")


def build_lie_trainability_certificate(
    cert_id: str,
    family: str,
    n_qubits: int,
    threshold_tau: Fraction,
    sector_specs: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    """Build a self-contained zero-float LIE-TRAIN-CERT-1.0 certificate."""
    if type(cert_id) is not str or not _CERT_ID_RE.match(cert_id):
        raise ValueError(f"cert_id must match [A-Za-z0-9._-]{{1,120}}, got {cert_id!r}")
    if family not in _VALID_FAMILIES:
        raise ValueError(f"Unsupported family: {family!r}")
    if type(n_qubits) is not int or n_qubits < 2:
        raise ValueError(f"n_qubits must be int >= 2, got {n_qubits!r}")
    if threshold_tau <= 0:
        raise ValueError("threshold_tau must be strictly positive")
    if not sector_specs:
        raise ValueError("sector_specs must be non-empty")

    sectors_out: List[Dict[str, Any]] = []
    total_var = Fraction(0, 1)
    total_rho_proj = Fraction(0, 1)
    u1_weighted_purity = Fraction(0, 1)
    so2n_degrees_seen: set = set()
    so2n_gate_energy_ref: Fraction | None = None

    for sp in sector_specs:
        if family == "so2n_majorana":
            k = sp["k_degree"]
            chir = sp.get("chirality", "full")
            if chir == "full" and (n_qubits % 2 == 0) and k == n_qubits:
                raise ValueError(
                    f"Even n={n_qubits}, k={k} is reducible (Lambda^n_+ + Lambda^n_-); "
                    "explicit chiral sectors ('plus'/'minus') are required"
                )
            if (2 * n_qubits - k) in so2n_degrees_seen and (2 * n_qubits - k) != k:
                raise ValueError(
                    f"Complementary SO(2n) sectors k={k} and 2n-k={2*n_qubits-k} are isomorphic via P_ferm "
                    "(non-multiplicity-free)"
                )
            so2n_degrees_seen.add(k)
            ge = sp["gate_energy"]
            if so2n_gate_energy_ref is None:
                so2n_gate_energy_ref = ge
            elif ge != so2n_gate_energy_ref:
                raise ValueError(f"Inconsistent gate_energy across SO(2n) sectors: {ge} != {so2n_gate_energy_ref}")
            rn = sp["rho_proj_norm_sq"]
            mn = sp["M_proj_norm_sq"]
            dim_v, cas, var = so2n_sector_variance_exact(n_qubits, k, ge, rn, mn, chir)
            sectors_out.append({
                "k_degree": k,
                "chirality": chir,
                "irrep_dim": dim_v,
                "gate_energy": canonical_fraction_str(ge),
                "casimir_ad_norm_sq": canonical_fraction_str(cas),
                "rho_proj_norm_sq": canonical_fraction_str(rn),
                "M_proj_norm_sq": canonical_fraction_str(mn),
                "sector_variance": canonical_fraction_str(var),
            })
            total_rho_proj += rn
            total_var += var
        elif family == "u1_hamming":
            w = sp["hamming_weight"]
            hn = sp["H_traceless_norm_sq"]
            rn = sp["rho_traceless_norm_sq"]
            mn = sp["M_traceless_norm_sq"]
            dw = math.comb(n_qubits, w)
            dim_adj, cas, var = u1_sector_variance_exact(n_qubits, w, hn, rn, mn)
            if dw >= 2:
                max_sec_rho = Fraction(dw - 1, dw)
                if rn > max_sec_rho:
                    raise ValueError(f"Unphysical U(1) sector w={w} purity: {rn} > 1 - 1/{dw}")
                u1_weighted_purity += rn / max_sec_rho
            sectors_out.append({
                "hamming_weight": w,
                "block_dim": dw,
                "irrep_dim": dim_adj,
                "H_traceless_norm_sq": canonical_fraction_str(hn),
                "casimir_ad_norm_sq": canonical_fraction_str(cas),
                "rho_traceless_norm_sq": canonical_fraction_str(rn),
                "M_traceless_norm_sq": canonical_fraction_str(mn),
                "sector_variance": canonical_fraction_str(var),
            })
            total_rho_proj += rn
            total_var += var
        else:
            hn = sp["H_traceless_norm_sq"]
            rn = sp["rho_traceless_norm_sq"]
            mn = sp["M_traceless_norm_sq"]
            dim_adj, cas, var = su_unstructured_variance_exact(n_qubits, hn, rn, mn)
            sectors_out.append({
                "hilbert_dim": 1 << n_qubits,
                "irrep_dim": dim_adj,
                "H_traceless_norm_sq": canonical_fraction_str(hn),
                "casimir_ad_norm_sq": canonical_fraction_str(cas),
                "rho_traceless_norm_sq": canonical_fraction_str(rn),
                "M_traceless_norm_sq": canonical_fraction_str(mn),
                "sector_variance": canonical_fraction_str(var),
            })
            total_rho_proj += rn
            total_var += var

    max_traceless_purity = Fraction((1 << n_qubits) - 1, 1 << n_qubits)
    if total_rho_proj > max_traceless_purity:
        raise ValueError(
            f"Unphysical state purity: sum ||Pi(rho)||_HS^2 = {total_rho_proj} > 1 - 1/2^{n_qubits}"
        )
    if family == "u1_hamming" and u1_weighted_purity > Fraction(1, 1):
        raise ValueError(
            f"Unphysical U(1) multi-sector purity: sum ||rho_w^o||^2 / (1 - 1/d_w) = {u1_weighted_purity} > 1"
        )

    if total_var == 0:
        verdict = "CERTIFIED_BLIND_SPOT"
    elif total_var >= threshold_tau:
        verdict = "CERTIFIED_TRAINABLE"
    else:
        verdict = "CERTIFIED_BARREN_PLATEAU"

    cert: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "cert_id": cert_id,
        "family": family,
        "n_qubits": n_qubits,
        "dla_dimension": expected_dla_dim(family, n_qubits),
        "threshold_tau": canonical_fraction_str(threshold_tau),
        "verdict": verdict,
        "exact_variance": canonical_fraction_str(total_var),
        "sectors": sectors_out,
        "physical_checks": {
            "total_rho_proj_norm_sq": canonical_fraction_str(total_rho_proj),
            "rho_purity_bound_ok": 1,
            "all_sectors_distinct": 1,
        },
    }
    _assert_strict_json_types(cert)
    return cert


def verify_lie_trainability_certificate_exact(cert: Dict[str, Any]) -> Dict[str, Any]:
    """Independently verify a LIE-TRAIN-CERT-1.0 certificate with ZERO floating-point trust."""
    if type(cert) is not dict:
        raise ValueError("Certificate must be a dict")
    _assert_strict_json_types(cert)

    if cert.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"Unsupported schema_version: {cert.get('schema_version')!r}")
    if set(cert.keys()) != _TOP_KEYS:
        raise ValueError(
            f"Top-level key mismatch: missing={sorted(_TOP_KEYS - set(cert.keys()))}, "
            f"extra={sorted(set(cert.keys()) - _TOP_KEYS)}"
        )

    cert_id = cert["cert_id"]
    if type(cert_id) is not str or not _CERT_ID_RE.match(cert_id):
        raise ValueError(f"cert_id must match [A-Za-z0-9._-]{{1,120}}, got {cert_id!r}")

    family = cert["family"]
    if family not in _VALID_FAMILIES:
        raise ValueError(f"Unsupported family: {family!r}")

    n = cert["n_qubits"]
    if type(n) is not int or n < 2 or n > 500:
        raise ValueError(f"n_qubits must be strict int in [2, 500], got {n!r}")

    if cert["dla_dimension"] != expected_dla_dim(family, n):
        raise ValueError("dla_dimension mismatch")

    tau = parse_canonical_fraction(cert["threshold_tau"], "threshold_tau")
    if tau <= 0:
        raise ValueError("threshold_tau must be strictly positive")

    verdict = cert["verdict"]
    if verdict not in _VALID_VERDICTS:
        raise ValueError(f"Invalid verdict: {verdict!r}")

    sectors = cert["sectors"]
    if type(sectors) is not list or len(sectors) == 0:
        raise ValueError("sectors must be a non-empty list")

    seen_keys = set()
    seen_chiral_k: Dict[int, set] = {}
    seen_so2n_degrees: set = set()
    so2n_gate_energy_ref: Fraction | None = None
    total_var = Fraction(0, 1)
    total_rho = Fraction(0, 1)
    u1_weighted_purity = Fraction(0, 1)

    for idx, sec in enumerate(sectors):
        if type(sec) is not dict:
            raise ValueError(f"sectors[{idx}] must be a dict")
        if family == "so2n_majorana":
            if set(sec.keys()) != _SO2N_SECTOR_KEYS:
                raise ValueError(f"SO(2n) sector[{idx}] key mismatch")
            k = sec["k_degree"]
            if type(k) is not int or not (1 <= k <= 2 * n - 1):
                raise ValueError(f"Invalid k_degree at sectors[{idx}]: {k!r}")
            chir = sec["chirality"]
            if chir not in ("full", "plus", "minus"):
                raise ValueError(f"Invalid chirality at sectors[{idx}]: {chir!r}")
            # H02 fix: for even n, Lambda^n = Lambda^n_+ + Lambda^n_- is reducible; forbid 'full'!
            if chir == "full" and (n % 2 == 0) and k == n:
                raise ValueError(
                    f"Even n={n}, k={k} is reducible under SO(2n); 'chirality=full' is forbidden "
                    "(must decompose explicitly into 'plus'/'minus' real irreps)"
                )
            # H01 fix: Lambda^k ~= Lambda^{2n-k} via P_ferm under SO(2n); forbid complementary pairs!
            if (2 * n - k) in seen_so2n_degrees and (2 * n - k) != k:
                raise ValueError(
                    f"Complementary SO(2n) sectors k={k} and 2n-k={2*n-k} are isomorphic via P_ferm "
                    "(violates multiplicity-free condition)"
                )
            seen_so2n_degrees.add(k)
            key = (k, chir)
            if key in seen_keys:
                raise ValueError(f"Duplicate SO(2n) sector {key}")
            seen_keys.add(key)
            seen_chiral_k.setdefault(k, set()).add(chir)
            if "full" in seen_chiral_k[k] and len(seen_chiral_k[k]) > 1:
                raise ValueError(f"Cannot mix 'full' and chiral sectors for same degree k={k}")
            ge = parse_canonical_fraction(sec["gate_energy"], f"sectors[{idx}].gate_energy")
            # H17/H24 fix: a single circuit has ONE gate H and ONE scalar energy E_H across all sectors!
            if so2n_gate_energy_ref is None:
                so2n_gate_energy_ref = ge
            elif ge != so2n_gate_energy_ref:
                raise ValueError(
                    f"Inconsistent gate_energy at sectors[{idx}]: {ge} != {so2n_gate_energy_ref}"
                )
            rn = parse_canonical_fraction(sec["rho_proj_norm_sq"], f"sectors[{idx}].rho_proj_norm_sq")
            mn = parse_canonical_fraction(sec["M_proj_norm_sq"], f"sectors[{idx}].M_proj_norm_sq")
            dim_act, cas_act, var_act = so2n_sector_variance_exact(n, k, ge, rn, mn, chir)
            if sec["irrep_dim"] != dim_act:
                raise ValueError(f"irrep_dim mismatch at sectors[{idx}]")
            if parse_canonical_fraction(sec["casimir_ad_norm_sq"], f"sectors[{idx}].casimir") != cas_act:
                raise ValueError(f"casimir_ad_norm_sq mismatch at sectors[{idx}]")
            if parse_canonical_fraction(sec["sector_variance"], f"sectors[{idx}].sector_variance") != var_act:
                raise ValueError(f"sector_variance mismatch at sectors[{idx}]")
            total_rho += rn
            total_var += var_act

        elif family == "u1_hamming":
            if set(sec.keys()) != _U1_SECTOR_KEYS:
                raise ValueError(f"U(1) sector[{idx}] key mismatch")
            w = sec["hamming_weight"]
            if type(w) is not int or not (0 <= w <= n):
                raise ValueError(f"Invalid hamming_weight at sectors[{idx}]: {w!r}")
            if w in seen_keys:
                raise ValueError(f"Duplicate U(1) hamming_weight {w}")
            seen_keys.add(w)
            dw_act = math.comb(n, w)
            if sec["block_dim"] != dw_act:
                raise ValueError(f"block_dim mismatch at sectors[{idx}]")
            hn = parse_canonical_fraction(sec["H_traceless_norm_sq"], f"sectors[{idx}].H_traceless_norm_sq")
            rn = parse_canonical_fraction(sec["rho_traceless_norm_sq"], f"sectors[{idx}].rho_traceless_norm_sq")
            mn = parse_canonical_fraction(sec["M_traceless_norm_sq"], f"sectors[{idx}].M_traceless_norm_sq")
            dim_act, cas_act, var_act = u1_sector_variance_exact(n, w, hn, rn, mn)
            # H03/H04 fix: sharp per-sector and multi-sector U(1) purity bounds!
            if dw_act >= 2:
                max_sec_rho = Fraction(dw_act - 1, dw_act)
                if rn > max_sec_rho:
                    raise ValueError(
                        f"Unphysical U(1) sector w={w} purity: {rn} > 1 - 1/{dw_act}"
                    )
                u1_weighted_purity += rn / max_sec_rho
            if sec["irrep_dim"] != dim_act:
                raise ValueError(f"irrep_dim mismatch at sectors[{idx}]")
            if parse_canonical_fraction(sec["casimir_ad_norm_sq"], f"sectors[{idx}].casimir") != cas_act:
                raise ValueError(f"casimir_ad_norm_sq mismatch at sectors[{idx}]")
            if parse_canonical_fraction(sec["sector_variance"], f"sectors[{idx}].sector_variance") != var_act:
                raise ValueError(f"sector_variance mismatch at sectors[{idx}]")
            total_rho += rn
            total_var += var_act

        else:  # su_unstructured
            if len(sectors) != 1:
                raise ValueError("su_unstructured must have exactly 1 adjoint sector")
            if set(sec.keys()) != _SU_SECTOR_KEYS:
                raise ValueError(f"SU(2^n) sector[{idx}] key mismatch")
            if sec["hilbert_dim"] != (1 << n):
                raise ValueError("hilbert_dim mismatch")
            hn = parse_canonical_fraction(sec["H_traceless_norm_sq"], "H_traceless_norm_sq")
            rn = parse_canonical_fraction(sec["rho_traceless_norm_sq"], "rho_traceless_norm_sq")
            mn = parse_canonical_fraction(sec["M_traceless_norm_sq"], "M_traceless_norm_sq")
            dim_act, cas_act, var_act = su_unstructured_variance_exact(n, hn, rn, mn)
            if sec["irrep_dim"] != dim_act:
                raise ValueError("irrep_dim mismatch")
            if parse_canonical_fraction(sec["casimir_ad_norm_sq"], "casimir") != cas_act:
                raise ValueError("casimir_ad_norm_sq mismatch")
            if parse_canonical_fraction(sec["sector_variance"], "sector_variance") != var_act:
                raise ValueError("sector_variance mismatch")
            total_rho += rn
            total_var += var_act

    if parse_canonical_fraction(cert["exact_variance"], "exact_variance") != total_var:
        raise ValueError("exact_variance sum mismatch")

    # Sharp physical purity bound: sum ||Pi_lambda(rho)||_HS^2 <= Tr(rho^2) - 1/2^n <= 1 - 1/2^n
    max_traceless_purity = Fraction((1 << n) - 1, 1 << n)
    if total_rho > max_traceless_purity:
        raise ValueError(
            f"Unphysical state purity: total_rho_proj_norm_sq = {total_rho} > 1 - 1/2^{n}"
        )
    if family == "u1_hamming" and u1_weighted_purity > Fraction(1, 1):
        raise ValueError(
            f"Unphysical U(1) multi-sector purity: sum ||rho_w^o||^2 / (1 - 1/d_w) = {u1_weighted_purity} > 1"
        )

    pc = cert["physical_checks"]
    if type(pc) is not dict or set(pc.keys()) != _PHYS_KEYS:
        raise ValueError("physical_checks schema mismatch")
    if parse_canonical_fraction(pc["total_rho_proj_norm_sq"], "total_rho_proj_norm_sq") != total_rho:
        raise ValueError("total_rho_proj_norm_sq mismatch")
    if pc["rho_purity_bound_ok"] != 1 or pc["all_sectors_distinct"] != 1:
        raise ValueError("physical_checks flags must be strict int 1")

    # Enforce strict verdict mutual exclusivity (Lean 4 Theorems 15-16)
    if verdict == "CERTIFIED_BLIND_SPOT":
        if total_var != 0:
            raise ValueError(f"CERTIFIED_BLIND_SPOT requires exact_variance == 0, got {total_var}")
    elif verdict == "CERTIFIED_TRAINABLE":
        if total_var < tau:
            raise ValueError(f"CERTIFIED_TRAINABLE requires exact_variance >= tau ({total_var} < {tau})")
    elif verdict == "CERTIFIED_BARREN_PLATEAU":
        if not (0 < total_var < tau):
            raise ValueError(f"CERTIFIED_BARREN_PLATEAU requires 0 < exact_variance < tau (got {total_var}, tau={tau})")

    return {
        "verified": True,
        "cert_id": cert_id,
        "family": family,
        "n_qubits": n,
        "verdict": verdict,
        "exact_variance": str(total_var),
    }
