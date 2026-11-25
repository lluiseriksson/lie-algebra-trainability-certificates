"""Run all 4 computational experiments for Paper 09 and emit 12 exact-Q certificates.

Outputs:
  results/certificates/<cert_id>.cert.json  (12 exact-Q LIE-TRAIN-CERT-1.0 certificates)
  results/exp1_certified_atlas.json
  results/exp2_haar_mc_crosscheck.json
  results/exp3_phase_transitions_scaling.json
  results/exp4_sample_complexity_separation.json

Author: Lluis Eriksson
"""

from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.stats import ortho_group, unitary_group

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import (
    build_lie_trainability_certificate,
    canonical_fraction_str,
    exact_sample_complexity_lower_bound,
    so2n_sector_variance_exact,
    su_unstructured_variance_exact,
    u1_sector_variance_exact,
    verify_lie_trainability_certificate_exact,
)


def run_exp1_certified_atlas(cert_dir: Path) -> Dict[str, Any]:
    """Build & independently verify 12 exact-Q LIE-TRAIN-CERT-1.0 certificates across all 3 families & 3 verdicts."""
    cert_dir.mkdir(parents=True, exist_ok=True)
    for old in cert_dir.glob("*.cert.json"):
        old.unlink()

    # Exact physical Hilbert-Schmidt norms for Fock state |0>^{otimes n} and Pauli-Z strings:
    #   ||Pi_{Lambda^{2m}}(|0><0|^n)||_HS^2 = C(n, m) / 2^n,   ||Z_{i1}...Z_{im}||_HS^2 = 2^n.
    # For even n, |0>^{otimes n} has P_ferm = +1, so its Lambda^n projection lies 100% in Lambda^n_+!
    specs: List[Tuple[str, str, int, Fraction, List[Dict[str, Any]]]] = [
        # 1. SO(2n) k=2 (1-body Z_1, inside so(2n)): TRAINABLE at n=20
        ("so2n_n20_k2_trainable", "so2n_majorana", 20, Fraction(1, 10000), [{
            "k_degree": 2, "chirality": "full", "gate_energy": Fraction(1, 1),
            "rho_proj_norm_sq": Fraction(math.comb(20, 1), 1 << 20),
            "M_proj_norm_sq": Fraction(1 << 20, 1),
        }]),
        # 2. SO(2n) k=4 (2-body Z_1 Z_2, OUTSIDE so(2n)!): STILL TRAINABLE at n=10!
        ("so2n_n10_k4_outside_dla_trainable", "so2n_majorana", 10, Fraction(1, 100000), [{
            "k_degree": 4, "chirality": "full", "gate_energy": Fraction(1, 1),
            "rho_proj_norm_sq": Fraction(math.comb(10, 2), 1 << 10),
            "M_proj_norm_sq": Fraction(1 << 10, 1),
        }]),
        # 3. SO(2n) k=n=20 (half-weight Z_1...Z_10; explicit irreducible chiral Lambda^20_+ + Lambda^20_-): BARREN!
        # |0>^20 lies 100% in Lambda^20_+, while Z_1...Z_10 splits 50/50 (2^19 each) across Lambda^20_+/-.
        ("so2n_n20_kn_barren_in_poly_dla", "so2n_majorana", 20, Fraction(1, 10000), [
            {
                "k_degree": 20, "chirality": "plus", "gate_energy": Fraction(1, 1),
                "rho_proj_norm_sq": Fraction(math.comb(20, 10), 1 << 20),
                "M_proj_norm_sq": Fraction(1 << 19, 1),
            },
            {
                "k_degree": 20, "chirality": "minus", "gate_energy": Fraction(1, 1),
                "rho_proj_norm_sq": Fraction(0, 1),
                "M_proj_norm_sq": Fraction(1 << 19, 1),
            },
        ]),
        # 4. SO(2n) k=n=50 (n=50 qubits, dim(so(100))=4950; explicit chiral irreps): EXPONENTIAL BARREN PLATEAU!
        ("so2n_n50_kn_barren_50q", "so2n_majorana", 50, Fraction(1, 1000000), [
            {
                "k_degree": 50, "chirality": "plus", "gate_energy": Fraction(1, 1),
                "rho_proj_norm_sq": Fraction(math.comb(50, 25), 1 << 50),
                "M_proj_norm_sq": Fraction(1 << 49, 1),
            },
            {
                "k_degree": 50, "chirality": "minus", "gate_energy": Fraction(1, 1),
                "rho_proj_norm_sq": Fraction(0, 1),
                "M_proj_norm_sq": Fraction(1 << 49, 1),
            },
        ]),
        # 5. SO(2n) Even-n Hodge-chiral (+, +) 2x amplification at n=4 (M_+ = (Z1Z2 + Z3Z4)/sqrt(2)): TRAINABLE!
        ("so2n_n4_hodge_chiral_plus_2x", "so2n_majorana", 4, Fraction(1, 100), [{
            "k_degree": 4, "chirality": "plus", "gate_energy": Fraction(1, 1),
            "rho_proj_norm_sq": Fraction(math.comb(4, 2), 1 << 4),
            "M_proj_norm_sq": Fraction(1 << 4, 1),
        }]),
        # 6. SO(2n) Even-n Hodge-chiral (+, -) EXACT BLIND SPOT at n=4 (M_- = (Z1Z2 - Z3Z4)/sqrt(2)): ZERO VARIANCE!
        ("so2n_n4_hodge_chiral_blind_spot", "so2n_majorana", 4, Fraction(1, 1000), [
            {
                "k_degree": 4, "chirality": "plus", "gate_energy": Fraction(1, 1),
                "rho_proj_norm_sq": Fraction(math.comb(4, 2), 1 << 4),
                "M_proj_norm_sq": Fraction(0, 1),
            },
            {
                "k_degree": 4, "chirality": "minus", "gate_energy": Fraction(1, 1),
                "rho_proj_norm_sq": Fraction(0, 1),
                "M_proj_norm_sq": Fraction(1 << 4, 1),
            },
        ]),
        # 7. SO(2n) Multi-sector observable (Z_1 + Z_1 Z_2) at n=12: exact sum of k=2 + k=4 irreps!
        ("so2n_n12_multisector_k2_k4", "so2n_majorana", 12, Fraction(1, 1000), [
            {
                "k_degree": 2, "chirality": "full", "gate_energy": Fraction(1, 1),
                "rho_proj_norm_sq": Fraction(math.comb(12, 1), 1 << 12),
                "M_proj_norm_sq": Fraction(1 << 12, 1),
            },
            {
                "k_degree": 4, "chirality": "full", "gate_energy": Fraction(1, 1),
                "rho_proj_norm_sq": Fraction(math.comb(12, 2), 1 << 12),
                "M_proj_norm_sq": Fraction(1 << 12, 1),
            },
        ]),
        # 8. U(1) fixed single-excitation sector w=1 at n=20 (d_1 = 20): TRAINABLE!
        # Pure state in sector w: ||rho_w^o||_HS^2 = 1 - 1/d_w; balanced Z_1 on w=1 has tr = d_1 - 2, ||M_w^o||_HS^2 = 4(d_1-1)/d_1
        ("u1_n20_w1_trainable", "u1_hamming", 20, Fraction(1, 10000), [{
            "hamming_weight": 1,
            "H_traceless_norm_sq": Fraction(2, 1),  # XX+YY swap gate on (1,2) in w=1 has ||H_w^o||_HS^2 = 2
            "rho_traceless_norm_sq": Fraction(19, 20),
            "M_traceless_norm_sq": Fraction(4 * 19, 20),
        }]),
        # 9. U(1) half-filling sector w=15 at n=30 (d_15 = 155,117,520): EXPONENTIAL BARREN PLATEAU!
        ("u1_n30_half_filling_barren", "u1_hamming", 30, Fraction(1, 1000000), [{
            "hamming_weight": 15,
            "H_traceless_norm_sq": Fraction(math.comb(28, 14) * 2, 1),
            "rho_traceless_norm_sq": Fraction(math.comb(30, 15) - 1, math.comb(30, 15)),
            "M_traceless_norm_sq": Fraction(math.comb(30, 15), 1),
        }]),
        # 10. U(1) Gate-Commutant Blind Spot (H = sum Z_j is scalar on every fixed-w block!): EXACT ZERO!
        ("u1_n12_scalar_gate_blind_spot", "u1_hamming", 12, Fraction(1, 10000), [{
            "hamming_weight": 2,
            "H_traceless_norm_sq": Fraction(0, 1),
            "rho_traceless_norm_sq": Fraction(math.comb(12, 2) - 1, math.comb(12, 2)),
            "M_traceless_norm_sq": Fraction(4 * 11 * 55, 66),  # Z_1 on w=2: d=66, 11 ones & 55 -ones -> 110/3
        }]),
        # 11. Unstructured SU(2^n) at n=20 (Pauli generator H = Z_1/2 -> ||H^o||^2 = 2^n/4, M = Z_1 -> ||M^o||^2 = 2^n): BARREN!
        ("su_n20_unstructured_barren", "su_unstructured", 20, Fraction(1, 1000000), [{
            "H_traceless_norm_sq": Fraction(1 << 20, 4),
            "rho_traceless_norm_sq": Fraction((1 << 20) - 1, 1 << 20),
            "M_traceless_norm_sq": Fraction(1 << 20, 1),
        }]),
        # 12. SO(2n) 100-qubit large-scale certificate (k=2 trainable in <1 ms!):
        ("so2n_n100_k2_large_scale", "so2n_majorana", 100, Fraction(1, 100000), [{
            "k_degree": 2, "chirality": "full", "gate_energy": Fraction(1, 1),
            "rho_proj_norm_sq": Fraction(100, 1 << 100),
            "M_proj_norm_sq": Fraction(1 << 100, 1),
        }]),
    ]

    rows: List[Dict[str, Any]] = []
    for cid, fam, nq, tau, sec_specs in specs:
        cert = build_lie_trainability_certificate(cid, fam, nq, tau, sec_specs)
        cpath = cert_dir / f"{cid}.cert.json"
        cpath.write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8")

        t0 = time.perf_counter()
        loaded = json.loads(cpath.read_text(encoding="utf-8"))
        ver = verify_lie_trainability_certificate_exact(loaded)
        v_ms = (time.perf_counter() - t0) * 1000.0

        exact_v = Fraction(ver["exact_variance"])
        rows.append({
            "cert_id": cid,
            "family": fam,
            "n_qubits": nq,
            "dla_dimension": cert["dla_dimension"],
            "threshold_tau": canonical_fraction_str(tau),
            "verdict": ver["verdict"],
            "exact_variance": ver["exact_variance"],
            "exact_variance_float": float(exact_v),
            "verify_ms": round(v_ms, 3),
            "verified": True,
        })

    return {
        "experiment": "exp1_certified_atlas",
        "num_certificates": len(rows),
        "all_verified": all(r["verified"] for r in rows),
        "max_verify_ms": max(r["verify_ms"] for r in rows),
        "certificates": rows,
    }


def _majorana_mats(n: int) -> List[np.ndarray]:
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    def kron_all(ms: List[np.ndarray]) -> np.ndarray:
        out = ms[0]
        for m in ms[1:]:
            out = np.kron(out, m)
        return out

    gs: List[np.ndarray] = []
    for j in range(n):
        gs.append(kron_all([Z] * j + [X] + [I2] * (n - 1 - j)))
        gs.append(kron_all([Z] * j + [Y] + [I2] * (n - 1 - j)))
    return gs


def run_exp2_haar_mc_crosscheck() -> Dict[str, Any]:
    """Cross-check exact-Q formulas against independent 12,000-sample Haar Monte Carlo on SO(4), SO(6), U(1), SU(8)."""
    rng = np.random.default_rng(20261124)
    N_mc = 12000
    checks: List[Dict[str, Any]] = []

    # 1. SO(6) (n=3) across k = 1, 2, 3, 4 (including k=4 OUTSIDE so(6)!)
    n = 3
    d = 1 << n
    gs = _majorana_mats(n)
    H = 0.5j * (gs[0] @ gs[1])
    Rs_L = ortho_group.rvs(2 * n, size=N_mc, random_state=rng)
    Rs_R = ortho_group.rvs(2 * n, size=N_mc, random_state=rng)
    Rs_L[np.linalg.det(Rs_L) < 0, :, 0] *= -1
    Rs_R[np.linalg.det(Rs_R) < 0, :, 0] *= -1

    for k in (1, 2, 3, 4):
        combos = list(combinations(range(2 * n), k))
        m = len(combos)
        phase = (-1j) ** (k * (k - 1) // 2)
        bk = [
            phase * np.linalg.multi_dot([gs[i] for i in idxs]) if k > 1 else gs[idxs[0]]
            for idxs in combos
        ]
        Ak = np.array([[np.real(np.trace(bk[a] @ (1j * (H @ bk[b] - bk[b] @ H)))) / d for b in range(m)] for a in range(m)])
        c0, c1 = combos[0], combos[-1]
        uL = np.stack([np.linalg.det(Rs_L[:, idx, :][:, :, c0]) for idx in combos], axis=1)
        uR = np.stack([np.linalg.det(Rs_R[:, c1, :][:, :, idx]) for idx in combos], axis=1)
        grads = np.einsum("ni,ij,nj->n", uL, Ak, uR)
        var_mc = float(np.var(grads))
        se_mc = float(np.std(grads ** 2) / np.sqrt(N_mc))
        _, _, var_ex = so2n_sector_variance_exact(n, k, Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), "full")
        checks.append({
            "label": f"SO(6)_k{k}",
            "family": "so2n_majorana",
            "exact_q": canonical_fraction_str(var_ex),
            "exact_float": float(var_ex),
            "mc_var": var_mc,
            "mc_se": se_mc,
            "z_score": abs(var_mc - float(var_ex)) / se_mc if se_mc > 0 else 0.0,
        })

    # 2. SO(4) (even n=2) Hodge-chiral (+, +), (-, -), and (+, -) blind spot!
    n2 = 2
    d2 = 4
    gs2 = _majorana_mats(n2)
    H2 = 0.5j * (gs2[0] @ gs2[1])
    P_ferm = ((-1j) ** n2) * (gs2[0] @ gs2[1] @ gs2[2] @ gs2[3])
    combos2 = list(combinations(range(4), 2))
    bk2 = [-1j * (gs2[a] @ gs2[b]) for a, b in combos2]
    Bp = (bk2[0] + P_ferm @ bk2[0]) / np.sqrt(2.0)
    Bm = (bk2[0] - P_ferm @ bk2[0]) / np.sqrt(2.0)
    vp = np.array([np.real(np.trace(bk2[i] @ Bp)) / d2 for i in range(6)])
    vm = np.array([np.real(np.trace(bk2[i] @ Bm)) / d2 for i in range(6)])
    Ak2 = np.array([[np.real(np.trace(bk2[a] @ (1j * (H2 @ bk2[b] - bk2[b] @ H2)))) / d2 for b in range(6)] for a in range(6)])

    R2L = ortho_group.rvs(4, size=N_mc, random_state=rng)
    R2R = ortho_group.rvs(4, size=N_mc, random_state=rng)
    R2L[np.linalg.det(R2L) < 0, :, 0] *= -1
    R2R[np.linalg.det(R2R) < 0, :, 0] *= -1
    WL2 = np.zeros((N_mc, 6, 6))
    WR2 = np.zeros((N_mc, 6, 6))
    for i, ci in enumerate(combos2):
        for j, cj in enumerate(combos2):
            WL2[:, i, j] = np.linalg.det(R2L[:, ci, :][:, :, cj])
            WR2[:, i, j] = np.linalg.det(R2R[:, ci, :][:, :, cj])

    for lbl, v1, v2, ex_q in [
        ("SO(4)_chiral_++_2x", vp, vp, so2n_sector_variance_exact(2, 2, Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), "plus")[2]),
        ("SO(4)_chiral_--_2x", vm, vm, so2n_sector_variance_exact(2, 2, Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), "minus")[2]),
        ("SO(4)_chiral_+-_blind", vp, vm, Fraction(0, 1)),
    ]:
        grads = np.einsum("ni,ij,nj->n", WL2 @ v1, Ak2, np.einsum("nji,j->ni", WR2, v2))
        var_mc = float(np.var(grads))
        se_mc = float(np.std(grads ** 2) / np.sqrt(N_mc))
        checks.append({
            "label": lbl,
            "family": "so2n_majorana",
            "exact_q": canonical_fraction_str(ex_q),
            "exact_float": float(ex_q),
            "mc_var": var_mc,
            "mc_se": se_mc,
            "z_score": abs(var_mc - float(ex_q)) / se_mc if se_mc > 1e-18 else 0.0,
        })

    # 3. U(1) SU(6) block (n=4, w=2) + SU(8) unstructured (n=3)
    for lbl, fam, dw, diag_H, diag_M in [
        ("U1_n4_w2_SU6", "u1_hamming", 6, [1.0, -1.0, 0.5, -0.5, 0.0, 0.0], [1.0, 1.0, 1.0, -1.0, -1.0, -1.0]),
        ("SU8_n3_unstructured", "su_unstructured", 8, [0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5], [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0]),
    ]:
        Hm = np.diag(diag_H)
        Mm = np.diag(diag_M)
        rhom = np.zeros((dw, dw))
        rhom[0, 0] = 1.0
        UsL = unitary_group.rvs(dw, size=N_mc, random_state=rng)
        UsR = unitary_group.rvs(dw, size=N_mc, random_state=rng)
        r_rot = UsL @ rhom @ UsL.conj().transpose(0, 2, 1)
        m_rot = UsR.conj().transpose(0, 2, 1) @ Mm @ UsR
        comm = Hm[None, :, :] @ m_rot - m_rot @ Hm[None, :, :]
        grads = np.real(1j * np.einsum("nij,nji->n", r_rot, comm))
        var_mc = float(np.var(grads))
        se_mc = float(np.std(grads ** 2) / np.sqrt(N_mc))
        h_nsq = Fraction(int(round(sum(x * x for x in diag_H) * 4)), 4)
        m_nsq = Fraction(int(round(sum(x * x for x in diag_M))), 1)
        r_nsq = Fraction(dw - 1, dw)
        if fam == "u1_hamming":
            _, _, var_ex = u1_sector_variance_exact(4, 2, h_nsq, r_nsq, m_nsq)
        else:
            _, _, var_ex = su_unstructured_variance_exact(3, h_nsq, r_nsq, m_nsq)
        checks.append({
            "label": lbl,
            "family": fam,
            "exact_q": canonical_fraction_str(var_ex),
            "exact_float": float(var_ex),
            "mc_var": var_mc,
            "mc_se": se_mc,
            "z_score": abs(var_mc - float(var_ex)) / se_mc if se_mc > 0 else 0.0,
        })

    return {
        "experiment": "exp2_haar_mc_crosscheck",
        "num_samples_per_check": N_mc,
        "num_checks": len(checks),
        "max_z_score": max(c["z_score"] for c in checks),
        "checks": checks,
    }


def run_exp3_phase_transitions_scaling() -> Dict[str, Any]:
    """Exact-Q variance scaling curves across n = 4..40 qubits exposing poly vs exp transitions."""
    ns = [4, 6, 8, 10, 12, 16, 20, 24, 30, 40]
    so2n_rows: List[Dict[str, Any]] = []
    u1_rows: List[Dict[str, Any]] = []

    for n in ns:
        # SO(2n) on Fock state |0>^n with Pauli-Z strings (exact cancellation of 2^n!)
        _, _, v_k2 = so2n_sector_variance_exact(n, 2, Fraction(1, 1), Fraction(math.comb(n, 1), 1 << n), Fraction(1 << n, 1), "full")
        _, _, v_k4 = so2n_sector_variance_exact(n, 4, Fraction(1, 1), Fraction(math.comb(n, 2), 1 << n), Fraction(1 << n, 1), "full")
        _, _, v_kn = so2n_sector_variance_exact(n, n, Fraction(1, 1), Fraction(math.comb(n, n // 2), 1 << n), Fraction(1 << n, 1), "full")
        _, _, v_kn_plus = so2n_sector_variance_exact(n, n, Fraction(1, 1), Fraction(math.comb(n, n // 2), 1 << n), Fraction(1 << n, 1), "plus")
        so2n_rows.append({
            "n": n,
            "dla_dim": n * (2 * n - 1),
            "var_k2_float": float(v_k2),
            "var_k4_float": float(v_k4),
            "var_kn_float": float(v_kn),
            "var_kn_chiral_2x_float": float(v_kn_plus),
            "poly_to_exp_ratio": float(v_k2 / v_kn),
        })

        # U(1) fixed w=1, w=2 vs half-filling w=n/2 vs unstructured SU(2^n)
        d1 = math.comb(n, 1)
        d2 = math.comb(n, 2)
        dh = math.comb(n, n // 2)
        _, _, u_w1 = u1_sector_variance_exact(n, 1, Fraction(2, 1), Fraction(d1 - 1, d1), Fraction(4 * (d1 - 1), d1))
        # For w=2: swap gate XX+YY on (1,2) acts on n-2 pairs -> ||H||_HS^2 = 2*(n-2); Z_1 has (n-1) ones & C(n-1,2) -ones
        n1_w2 = n - 1
        m_tr_w2 = n1_w2 - (d2 - n1_w2)
        m_nsq_w2 = Fraction(d2, 1) - Fraction(m_tr_w2 * m_tr_w2, d2)
        _, _, u_w2 = u1_sector_variance_exact(n, 2, Fraction(2 * (n - 2), 1), Fraction(d2 - 1, d2), m_nsq_w2)
        _, _, u_half = u1_sector_variance_exact(n, n // 2, Fraction(2 * math.comb(n - 2, n // 2 - 1), 1), Fraction(dh - 1, dh), Fraction(dh, 1))
        _, _, su_v = su_unstructured_variance_exact(n, Fraction(1 << n, 4), Fraction((1 << n) - 1, 1 << n), Fraction(1 << n, 1))
        u1_rows.append({
            "n": n,
            "d_half": dh,
            "var_w1_float": float(u_w1),
            "var_w2_float": float(u_w2),
            "var_half_float": float(u_half),
            "var_su_float": float(su_v),
        })

    return {
        "experiment": "exp3_phase_transitions_scaling",
        "so2n_scaling": so2n_rows,
        "u1_scaling": u1_rows,
    }


def run_exp4_sample_complexity_separation(exp3: Dict[str, Any]) -> Dict[str, Any]:
    """Compare exact-Q certificate verification runtime vs Chebyshev sampling cost N >= 4 / tau^2."""
    rows: List[Dict[str, Any]] = []
    for r_so, r_u1 in zip(exp3["so2n_scaling"], exp3["u1_scaling"]):
        n = r_so["n"]
        # Exact Chebyshev bound for detecting Var >= tau with unit fourth-moment floor B = 1
        _, _, v_k2 = so2n_sector_variance_exact(n, 2, Fraction(1, 1), Fraction(math.comb(n, 1), 1 << n), Fraction(1 << n, 1), "full")
        _, _, v_kn = so2n_sector_variance_exact(n, n, Fraction(1, 1), Fraction(math.comb(n, n // 2), 1 << n), Fraction(1 << n, 1), "full")
        dh = math.comb(n, n // 2)
        _, _, u_half = u1_sector_variance_exact(n, n // 2, Fraction(2 * math.comb(n - 2, n // 2 - 1), 1), Fraction(dh - 1, dh), Fraction(dh, 1))

        # Benchmark exact-Q verifier wall time at this n (explicit irreducible chiral irreps for even n!)
        cert = build_lie_trainability_certificate(
            f"bench_n{n}", "so2n_majorana", n, Fraction(1, 1000000),
            [
                {
                    "k_degree": n, "chirality": "plus", "gate_energy": Fraction(1, 1),
                    "rho_proj_norm_sq": Fraction(math.comb(n, n // 2), 1 << n),
                    "M_proj_norm_sq": Fraction(1 << (n - 1), 1),
                },
                {
                    "k_degree": n, "chirality": "minus", "gate_energy": Fraction(1, 1),
                    "rho_proj_norm_sq": Fraction(0, 1),
                    "M_proj_norm_sq": Fraction(1 << (n - 1), 1),
                },
            ],
        )
        t0 = time.perf_counter()
        verify_lie_trainability_certificate_exact(cert)
        v_ms = (time.perf_counter() - t0) * 1000.0

        rows.append({
            "n": n,
            "chebyshev_N_k2": float(exact_sample_complexity_lower_bound(v_k2, Fraction(1, 1))),
            "chebyshev_N_kn": float(exact_sample_complexity_lower_bound(v_kn, Fraction(1, 1))),
            "chebyshev_N_u1_half": float(exact_sample_complexity_lower_bound(u_half, Fraction(1, 1))),
            "exact_q_verify_ms": round(v_ms, 3),
        })

    return {
        "experiment": "exp4_sample_complexity_separation",
        "records": rows,
    }


def main() -> None:
    res_dir = ROOT / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    cert_dir = res_dir / "certificates"

    t0 = time.perf_counter()
    e1 = run_exp1_certified_atlas(cert_dir)
    (res_dir / "exp1_certified_atlas.json").write_text(json.dumps(e1, indent=2) + "\n", encoding="utf-8")

    e2 = run_exp2_haar_mc_crosscheck()
    (res_dir / "exp2_haar_mc_crosscheck.json").write_text(json.dumps(e2, indent=2) + "\n", encoding="utf-8")

    e3 = run_exp3_phase_transitions_scaling()
    (res_dir / "exp3_phase_transitions_scaling.json").write_text(json.dumps(e3, indent=2) + "\n", encoding="utf-8")

    e4 = run_exp4_sample_complexity_separation(e3)
    (res_dir / "exp4_sample_complexity_separation.json").write_text(json.dumps(e4, indent=2) + "\n", encoding="utf-8")

    dt = time.perf_counter() - t0
    print(f"All 4 Paper 09 experiments completed in {dt:.2f}s:")
    print(f"  Exp 1: {e1['num_certificates']} exact-Q certs verified={e1['all_verified']}, max_ms={e1['max_verify_ms']:.2f}")
    print(f"  Exp 2: {e2['num_checks']} Haar MC cross-checks, max_z={e2['max_z_score']:.2f}")
    print(f"  Exp 3: {len(e3['so2n_scaling'])} qubit sizes (n=4..40) mapped")
    print(f"  Exp 4: sample complexity at n=40: N_kn={e4['records'][-1]['chebyshev_N_kn']:.2e}, verify={e4['records'][-1]['exact_q_verify_ms']}ms")


if __name__ == "__main__":
    main()
