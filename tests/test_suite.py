"""Comprehensive test & 48-attack adversarial red-team suite for Paper 09.

Author: Lluis Eriksson
"""

from __future__ import annotations

import copy
import json
import math
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import pytest

from src import (
    so2n_sector_variance_exact,
    su_unstructured_variance_exact,
    u1_sector_variance_exact,
    verify_lie_trainability_certificate_exact,
)

ROOT = Path(__file__).resolve().parent.parent
CERT_DIR = ROOT / "results" / "certificates"


def test_exact_algebra_identities_across_families() -> None:
    """Verify exact binomial Casimir identity, Hodge-chiral 2x/0x laws, and U(1)/SU formulas."""
    for n in (2, 3, 4, 6, 10):
        for k in range(1, 2 * n):
            dim_v, cas, var = so2n_sector_variance_exact(n, k, Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), "full")
            assert dim_v == math.comb(2 * n, k)
            assert cas == Fraction(2 * math.comb(2 * n - 2, k - 1), 1)
            # Exact Theorem 3.2 simplified rational ratio: k*(2n-k) / (n*(2n-1)*C(2n, k))
            assert var == Fraction(k * (2 * n - k), n * (2 * n - 1) * math.comb(2 * n, k))

        if n % 2 == 0:
            _, _, v_full = so2n_sector_variance_exact(n, n, Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), "full")
            _, _, v_plus = so2n_sector_variance_exact(n, n, Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), "plus")
            _, _, v_minus = so2n_sector_variance_exact(n, n, Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), "minus")
            assert v_plus == 2 * v_full
            assert v_minus == 2 * v_full

    # U(1) blind spot when H_traceless == 0
    _, cas0, var0 = u1_sector_variance_exact(8, 2, Fraction(0, 1), Fraction(27, 28), Fraction(28, 1))
    assert cas0 == 0 and var0 == 0

    # SU(2^n) exact McClean-style scaling
    dim_su, _, var_su = su_unstructured_variance_exact(4, Fraction(4, 1), Fraction(15, 16), Fraction(16, 1))
    assert dim_su == 255
    assert var_su == Fraction(2 * 16 * 4 * 15, 255 * 255)


def test_all_12_shipped_certificates_verify() -> None:
    certs = sorted(CERT_DIR.glob("*.cert.json"))
    assert len(certs) == 12
    families = set()
    verdicts = set()
    for cp in certs:
        data = json.loads(cp.read_text(encoding="utf-8"))
        res = verify_lie_trainability_certificate_exact(data)
        assert res["verified"] is True
        families.add(res["family"])
        verdicts.add(res["verdict"])
    assert families == {"so2n_majorana", "u1_hamming", "su_unstructured"}
    assert verdicts == {"CERTIFIED_TRAINABLE", "CERTIFIED_BARREN_PLATEAU", "CERTIFIED_BLIND_SPOT"}


def test_48_live_adversarial_red_team_attacks_rejected() -> None:
    """Every single adversarial corruption across all 3 families & 3 verdicts MUST be rejected."""
    c_so = json.loads((CERT_DIR / "so2n_n20_k2_trainable.cert.json").read_text(encoding="utf-8"))
    c_ch = json.loads((CERT_DIR / "so2n_n4_hodge_chiral_blind_spot.cert.json").read_text(encoding="utf-8"))
    c_ms = json.loads((CERT_DIR / "so2n_n12_multisector_k2_k4.cert.json").read_text(encoding="utf-8"))
    c_u1 = json.loads((CERT_DIR / "u1_n20_w1_trainable.cert.json").read_text(encoding="utf-8"))
    c_su = json.loads((CERT_DIR / "su_n20_unstructured_barren.cert.json").read_text(encoding="utf-8"))

    attacks = []

    def add(name, fn, base=c_so):
        cand = copy.deepcopy(base)
        fn(cand)
        attacks.append((name, cand))

    # 1-10: Schema, type, and key-set attacks
    add("01_bad_schema", lambda c: c.__setitem__("schema_version", "FAKE-1.0"))
    add("02_missing_key", lambda c: c.pop("exact_variance"))
    add("03_extra_key", lambda c: c.__setitem__("rogue", 1))
    add("04_blank_cert_id", lambda c: c.__setitem__("cert_id", "   "))
    add("05_bool_n_qubits", lambda c: c.__setitem__("n_qubits", True))
    add("06_float_variance", lambda c: c.__setitem__("exact_variance", 0.01))
    add("07_decimal_tau", lambda c: c.__setitem__("threshold_tau", Decimal("0.01")))
    add("08_tuple_sectors", lambda c: c.__setitem__("sectors", tuple(c["sectors"])))
    add("09_non_str_key", lambda c: c.__setitem__(99, "bad"))
    add("10_bad_family", lambda c: c.__setitem__("family", "e8_exceptional"))

    # 11-18: Canonical rational parser attacks
    add("11_unreduced_frac", lambda c: c.__setitem__("threshold_tau", "2/20000"))
    add("12_denom_one", lambda c: c["sectors"][0].__setitem__("gate_energy", "1/1"))
    add("13_neg_zero", lambda c: c["sectors"][0].__setitem__("M_proj_norm_sq", "-0"), base=c_ch)
    add("14_leading_plus", lambda c: c["sectors"][0].__setitem__("gate_energy", "+1"))
    add("15_leading_zero", lambda c: c["sectors"][0].__setitem__("gate_energy", "01"))
    add("16_trailing_newline", lambda c: c.__setitem__("threshold_tau", "1/10000\n"))
    add("17_unicode_digit", lambda c: c["sectors"][0].__setitem__("gate_energy", "\u0661"))
    add("18_zero_tau", lambda c: c.__setitem__("threshold_tau", "0"))

    # 19-28: SO(2n) algebra, Hodge-chiral, and physical purity attacks
    add("19_wrong_dla_dim", lambda c: c.__setitem__("dla_dimension", 999))
    add("20_wrong_irrep_dim", lambda c: c["sectors"][0].__setitem__("irrep_dim", 781))
    add("21_wrong_casimir", lambda c: c["sectors"][0].__setitem__("casimir_ad_norm_sq", "77"))
    add("22_wrong_sector_var", lambda c: c["sectors"][0].__setitem__("sector_variance", "1"))
    add("23_chiral_on_odd_n", lambda c: (c.__setitem__("n_qubits", 3), c.__setitem__("dla_dimension", 15), c["sectors"][0].__setitem__("k_degree", 3)), base=c_ch)
    add("24_chiral_on_k_ne_n", lambda c: c["sectors"][0].__setitem__("chirality", "plus"))
    add("25_duplicate_sector", lambda c: c["sectors"].append(copy.deepcopy(c["sectors"][0])))
    add("26_mix_full_and_chiral", lambda c: c["sectors"].__setitem__(1, {**c["sectors"][1], "chirality": "full"}), base=c_ch)
    add("27_unphysical_purity_gt_1", lambda c: (
        c["sectors"][0].__setitem__("rho_proj_norm_sq", "2"),
        c["physical_checks"].__setitem__("total_rho_proj_norm_sq", "2"),
    ))
    add("28_negative_gate_energy", lambda c: c["sectors"][0].__setitem__("gate_energy", "-1"))

    # 29-34: U(1) and SU(2^n) attacks
    add("29_u1_wrong_block_dim", lambda c: c["sectors"][0].__setitem__("block_dim", 21), base=c_u1)
    add("30_u1_wrong_irrep_dim", lambda c: c["sectors"][0].__setitem__("irrep_dim", 400), base=c_u1)
    add("31_u1_out_of_range_w", lambda c: c["sectors"][0].__setitem__("hamming_weight", 25), base=c_u1)
    add("32_su_wrong_hilbert_dim", lambda c: c["sectors"][0].__setitem__("hilbert_dim", 1024), base=c_su)
    add("33_su_two_sectors", lambda c: c["sectors"].append(copy.deepcopy(c["sectors"][0])), base=c_su)
    add("34_phys_flag_bool", lambda c: c["physical_checks"].__setitem__("rho_purity_bound_ok", True))

    # 35-40: Verdict mutual-exclusivity & sharp traceless purity attacks
    add("35_trainable_claimed_barren", lambda c: c.__setitem__("verdict", "CERTIFIED_BARREN_PLATEAU"))
    add("36_trainable_claimed_blind", lambda c: c.__setitem__("verdict", "CERTIFIED_BLIND_SPOT"))
    add("37_blind_claimed_trainable", lambda c: c.__setitem__("verdict", "CERTIFIED_TRAINABLE"), base=c_ch)
    add("38_barren_claimed_trainable", lambda c: c.__setitem__("verdict", "CERTIFIED_TRAINABLE"), base=c_su)
    add("39_sharp_traceless_purity_eq_1", lambda c: (
        c["sectors"][0].__setitem__("rho_traceless_norm_sq", "1"),
        c["physical_checks"].__setitem__("total_rho_proj_norm_sq", "1"),
    ), base=c_su)
    add("40_u1_duplicate_weight", lambda c: c["sectors"].append(copy.deepcopy(c["sectors"][0])), base=c_u1)

    # 41-48: Hostile representation-theoretic, physical-purity, consistency, and hygiene attacks
    # 41 (H01): Complementary SO(2n) sectors k and 2n-k (non-multiplicity-free via P_ferm!)
    attacks.append(("41_complementary_so2n_k_2nk", {
        "schema_version": "LIE-TRAIN-CERT-1.0",
        "cert_id": "forged_so6_comp",
        "family": "so2n_majorana",
        "n_qubits": 3,
        "dla_dimension": 15,
        "threshold_tau": "1/100",
        "verdict": "CERTIFIED_TRAINABLE",
        "exact_variance": "8/75",
        "sectors": [
            {
                "k_degree": 2, "chirality": "full", "irrep_dim": 15, "gate_energy": "1",
                "casimir_ad_norm_sq": "8", "rho_proj_norm_sq": "3/8", "M_proj_norm_sq": "4",
                "sector_variance": "4/75",
            },
            {
                "k_degree": 4, "chirality": "full", "irrep_dim": 15, "gate_energy": "1",
                "casimir_ad_norm_sq": "8", "rho_proj_norm_sq": "3/8", "M_proj_norm_sq": "4",
                "sector_variance": "4/75",
            },
        ],
        "physical_checks": {"total_rho_proj_norm_sq": "3/4", "rho_purity_bound_ok": 1, "all_sectors_distinct": 1},
    }))
    # 42 (H02): Even-n reducible k=n with chirality="full" (hides chiral blind spot!)
    attacks.append(("42_even_n_reducible_full_chirality", {
        "schema_version": "LIE-TRAIN-CERT-1.0",
        "cert_id": "forged_so8_full_k4",
        "family": "so2n_majorana",
        "n_qubits": 4,
        "dla_dimension": 28,
        "threshold_tau": "1/100",
        "verdict": "CERTIFIED_TRAINABLE",
        "exact_variance": "12/245",
        "sectors": [{
            "k_degree": 4, "chirality": "full", "irrep_dim": 70, "gate_energy": "1",
            "casimir_ad_norm_sq": "40", "rho_proj_norm_sq": "3/8", "M_proj_norm_sq": "16",
            "sector_variance": "12/245",
        }],
        "physical_checks": {"total_rho_proj_norm_sq": "3/8", "rho_purity_bound_ok": 1, "all_sectors_distinct": 1},
    }))
    # 43 (H17): Inconsistent gate_energy across multi-sector SO(2n) certificate (arithmetically consistent otherwise!)
    def _attack_43(c):
        c["sectors"][1]["gate_energy"] = "2"
        cas1 = Fraction(c["sectors"][1]["casimir_ad_norm_sq"]) * 2
        v1 = Fraction(c["sectors"][1]["sector_variance"]) * 2
        v0 = Fraction(c["sectors"][0]["sector_variance"])
        c["sectors"][1]["casimir_ad_norm_sq"] = str(cas1)
        c["sectors"][1]["sector_variance"] = str(v1)
        c["exact_variance"] = str(v0 + v1)

    add("43_inconsistent_so2n_gate_energy_multisector", _attack_43, base=c_ms)
    # 44 (H24): Inconsistent gate_energy across chiral plus/minus sectors
    add("44_inconsistent_so2n_gate_energy_chiral", lambda c: (
        c["sectors"][1].__setitem__("gate_energy", "100"),
        c["sectors"][1].__setitem__("casimir_ad_norm_sq", "2000"),
    ), base=c_ch)
    # 45 (H03): U(1) per-sector purity violation ||rho_w^o||^2 > 1 - 1/d_w
    attacks.append(("45_u1_per_sector_purity_exceeded", {
        "schema_version": "LIE-TRAIN-CERT-1.0",
        "cert_id": "forged_u1_sector_purity",
        "family": "u1_hamming",
        "n_qubits": 4,
        "dla_dimension": 65,
        "threshold_tau": "1/100",
        "verdict": "CERTIFIED_TRAINABLE",
        "exact_variance": "1/5",
        "sectors": [{
            "hamming_weight": 1, "block_dim": 4, "irrep_dim": 15,
            "H_traceless_norm_sq": "2", "casimir_ad_norm_sq": "16",
            "rho_traceless_norm_sq": "15/16", "M_traceless_norm_sq": "3",
            "sector_variance": "1/5",
        }],
        "physical_checks": {"total_rho_proj_norm_sq": "15/16", "rho_purity_bound_ok": 1, "all_sectors_distinct": 1},
    }))
    # 46 (H04): U(1) multi-sector weighted purity violation sum ||rho_w^o||^2 / (1 - 1/d_w) = 71/60 > 1
    attacks.append(("46_u1_multisector_weighted_purity_exceeded", {
        "schema_version": "LIE-TRAIN-CERT-1.0",
        "cert_id": "forged_u1_multisector_purity",
        "family": "u1_hamming",
        "n_qubits": 4,
        "dla_dimension": 65,
        "threshold_tau": "1/100",
        "verdict": "CERTIFIED_TRAINABLE",
        "exact_variance": "31/147",
        "sectors": [
            {
                "hamming_weight": 1, "block_dim": 4, "irrep_dim": 15,
                "H_traceless_norm_sq": "2", "casimir_ad_norm_sq": "16",
                "rho_traceless_norm_sq": "7/16", "M_traceless_norm_sq": "3",
                "sector_variance": "7/75",
            },
            {
                "hamming_weight": 2, "block_dim": 6, "irrep_dim": 35,
                "H_traceless_norm_sq": "4", "casimir_ad_norm_sq": "48",
                "rho_traceless_norm_sq": "1/2", "M_traceless_norm_sq": "6",
                "sector_variance": "144/1225",
            },
        ],
        "physical_checks": {"total_rho_proj_norm_sq": "15/16", "rho_purity_bound_ok": 1, "all_sectors_distinct": 1},
    }))
    # 47 (H05): Nonzero traceless norms on 1D U(1) block w=0 (d_0 = 1)
    attacks.append(("47_u1_1d_block_nonzero_traceless", {
        "schema_version": "LIE-TRAIN-CERT-1.0",
        "cert_id": "forged_u1_1d_block",
        "family": "u1_hamming",
        "n_qubits": 4,
        "dla_dimension": 65,
        "threshold_tau": "1/100",
        "verdict": "CERTIFIED_BLIND_SPOT",
        "exact_variance": "0",
        "sectors": [{
            "hamming_weight": 0, "block_dim": 1, "irrep_dim": 0,
            "H_traceless_norm_sq": "999", "casimir_ad_norm_sq": "0",
            "rho_traceless_norm_sq": "15/16", "M_traceless_norm_sq": "999",
            "sector_variance": "0",
        }],
        "physical_checks": {"total_rho_proj_norm_sq": "15/16", "rho_purity_bound_ok": 1, "all_sectors_distinct": 1},
    }))
    # 48 (H28): Control character / newline injected in cert_id
    add("48_cert_id_embedded_newline", lambda c: c.__setitem__("cert_id", "bad\nid"))

    assert len(attacks) == 48
    for name, cand in attacks:
        with pytest.raises(ValueError, match="."):
            verify_lie_trainability_certificate_exact(cand)
