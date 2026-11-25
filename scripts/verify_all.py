"""Gate-2 end-to-end verification harness for Paper 09.

Author: Lluis Eriksson
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import verify_lie_trainability_certificate_exact

EXPECTED_LEAN_THEOREMS = [
    "LieTrainability.nchoose_absorption",
    "LieTrainability.nchoose_symm",
    "LieTrainability.majorana_casimir_binomial_identity",
    "LieTrainability.pow2_sum_nchoose",
    "LieTrainability.nchoose_le_succ_left",
    "LieTrainability.nchoose_le_central",
    "LieTrainability.central_binom_exponential_bound",
    "LieTrainability.hodge_chiral_doubling",
    "LieTrainability.hodge_chiral_orthogonality_zero",
    "LieTrainability.hodge_chiral_balanced_recovers_base",
    "LieTrainability.u1_traceless_norm_zero_iff_scalar",
    "LieTrainability.u1_gate_commutant_blind_spot",
    "LieTrainability.u1_active_sector_variance_pos",
    "LieTrainability.irrep_variance_sum_nonneg",
    "LieTrainability.mutual_exclusivity_trainable_barren",
    "LieTrainability.blind_spot_excludes_trainable",
    "LieTrainability.chebyshev_sample_complexity_bound",
    "LieTrainability.matchgate_poly_vs_exp_separation",
]

FORBIDDEN_AXIOMS = ("sorryAx", "Classical.choice", "Lean.ofReduceBool", "Lean.trustCompiler")

MANIFEST_FILES = [
    "LICENSE",
    "README.md",
    "pyproject.toml",
    "lean/lean-toolchain",
    "lean/lakefile.toml",
    "lean/LieTrainability.lean",
    "lean/AuditAxioms.lean",
    "src/__init__.py",
    "src/lie_algebra_exact.py",
    "src/lie_train_certifier.py",
    "scripts/run_all_experiments.py",
    "scripts/make_figures.py",
    "scripts/verify_all.py",
    "tests/test_suite.py",
    "results/exp1_certified_atlas.json",
    "results/exp2_haar_mc_crosscheck.json",
    "results/exp3_phase_transitions_scaling.json",
    "results/exp4_sample_complexity_separation.json",
    "results/certificates/so2n_n20_k2_trainable.cert.json",
    "results/certificates/so2n_n10_k4_outside_dla_trainable.cert.json",
    "results/certificates/so2n_n20_kn_barren_in_poly_dla.cert.json",
    "results/certificates/so2n_n50_kn_barren_50q.cert.json",
    "results/certificates/so2n_n4_hodge_chiral_plus_2x.cert.json",
    "results/certificates/so2n_n4_hodge_chiral_blind_spot.cert.json",
    "results/certificates/so2n_n12_multisector_k2_k4.cert.json",
    "results/certificates/u1_n20_w1_trainable.cert.json",
    "results/certificates/u1_n30_half_filling_barren.cert.json",
    "results/certificates/u1_n12_scalar_gate_blind_spot.cert.json",
    "results/certificates/su_n20_unstructured_barren.cert.json",
    "results/certificates/so2n_n100_k2_large_scale.cert.json",
    "paper/figures/fig1_phase_transitions.png",
    "paper/figures/fig2_haar_mc_hodge_chiral.png",
    "paper/figures/fig3_sample_complexity_separation.png",
    "paper/main.tex",
    "paper/references.bib",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def check_lean() -> None:
    lean_dir = ROOT / "lean"
    r_build = subprocess.run(["lake", "build"], cwd=lean_dir, capture_output=True, text=True, check=False)
    if r_build.returncode != 0:
        raise RuntimeError(f"lake build failed:\n{r_build.stdout}\n{r_build.stderr}")
    if "warning:" in (r_build.stdout + r_build.stderr).lower():
        raise RuntimeError(f"Lean build emitted warnings:\n{r_build.stdout}\n{r_build.stderr}")

    r_audit = subprocess.run(
        ["lake", "env", "lean", "AuditAxioms.lean"],
        cwd=lean_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if r_audit.returncode != 0:
        raise RuntimeError(f"AuditAxioms.lean failed:\n{r_audit.stdout}\n{r_audit.stderr}")
    out = r_audit.stdout + "\n" + r_audit.stderr
    for thm in EXPECTED_LEAN_THEOREMS:
        if f"'{thm}'" not in out:
            raise RuntimeError(f"Missing Lean theorem in axiom audit: {thm}")
    for bad in FORBIDDEN_AXIOMS:
        if bad in out:
            raise RuntimeError(f"Forbidden Lean axiom detected: {bad}\n{out}")
    print(f"[PASS] Lean 4 kernel: {len(EXPECTED_LEAN_THEOREMS)}/{len(EXPECTED_LEAN_THEOREMS)} theorems, 0 warnings, 0 sorry, 0 Classical.choice")


def check_certificates() -> None:
    cert_dir = ROOT / "results" / "certificates"
    certs = sorted(cert_dir.glob("*.cert.json"))
    if len(certs) != 12:
        raise RuntimeError(f"Expected 12 certificate files, found {len(certs)}")
    for cp in certs:
        data = json.loads(cp.read_text(encoding="utf-8"))
        res = verify_lie_trainability_certificate_exact(data)
        if not res.get("verified"):
            raise RuntimeError(f"Certificate failed verification: {cp.name}")
    print(f"[PASS] Exact-Q certificates: {len(certs)}/{len(certs)} independently verified (0 floats)")


def check_pytest() -> None:
    r = subprocess.run(
        [sys.executable, "-m", "pytest", "-q"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if r.returncode != 0:
        raise RuntimeError(f"pytest failed:\n{r.stdout}\n{r.stderr}")
    print(f"[PASS] pytest: {r.stdout.strip().splitlines()[-1]}")


def write_or_check_manifest(write_mode: bool) -> None:
    manifest_path = ROOT / "SHA256SUMS.json"
    computed: Dict[str, str] = {}
    missing: List[str] = []
    for rel in MANIFEST_FILES:
        p = ROOT / rel
        if not p.exists():
            missing.append(rel)
        else:
            computed[rel] = sha256_file(p)
    if missing:
        raise RuntimeError(f"Missing manifest files: {missing}")

    if write_mode:
        manifest_path.write_text(json.dumps(computed, indent=2) + "\n", encoding="utf-8")
        print(f"[PASS] Wrote SHA256SUMS.json ({len(computed)} files)")
        return

    if not manifest_path.exists():
        raise RuntimeError("SHA256SUMS.json not found; run with --write-manifest first")
    recorded = json.loads(manifest_path.read_text(encoding="utf-8"))
    if recorded != computed:
        diffs = [k for k in set(recorded) | set(computed) if recorded.get(k) != computed.get(k)]
        raise RuntimeError(f"SHA256SUMS.json mismatch on: {diffs}")
    print(f"[PASS] SHA256SUMS.json: {len(computed)}/{len(computed)} hashes match")


def main() -> None:
    parser = argparse.ArgumentParser(description="Gate-2 verification harness for Paper 09")
    parser.add_argument("--write-manifest", action="store_true")
    parser.add_argument("--skip-manifest", action="store_true")
    args = parser.parse_args()

    check_lean()
    check_certificates()
    check_pytest()
    if not args.skip_manifest:
        write_or_check_manifest(args.write_manifest)
    print("=== ALL PAPER 09 GATE-2 VERIFICATION CHECKS PASSED ===")


if __name__ == "__main__":
    main()
