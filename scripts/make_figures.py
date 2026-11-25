"""Generate 3 publication-quality figures for Paper 09 from exact JSON outputs.

Author: Lluis Eriksson
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES_DIR = ROOT / "results"
FIG_DIR = ROOT / "paper" / "figures"


def style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle="--", alpha=0.35, linewidth=0.7)


def make_fig1(exp3: dict) -> None:
    """Fig 1: Exact-Q trainability vs barren-plateau phase transitions for SO(2n) and U(1)."""
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.6))

    # Panel (a): Matchgate SO(2n) — poly k=2,4 (outside DLA!) vs exp k=n (inside poly DLA!)
    ax = axes[0]
    style_axes(ax)
    rows = exp3["so2n_scaling"]
    ns = [r["n"] for r in rows]
    ax.semilogy(ns, [r["var_k2_float"] for r in rows], "o-", color="#1b7837", lw=2.1, ms=5.5,
                label=r"$k=2$ ($Z_1 \in i\mathfrak{so}(2n)$): $\Theta(n^{-2})$ trainable")
    ax.semilogy(ns, [r["var_k4_float"] for r in rows], "s-", color="#2166ac", lw=2.1, ms=5.5,
                label=r"$k=4$ ($Z_1 Z_2 \notin i\mathfrak{so}(2n)$): $\Theta(n^{-3})$ trainable")
    ax.semilogy(ns, [r["var_kn_chiral_2x_float"] for r in rows], "^--", color="#d95f02", lw=1.7, ms=5.2,
                label=r"$k=n$ Hodge-chiral $\Lambda^n_+$ ($2\times$ exact boost)")
    ax.semilogy(ns, [r["var_kn_float"] for r in rows], "D-", color="#b2182b", lw=2.1, ms=5.2,
                label=r"$k=n$ ($Z_1\cdots Z_{n/2}$): $\Theta(2^{-n})$ barren (quadratic DLA)")
    ax.set_xlabel(r"Number of qubits $n$ ($\dim \mathfrak{so}(2n) = n(2n-1)$)", fontsize=9.8)
    ax.set_ylabel(r"Exact gradient variance $\operatorname{Var}_{SO(2n)}(\partial_\theta C)$", fontsize=9.8)
    ax.set_title(r"(a) Matchgate $SO(2n)$ Majorana-Degree Transition", fontsize=10.3, fontweight="bold")
    ax.set_ylim(1e-13, 1.0)
    ax.legend(loc="lower left", fontsize=7.6, frameon=True, framealpha=0.95)

    # Panel (b): U(1) Hamming-weight sectors vs unstructured SU(2^n)
    ax = axes[1]
    style_axes(ax)
    urows = exp3["u1_scaling"]
    ns_u = [r["n"] for r in urows]
    ax.semilogy(ns_u, [r["var_w1_float"] for r in urows], "o-", color="#1b7837", lw=2.1, ms=5.5,
                label=r"$U(1)$ fixed $w=1$: $\Theta(n^{-3})$ trainable")
    ax.semilogy(ns_u, [r["var_w2_float"] for r in urows], "s-", color="#2166ac", lw=2.1, ms=5.5,
                label=r"$U(1)$ fixed $w=2$: $\Theta(n^{-4})$ trainable")
    ax.semilogy(ns_u, [r["var_half_float"] for r in urows], "D-", color="#d95f02", lw=2.1, ms=5.2,
                label=r"$U(1)$ half-filling $w=n/2$: $\Theta(\sqrt{n}\,2^{-n})$ barren")
    ax.semilogy(ns_u, [r["var_su_float"] for r in urows], "X--", color="#7b3294", lw=1.9, ms=5.5,
                label=r"Unstructured $SU(2^n)$: $\Theta(2^{-n})$ barren")
    ax.set_xlabel(r"Number of qubits $n$", fontsize=9.8)
    ax.set_ylabel(r"Exact gradient variance $\operatorname{Var}(\partial_\theta C)$", fontsize=9.8)
    ax.set_title(r"(b) $U(1)$ Filling-Fraction Transition vs. $SU(2^n)$", fontsize=10.3, fontweight="bold")
    ax.set_ylim(1e-13, 1.0)
    ax.legend(loc="lower left", fontsize=7.6, frameon=True, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_phase_transitions.png", dpi=260, bbox_inches="tight")
    plt.close(fig)


def make_fig2(exp2: dict) -> None:
    """Fig 2: Exact-Q theory vs 12,000-sample Haar Monte Carlo + Hodge-chiral 2x / 0x selection."""
    fig, ax = plt.subplots(figsize=(10.6, 4.4))
    style_axes(ax)

    checks = exp2["checks"]
    labels = [
        r"$SO(6)$" + "\n$k=1$",
        r"$SO(6)$" + "\n$k=2$",
        r"$SO(6)$" + "\n$k=3$",
        r"$SO(6)$" + "\n$k=4$",
        r"$SO(4)$" + "\n" + r"$(+,+)$ $2\times$",
        r"$SO(4)$" + "\n" + r"$(-,-)$ $2\times$",
        r"$SO(4)$" + "\n" + r"$(+,-)$ blind",
        r"$U(1)$" + "\n$SU(6)$",
        r"$SU(8)$" + "\nfull",
    ]
    x = np.arange(len(checks))
    w = 0.36
    ex = [c["exact_float"] for c in checks]
    mc = [c["mc_var"] for c in checks]
    se = [2.0 * c["mc_se"] for c in checks]

    ax.bar(x - w / 2, ex, width=w, color="#2166ac", edgecolor="#111111", linewidth=0.75,
           label=r"Exact-$\mathbb{Q}$ adjoint-irrep formula (Theorems 3.1–3.3)")
    ax.bar(x + w / 2, mc, width=w, yerr=se, capsize=3.5, color="#fdae61", edgecolor="#111111", linewidth=0.75,
           error_kw=dict(ecolor="#222222", lw=1.0),
           label=r"Independent Haar Monte Carlo ($N = 12{,}000$, $\pm 2\sigma$)")

    ax.annotate("Exact chiral\nblind spot (0)", xy=(6.0, 0.004), xytext=(6.0, 0.068),
                ha="center", fontsize=8.0, color="#b2182b", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#b2182b", lw=1.0))

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.2)
    ax.set_ylabel(r"Gradient variance $\operatorname{Var}(\partial_\theta C)$", fontsize=9.8)
    ax.set_ylim(0.0, 0.275)
    ax.set_title(r"Exact-$\mathbb{Q}$ Theory vs. Haar Monte Carlo Across $SO(2n)$, Hodge-Chiral, $U(1)$, and $SU(2^n)$",
                 fontsize=10.2, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8.1, frameon=True, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_haar_mc_hodge_chiral.png", dpi=260, bbox_inches="tight")
    plt.close(fig)


def make_fig3(exp4: dict) -> None:
    """Fig 3: Exponential sampling-cost wall vs O(1) exact-Q algebraic certification."""
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4))

    rows = exp4["records"]
    ns = [r["n"] for r in rows]

    ax = axes[0]
    style_axes(ax)
    ax.semilogy(ns, [r["chebyshev_N_kn"] for r in rows], "D-", color="#b2182b", lw=2.1, ms=5.5,
                label=r"Sampling bound $N_{\mathrm{Cheb}}$: $SO(2n)$ $k=n$ ($\Omega(4^n)$)")
    ax.semilogy(ns, [r["chebyshev_N_u1_half"] for r in rows], "s-", color="#d95f02", lw=2.0, ms=5.2,
                label=r"Sampling bound $N_{\mathrm{Cheb}}$: $U(1)$ half-filling ($\Omega(4^n/n)$)")
    ax.semilogy(ns, [r["chebyshev_N_k2"] for r in rows], "o-", color="#2166ac", lw=2.0, ms=5.2,
                label=r"Sampling bound $N_{\mathrm{Cheb}}$: $SO(2n)$ $k=2$ ($\Omega(n^4)$)")
    ax.axhline(1.0, color="#1b7837", linestyle="--", lw=1.8,
               label=r"Exact-$\mathbb{Q}$ algebraic certificate (0 shots; baseline $y=1$)")
    ax.set_xlabel(r"Number of qubits $n$", fontsize=9.8)
    ax.set_ylabel(r"Required circuit evaluations $N \geq 4/\tau^2$", fontsize=9.8)
    ax.set_title(r"(a) Sampling Detector Cost vs. Zero-Shot $\mathbb{Q}$-Certificate", fontsize=10.1, fontweight="bold")
    ax.set_ylim(0.3, 1e27)
    ax.legend(loc="upper left", fontsize=7.5, frameon=True, framealpha=0.95)

    ax = axes[1]
    style_axes(ax)
    ms = [r["exact_q_verify_ms"] for r in rows]
    ax.bar(range(len(ns)), ms, color="#1b7837", edgecolor="#111111", linewidth=0.75, width=0.58)
    for idx_b, val_b in enumerate(ms):
        ax.text(idx_b, val_b + 0.003, f"{val_b:.2f}", ha="center", va="bottom", fontsize=7.6)
    ax.set_xticks(range(len(ns)))
    ax.set_xticklabels([str(n) for n in ns], fontsize=8.5)
    ax.set_xlabel(r"Number of qubits $n$", fontsize=9.8)
    ax.set_ylabel(r"Exact-$\mathbb{Q}$ verification time (ms)", fontsize=9.8)
    ax.set_title(r"(b) Zero-Float Exact-$\mathbb{Q}$ Verifier Runtime ($n=4\dots 40$)", fontsize=10.1, fontweight="bold")
    ax.set_ylim(0.0, max(0.12, max(ms) * 1.45))

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig3_sample_complexity_separation.png", dpi=260, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    exp2 = json.loads((RES_DIR / "exp2_haar_mc_crosscheck.json").read_text(encoding="utf-8"))
    exp3 = json.loads((RES_DIR / "exp3_phase_transitions_scaling.json").read_text(encoding="utf-8"))
    exp4 = json.loads((RES_DIR / "exp4_sample_complexity_separation.json").read_text(encoding="utf-8"))

    make_fig1(exp3)
    make_fig2(exp2)
    make_fig3(exp4)
    print(f"Generated 3 Paper 09 figures in {FIG_DIR}")


if __name__ == "__main__":
    main()
