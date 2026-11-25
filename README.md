# Exact-Q Adjoint-Irrep Variance Laws and Machine-Checked Trainability Certificates

Machine-checked Lean 4 formalization, zero-floating-point-trust exact-$\mathbb{Q}$ certificate verifier (`LIE-TRAIN-CERT-1.0`), and reproducible computational benchmarks for:

> **Exact-$\mathbb{Q}$ Adjoint-Irrep Variance Laws, Majorana Hodge-Chiral Selection Rules, and Machine-Checked Trainability Certificates for Compact Lie-Group Quantum Circuits**
> — *Lluis Eriksson* (November 2026)

Directly advancing Xie, Szpankowski, and Grama (`arXiv:2611.18378`, Nov 2026) and Ragone et al. (*Nature Communications* **15**, 7172, 2024).

## Key Results

1. **Exact Rational Adjoint-Irrep Casimir–Weingarten Identity (Theorem 3.1):**
   Tightens Xie–Szpankowski–Grama's sufficient trainability condition into an exact closed-form rational identity $\operatorname{Var}_{G\times G}(\partial_\theta C) = \sum_\lambda \frac{\|\operatorname{ad}_H|_{V_\lambda}\|_{\mathrm{HS}}^2}{\dim(V_\lambda)^2} \|\Pi_{V_\lambda}(\rho)\|_{\mathrm{HS}}^2 \|\Pi_{V_\lambda}(M)\|_{\mathrm{HS}}^2$ over irreducible $\operatorname{Ad}_G$-representations (including observables $M \notin i\mathfrak{g}$ outside the dynamical Lie algebra).
2. **Matchgate $SO(2n)$ Majorana-Degree Phase Transition & Even-$n$ Hodge-Chiral Selection Law (Theorem 3.2):**
   Proves $\|\operatorname{ad}_H|_{\Lambda^k}\|_{\mathrm{HS}}^2 = 2 E_H \binom{2n-2}{k-1}$, yielding $\Theta(n^{-2})$ trainability at $k=2$ ($Z_1 \in i\mathfrak{so}(2n)$), $\Theta(n^{-3})$ trainability at $k=4$ ($Z_1 Z_2 \notin i\mathfrak{so}(2n)$, **outside** $i\mathfrak{so}(2n)$), and a provable $\Theta(2^{-n})$ exponential barren plateau at $k=n$ within a circuit whose DLA $\dim(\mathfrak{so}(2n)) = n(2n-1)$ is strictly quadratic. Characterizes the complementary-degree parity isomorphism $\Lambda^k \cong \Lambda^{2n-k}$ ($k < n$) and, for even $n$, the real Hodge-star split $\Lambda^n = \Lambda^n_+ \oplus \Lambda^n_-$ giving exact $2\times$ amplification for aligned chiralities ($M_+ = (Z_1\cdots Z_{n/2} + Z_{n/2+1}\cdots Z_n)/\sqrt{2}$ on $|0\rangle^{\otimes n}$) and an exact $0$-variance chiral blind spot for opposite chiralities ($M_-$).
3. **$U(1)$ Hamming-Weight Filling-Fraction Transition & Gate-Commutant Blind-Spot Theorem (Theorem 3.3):**
   Exact rational variance in every Hamming-weight block $w$ ($d_w = \binom{n}{w}$), separating $\Theta(n^{-3})$ ($w=1$) and $\Theta(n^{-4})$ ($w=2$) trainability from $\Theta(\sqrt{n}\, 2^{-n})$ half-filling barren plateaus, plus exact zero variance when $H$ acts as a scalar on active blocks (e.g., $H = \sum Z_j$).
4. **Exponential Sample-Complexity Separation & Zero-Float `LIE-TRAIN-CERT-1.0` Verifier (Theorem 3.4):**
   Certifies `CERTIFIED_TRAINABLE`, `CERTIFIED_BARREN_PLATEAU`, and `CERTIFIED_BLIND_SPOT` in $< 0.5$\,ms up to $n=100$ qubits with zero quantum shots, vs. $\Omega(\tau^{-2}) \sim 9.49\times 10^{24}$ finite-shot pairs at $n=40$ for sampling estimators.

## Verification & Reproducibility

```bash
python scripts/run_all_experiments.py
python scripts/make_figures.py
python scripts/verify_all.py
```

- Lean 4 (`v4.35.0`): 18/18 theorems, 0 warnings, 0 `sorry`, 0 `Classical.choice`.
- Pure-stdlib exact-$\mathbb{Q}$ certificates: 12/12 verified (0 floats).
- Test suite: 48/48 live adversarial red-team attacks rejected.

## License

MIT License (see `LICENSE`).
