"""Paper 09: Exact-Q Adjoint-Irrep Variance Laws and Machine-Checked Trainability Certificates."""

from .lie_algebra_exact import (
    exact_sample_complexity_lower_bound,
    so2n_casimir_ad_norm_sq,
    so2n_majorana_irrep_dim,
    so2n_sector_variance_exact,
    su_unstructured_variance_exact,
    u1_sector_variance_exact,
)
from .lie_train_certifier import (
    SCHEMA_VERSION,
    build_lie_trainability_certificate,
    canonical_fraction_str,
    expected_dla_dim,
    parse_canonical_fraction,
    verify_lie_trainability_certificate_exact,
)

__all__ = [
    "SCHEMA_VERSION",
    "build_lie_trainability_certificate",
    "canonical_fraction_str",
    "exact_sample_complexity_lower_bound",
    "expected_dla_dim",
    "parse_canonical_fraction",
    "so2n_casimir_ad_norm_sq",
    "so2n_majorana_irrep_dim",
    "so2n_sector_variance_exact",
    "su_unstructured_variance_exact",
    "u1_sector_variance_exact",
    "verify_lie_trainability_certificate_exact",
]
