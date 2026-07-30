"""GWAS preparation, GenoPred orchestration and PRS validation."""

from .gwas import TraitSpec, load_trait_specs, standardise_trait

__all__ = ["TraitSpec", "load_trait_specs", "standardise_trait"]
