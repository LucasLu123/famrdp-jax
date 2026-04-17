"""Core data types for famrdp-jax.

All array-carrying types (Block, Metrics, State) are registered as JAX pytrees
so that jax.tree_util.tree_map / jit / vmap work transparently.

Non-array config types (BlockTopology, GasModel, SchemeChoice, Config) use
frozen dataclasses; they are treated as static / aux data.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Optional, Tuple

import jax
import jax.numpy as jnp

from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.core.dtypes import check_array_f64


# ---------------------------------------------------------------------------
# Block
# ---------------------------------------------------------------------------

@dataclass
class Block:
    """Per-block conservative-variable array + coordinate array.

    q   : shape (5, ni, nj, nk) – conservative vars [rho, rho*u, rho*v, rho*w, rho*E]
    xyz : shape (3, ni, nj, nk) – Cartesian coordinates
    Both must be float64.
    """
    q:   jnp.ndarray
    xyz: jnp.ndarray

    def __post_init__(self):
        check_array_f64({"q": self.q, "xyz": self.xyz})

    # pytree protocol
    def tree_flatten(self):
        leaves = (self.q, self.xyz)
        aux = {}
        return leaves, aux

    @classmethod
    def tree_unflatten(cls, aux, leaves):
        # Bypass __post_init__ validation on unflatten (arrays already validated)
        obj = object.__new__(cls)
        object.__setattr__(obj, "q",   leaves[0])
        object.__setattr__(obj, "xyz", leaves[1])
        return obj


jax.tree_util.register_pytree_node_class(Block)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

@dataclass
class Metrics:
    """Grid metrics for one block.

    jac  : Jacobian,          shape (ni, nj, nk)
    kxyz : metric derivatives, shape (3, 3, ni, nj, nk)  – kxyz[i,j] = d xi_i / d x_j
    vol  : cell volume,        shape (ni, nj, nk)
    All must be float64.
    """
    jac:  jnp.ndarray
    kxyz: jnp.ndarray
    vol:  jnp.ndarray

    def __post_init__(self):
        check_array_f64({"jac": self.jac, "kxyz": self.kxyz, "vol": self.vol})

    def tree_flatten(self):
        leaves = (self.jac, self.kxyz, self.vol)
        aux = {}
        return leaves, aux

    @classmethod
    def tree_unflatten(cls, aux, leaves):
        obj = object.__new__(cls)
        object.__setattr__(obj, "jac",  leaves[0])
        object.__setattr__(obj, "kxyz", leaves[1])
        object.__setattr__(obj, "vol",  leaves[2])
        return obj


jax.tree_util.register_pytree_node_class(Metrics)


# ---------------------------------------------------------------------------
# State
# ---------------------------------------------------------------------------

@dataclass
class State:
    """Global simulation state.

    blocks  : tuple of Block   (one per block)
    metrics : tuple of Metrics (one per block)
    t       : simulation time  (scalar float, treated as leaf)
    step    : time-step index  (stored in aux, not a leaf)
    """
    blocks:  tuple
    metrics: tuple
    t:       float = 0.0
    step:    int = 0

    def tree_flatten(self):
        # Flatten blocks and metrics as sequences; t is a numeric leaf.
        # step is metadata (int) – kept in aux so that tree_map doesn't touch it.
        leaves = (*self.blocks, *self.metrics, self.t)
        aux = {
            "n_blocks":  len(self.blocks),
            "n_metrics": len(self.metrics),
            "step":      self.step,
        }
        return leaves, aux

    @classmethod
    def tree_unflatten(cls, aux, leaves):
        n_blocks  = aux["n_blocks"]
        n_metrics = aux["n_metrics"]
        blocks  = tuple(leaves[:n_blocks])
        metrics = tuple(leaves[n_blocks : n_blocks + n_metrics])
        t       = leaves[n_blocks + n_metrics]
        obj = object.__new__(cls)
        object.__setattr__(obj, "blocks",  blocks)
        object.__setattr__(obj, "metrics", metrics)
        object.__setattr__(obj, "t",       t)
        object.__setattr__(obj, "step",    aux["step"])
        return obj


jax.tree_util.register_pytree_node_class(State)


# ---------------------------------------------------------------------------
# Non-array config / topology types  (frozen dataclasses, not pytrees)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BlockTopology:
    """Connectivity and BC info for one block.

    block_id  : index of this block in the global list
    neighbors : dict[Face -> None | (block_id, Face, orientation_flag)]
    bc_type   : dict[Face -> BCType]
    """
    block_id:  int
    neighbors: dict
    bc_type:   dict


@dataclass(frozen=True)
class GasModel:
    """Ideal-gas / Sutherland parameters."""
    gamma:   float = 1.4
    R:       float = 287.058       # J/(kg·K)
    mu0:     float = 1.71608e-5    # Pa·s  at T0
    T0:      float = 273.15        # K
    Ts:      float = 110.4         # K  (Sutherland constant)
    Pr_lam:  float = 0.72


@dataclass(frozen=True)
class SchemeChoice:
    """Numerical scheme selections."""
    reconstruct: str
    flux:        str
    rk_order:    int
    limiter:     str


@dataclass(frozen=True)
class Config:
    """Top-level simulation configuration."""
    gas:       GasModel
    scheme:    SchemeChoice
    bc_params: dict
    topology:  tuple
    ghost:     int   = 2
    nan_check: bool  = True
    debug:     bool  = False
    cfl:       float = 1.0
