"""Fortran namelist (param.inp) parser and Config builders.

Uses f90nml if available; falls back to a minimal regex-based parser
that handles the subset of namelist syntax used by famrdp.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict

from famrdp_jax.core.types import GasModel

# ---------------------------------------------------------------------------
# Optional f90nml
# ---------------------------------------------------------------------------
try:
    import f90nml
    HAS_F90NML = True
except ImportError:
    HAS_F90NML = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_param_inp(path) -> Dict[str, Dict[str, Any]]:
    """Parse a Fortran namelist file and return a nested dict.

    Keys are lower-cased group and variable names.
    Values are Python int / float / str as appropriate.
    """
    path = Path(path)
    if HAS_F90NML:
        nml = f90nml.read(str(path))
        return {k.lower(): dict(v) for k, v in nml.items()}
    else:
        return _parse_namelist(path)


def build_gas_model(raw: Dict[str, Dict[str, Any]]) -> GasModel:
    """Construct a GasModel from the parsed namelist dict."""
    g = raw["gasmodel"]
    R = 8314.46 / g["wgas"]
    return GasModel(
        gamma=g["gamma"],
        R=R,
        mu0=g["mu0sth"],
        T0=g["t0sth"],
        Ts=g["tssth"],
        Pr_lam=g["prlam"],
    )


# ---------------------------------------------------------------------------
# Fallback parser (no f90nml dependency)
# ---------------------------------------------------------------------------

def _cast(val_str: str) -> Any:
    """Try to cast a namelist value string to int, float, or str."""
    s = val_str.strip()
    if s.startswith('"') or s.startswith("'"):
        return s.strip("\"'")
    # Fortran logical
    if s.lower() in (".true.", "t"):
        return True
    if s.lower() in (".false.", "f"):
        return False
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    return s


def _parse_namelist(path: Path) -> Dict[str, Dict[str, Any]]:
    """Minimal Fortran namelist parser for famrdp param.inp files."""
    text = path.read_text(encoding="utf-8")

    result: Dict[str, Dict[str, Any]] = {}
    # Match each &groupname ... / block (non-greedy, DOTALL)
    for m in re.finditer(r"&(\w+)(.*?)/", text, re.DOTALL):
        group = m.group(1).lower()
        body = m.group(2)
        entries: Dict[str, Any] = {}
        # Each key = value pair, value ends before comma or newline
        for kv in re.finditer(r"(\w+)\s*=\s*([^,\n/]+)", body):
            key = kv.group(1).lower()
            val_str = kv.group(2).strip().rstrip(",").strip()
            entries[key] = _cast(val_str)
        result[group] = entries
    return result
