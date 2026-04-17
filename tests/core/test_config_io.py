from pathlib import Path
from famrdp_jax.core.config_io import load_param_inp, build_gas_model

FIXTURE = Path(__file__).parent.parent / "fixtures" / "minimal.inp"

def test_load_param_inp_returns_nested_dict():
    raw = load_param_inp(FIXTURE)
    assert raw["general"]["casname"].strip() == "test"
    assert raw["inflow"]["moo"] == 0.3
    assert raw["method"]["nflux"] == 4

def test_build_gas_model_from_param():
    raw = load_param_inp(FIXTURE)
    gas = build_gas_model(raw)
    assert gas.gamma == 1.4
    assert gas.mu0 == 1.71608e-5
