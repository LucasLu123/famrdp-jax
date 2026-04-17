from pathlib import Path
import pytest

ROOT = Path(__file__).parent.parent.parent
REF_DIR = ROOT / "validation/references/test3/step"
PARAM = ROOT / "param.inp"

@pytest.mark.skipif(not REF_DIR.exists() or not PARAM.exists(), reason="reference/param 未就位")
def test_ten_steps_match():
    pass  # 实际测试逻辑留待 reference 生成后填写
