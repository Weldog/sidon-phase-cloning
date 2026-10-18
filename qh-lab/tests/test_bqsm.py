import numpy as np, importlib.util, pathlib

def test_shift_masking_gives_no_bqsm_protection():
    spec = importlib.util.spec_from_file_location("atk", pathlib.Path(__file__).parent.parent / "bqsm_mask_attack.py")
    src = (pathlib.Path(__file__).parent.parent / "bqsm_mask_attack.py").read_text().split("out = []")[0]
    ns = {}; exec(src, ns)
    rng = np.random.default_rng(3)
    K = np.sort(rng.choice(1024, size=8, replace=False))
    assert ns["run"](K, 1024, 64, trials=60) > 0.85      # zero quantum memory, near-certain key recovery
