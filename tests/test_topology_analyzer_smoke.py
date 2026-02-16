import sys
from pathlib import Path
import numpy as np


def test_topology_analyzer_compute_all_metrics_smoke():
    repo = Path(__file__).resolve().parents[1]
    src = repo / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))

    from fungalmorphospace.analysis.topology_analyzer import TopologyAnalyzer

    # Synthetic periodic pattern (smooth) + noise
    n = 128
    x = np.linspace(0, 2 * np.pi, n, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    pattern = np.sin(6 * X) + np.sin(6 * Y) + 0.1 * np.random.RandomState(0).normal(size=(n, n))

    analyzer = TopologyAnalyzer(pattern, dx=1.0)
    metrics = analyzer.compute_all_metrics()

    # Keys that must exist
    required = [
        "wavelength_fft",
        "wavelength_fft_qc_pass",
        "wavelength_autocorr",
        "wavelength_autocorr_qc_pass",
        "n_components",
        "euler_characteristic",
    ]
    for k in required:
        assert k in metrics

    # Types sanity
    assert isinstance(metrics["wavelength_fft_qc_pass"], (bool, np.bool_))
    assert isinstance(metrics["wavelength_autocorr_qc_pass"], (bool, np.bool_))
