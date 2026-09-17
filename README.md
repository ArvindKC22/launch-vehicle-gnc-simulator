# Launch Vehicle GNC Simulator

A C++17 six-degree-of-freedom launch-vehicle guidance, navigation, and control (GNC) simulation environment. The implementation couples rigid-body flight dynamics, dispersed sensors, delayed measurements, a UKF-based navigation estimator, attitude control, and Monte Carlo analysis.

## Scope

The repository contains the reproducible source code, configuration, unit tests, compact Monte Carlo outputs, figures, validation notes, and manuscript artifacts for the project. Large per-step trace files and compiled build products are intentionally excluded from GitHub; the public archival release is available at [Zenodo DOI: 10.5281/zenodo.22801065](https://doi.org/10.5281/zenodo.22801065).

## Build and test

Dependencies:

- CMake 3.20 or later
- C++17 compiler
- Eigen3
- GoogleTest
- OpenMP (optional)

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Monte Carlo execution

```bash
./build/run_monte_carlo results/mc_results.csv . 1000 600 0.01 4
python3 analysis/plot_results.py results/mc_results.csv figures/example
```

The YAML files in `config/` define vehicle and dispersion parameters. Command-line settings may override defaults. See `docs/validation_report.md` for campaign interpretation, limitations, and reproducibility notes.

## Repository contents

- `src/`, `include/`: C++ implementation of dynamics, sensors, estimation, control, and Monte Carlo execution
- `config/`: vehicle and dispersion configuration
- `tests/`: automated verification suite
- `analysis/`: result plotting and external trajectory comparison utilities
- `results/`, `figures/`: compact outputs and figures from validation campaigns
- `docs/`: validation report, technical report, and submission manuscript

## Status and limitations

The simulator is a research and educational tool, not flight-certified software. The included validation describes the tested scenarios and remaining requirements, including independent trajectory validation against OpenRocket or RASAero and further covariance tuning.

## Citation

If the archived dataset or release is used, cite:

> Chandirakala, Arvind Kanagasabapathi. *Launch Vehicle GNC Simulator: Reproducible Software, Configuration, and Monte Carlo Validation Dataset*. Zenodo. https://doi.org/10.5281/zenodo.22801065

## Licensing

- The software source code is available under the [MIT License](LICENSE).
- Monte Carlo datasets, figures, configuration datasets, and documentation are available under [Creative Commons Attribution 4.0 International](LICENSE-DATA.md).
- The Zenodo archive carries the same licensing information.
