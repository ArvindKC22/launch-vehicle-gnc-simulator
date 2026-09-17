# Validation status

## Automated verification

The current build passes all 14 CTest cases, including dynamics invariants,
sensor timing and latency, actuator limits, MRP shadow-set handling, EKF/UKF
interfaces, asynchronous IMU handling, GPS latency compensation, and the
linear UKF analytic benchmark, and an explicit specific-force navigation
propagation check.

## Navigation and covariance

The UKF now propagates position, velocity, and attitude from the latest IMU
specific-force and angular-rate samples. Accelerometer and gyro biases are
explicit error states, sensor turn-on/random-walk bias is modeled, propellant
depletion is retained in the IMU branch for powered/coast scheduling, delayed
GPS is evaluated against a state history and propagated to the current epoch,
and covariance uses a Joseph update.

The deterministic navigation-only mode is useful for checking propagation:

```text
./build/run_monte_carlo /tmp/navigation.csv . 1 30 0.002 1 0 /tmp/navigation_trace.csv 0.5 1
```

The dispersed closed-loop filter still needs additional tuning before it is
flight-quality. The reproducible 1,000-run, 60-second campaign in
`results/formal_1000.csv` produced 70/1,000 impacts within the short window;
the summary had mean per-run NEES about 160, mean GPS NIS about 3.48, and mean
barometer NIS about 0.91. A separate 20-run, 600-second check produced 20/20
valid impacts with mean NEES about 16.9, GPS NIS about 0.78, and barometer NIS
about 1.04. These results show substantial improvement but are not yet a
flight-certification claim.

Trace audit shows the rare NEES outliers are concentrated in the 0--30 s
powered phase. Powered gyro-model sigmas of 0.10, 0.25, and 0.50 rad/s were
compared; 0.50 rad/s was retained because it reduced the full-flight NEES
without covariance failures. The final 1,000-run, 600-second, 10 ms campaign
produced 1,000/1,000 valid impacts, mean per-run NEES 12.2 (median 4.87),
mean GPS NIS 0.73, and mean barometer NIS 0.97.

The completed full-flight 1,000-run, 10 ms campaign is preserved in
`results/formal_1000_full.csv` and `results/formal_1000_full_trace.csv`.
It produced 1,000/1,000 valid impacts, mean per-run NEES about 55.0, mean GPS
NIS about 1.51, and mean barometer NIS about 0.97.

A 1,000-run, 600-second, 20 ms timestep sensitivity campaign was also
executed. It completed all runs but produced mean NEES about 332 and mean GPS
NIS about 7.53; it is evidence that timestep convergence remains required for
the closed-loop controller/navigation implementation.

## YAML dispersion configuration

`config/dispersions.yaml` controls run count, seed, duration, timestep,
thread count, closed-loop mode, consistency cadence, and all current plant and
sensor dispersion sigmas. CLI values remain available as explicit overrides.

## External trajectory comparison

The repository includes `analysis/compare_trajectory.py`. Export a single-row
CSV from OpenRocket or RASAero containing apogee and range, then run:

```text
.venv/bin/python analysis/compare_trajectory.py simulator.csv reference.csv
```

No OpenRocket/RASAero reference export is present in this repository, so an
external agreement claim cannot yet be made.

## Formal campaign artifacts

- `results/formal_1000.csv`
- `results/formal_1000_trace.csv`
- `results/formal_1000_full.csv`
- `results/formal_1000_full_trace.csv`
- `results/formal_1000_final.csv`
- `results/formal_1000_final_trace.csv`
- `figures/formal_1000/apogee_distribution.png`
- `figures/formal_1000/landing_ellipse.png`
- `figures/formal_1000/nees_nis_consistency.png`
- `figures/formal_1000/tornado.png`
- `figures/formal_1000_final/apogee_distribution.png`
- `figures/formal_1000_final/landing_ellipse.png`
- `figures/formal_1000_final/nees_nis_consistency.png`
- `figures/formal_1000_final/tornado.png`
