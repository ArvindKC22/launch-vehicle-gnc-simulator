# Module 1 — dynamics core

## State and conventions

The state is

`x = [r_I, v_I, q_BI, omega_B, m_prop]`,

where `r_I` and `v_I` are inertial position and velocity, `q_BI` is a unit
quaternion rotating body vectors into inertial coordinates, `omega_B` is body
angular rate, and `m_prop` is remaining propellant mass. Total mass is dry mass
plus propellant mass. The current implementation uses a linear interpolation
between full and dry inertia; this is an explicit first-order model until a
tank/propellant geometry model is added. Thrust-vector moment is computed from
the CG-to-nozzle application-point vector. Aerodynamic moment tables and a
detailed moving-CG model remain explicit follow-on additions.

Forces are evaluated in body axes and rotated into inertial axes. Gravity is a
central inverse-square field. Aerodynamic force uses interpolated drag versus
Mach and lift versus angle of attack. Motor thrust is interpolated in time and
reduced by nozzle-exit-area times ambient pressure. Wind is interpolated versus
altitude. The atmosphere is the 1976 standard-atmosphere model through 20 km,
with the 20–200 km continuation documented as a temporary low-fidelity model.

## Quaternion normalization

The quaternion is normalized at every derivative evaluation and after each RK4
step. This prevents roundoff drift from changing the rotation matrix while
preserving fourth-order integration of the underlying differential equation.
The state is never renormalized component-by-component; the complete
four-vector is projected back to the unit 3-sphere. Later, the estimator will
use a multiplicative attitude error and will not average quaternions directly.

## Current validation

`tests/dynamics_test.cpp` validates a 600 s unpowered, drag-free ballistic arc
using specific orbital energy

`E = |v|²/2 − μ/|r|`.

With a fixed 0.1 s RK4 step, the acceptance bound is `1e-9` relative energy
error. The bound is stated against the chosen timestep and central-gravity
model; it is not a claim of physical-model accuracy. The test also verifies
unit quaternion preservation.

## Deliberate limits before integration with the rest of the simulator

Thrust misalignment, CG-dependent moments, six-component aerodynamic moments,
and a complete high-altitude atmosphere are not silently invented here. The
force/moment interfaces are kept explicit so those models can be added with
their own validation tests in subsequent work.
