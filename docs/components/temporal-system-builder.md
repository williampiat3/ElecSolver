# TemporalSystemBuilder

This class models **time-dependent** systems using resistors, capacitors, coils, and mutuals.

## Features

- Supports tension and intensity sources
- Models inductive and resistive mutuals
- Detects and couples multiple subsystems
- Accepts resistances, capacities, and coils
- Constructs sparse linear systems in COO format

## Example

We would like to study the following system:

![Temporal system](../img/schema2.png)

With `R=1`, `L=0.1`, `C=2` this gives:

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import spsolve
from ElecSolver import TemporalSystemBuilder

## Defining resistances
res_coords = np.array([[0, 2], [1, 3]], dtype=int)
res_data = np.array([1, 1], dtype=float)

## Defining coils
coil_coords = np.array([[1, 0], [3, 2]], dtype=int)
coil_data = np.array([0.1, 0.1], dtype=float)

## Defining capacities
capa_coords = np.array([[1, 3], [2, 0]], dtype=int)
capa_data = np.array([2, 2], dtype=float)

## Defining empty mutuals here
mutuals_coords = np.array([[], []], dtype=int)
mutuals_data = np.array([], dtype=float)

res_mutuals_coords = np.array([[], []], dtype=int)
res_mutuals_data = np.array([], dtype=float)

## Initializing system
elec_sys = TemporalSystemBuilder(
    coil_coords,
    coil_data,
    res_coords,
    res_data,
    capa_coords,
    capa_data,
    mutuals_coords,
    mutuals_data,
    res_mutuals_coords,
    res_mutuals_data,
)

## Add source
elec_sys.add_current_source(10, 1, 0)

## Setting ground at point 0
elec_sys.set_ground(0)

## Build system
elec_sys.build_system()

# Getting initial condition system
S_i, b = elec_sys.get_init_system()
sol = spsolve(S_i.tocsr(), b)

# Get system (S1 is real part, S2 derivative part)
S1, S2, rhs = elec_sys.get_system()

## Solving using implicit Euler scheme
dt = 0.08
vals_res1 = []
vals_res2 = []

for _ in range(50):
    temporal_response = elec_sys.build_intensity_and_voltage_from_vector(sol)
    vals_res1.append(temporal_response.intensities_res[1])
    vals_res2.append(temporal_response.intensities_res[0])
    sol = spsolve(S2 + dt * S1, b * dt + S2 @ sol)

plt.xlabel("Time")
plt.ylabel("Intensity")
plt.plot(vals_res1, label="intensity res 1")
plt.plot(vals_res2, label="intensity res 2")
plt.legend()
plt.savefig("intensities_res.png")
```

This outputs the following graph that displays the intensity passing through the resistances:

![Intensity through resistances](../img/intensities_res.png)
