# Extra uses: Hydraulic or Thermal system modeling

ElecSolver can also be used to model mass flow or thermal flux in hydraulic and thermal networks where a pressure or temperature difference can be assimilated to a tension source.

Since electric potentials are always computed relative to the ground node, you may need to rescale the resulting potentials.

We consider the following hydraulic problem:

![Hydraulic system](img/hydraulic.png)

Taking `R=1` gives:

```python
import numpy as np
from scipy.sparse.linalg import spsolve
from ElecSolver import TemporalSystemBuilder

## Defining resistances
R = 1
res_coords = np.array([[0, 2, 1, 0, 1, 3], [1, 3, 3, 2, 2, 0]], dtype=int)
res_data = R * np.array([2, 3, 1, 1, 1, 1], dtype=float)

## Here we are not using coils, capacities or mutuals
coil_coords = np.array([[], []], dtype=int)
coil_data = np.array([], dtype=float)
capa_coords = np.array([[], []], dtype=int)
capa_data = np.array([], dtype=float)
mutuals_coords = np.array([[], []], dtype=int)
mutuals_data = np.array([], dtype=float)
res_mutuals_coords = np.array([[], []], dtype=int)
res_mutuals_data = np.array([], dtype=float)

## Initializing system
hydraulic_sys = TemporalSystemBuilder(
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

## Enforcing a pressure delta of 10 Pa
hydraulic_sys.add_voltage_source(10, 1, 0)
hydraulic_sys.set_ground(0)
hydraulic_sys.build_system()

S1, S2, rhs = hydraulic_sys.get_system()
sol = spsolve(S1.tocsr(), rhs)
solution = hydraulic_sys.build_intensity_and_voltage_from_vector(sol)

pressure_input = 10000
pressure_node = 0
potentials = solution.potentials - solution.potentials[pressure_node] + pressure_input

print("Pressures in the system:", potentials)
print("Debit through the system", solution.intensities_sources[0])
```
