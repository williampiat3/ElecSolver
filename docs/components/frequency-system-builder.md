# FrequencySystemBuilder

This class handles **frequency-domain** analysis of linear electric systems.

## Features

- Supports tension and intensity sources
- Models inductive and resistive mutuals
- Detects and couples multiple subsystems
- Accepts arbitrary complex impedances and mutuals
- Constructs sparse linear systems in COO format

!!! tip

    Some solvers do not support complex-valued systems. Use `cast_complex_system_in_real_system` from `utils.py` to convert an `n`-dimensional complex system into a `2n`-dimensional real system.

## Example

We would like to study the following system:

![Multiple system](../img/schema.png)

This can be defined in the following manner. We took `R=1`, `L=1` and `M=2`.

```python
import numpy as np
from scipy.sparse.linalg import spsolve
from ElecSolver import FrequencySystemBuilder


# Complex and sparse impedance matrix
# notice coil impedence between points 0 and 2, and coil impedence between 3 and 4
impedence_coords = np.array([[0, 0, 1, 3], [1, 2, 2, 4]], dtype=int)
impedence_data = np.array([1, 1j, 1, 1j], dtype=complex)

# Mutual inductance or coupling
# The indexes here are the impedence indexes in impedence_data
# The coupling is inductive
mutuals_coords = np.array([[1], [3]], dtype=int)
mutuals_data = np.array([2.0j], dtype=complex)

electric_sys = FrequencySystemBuilder(
    impedence_coords,
    impedence_data,
    mutuals_coords,
    mutuals_data,
)

# Add source (current source here)
electric_sys.add_current_source(intensity=10, input_node=2, output_node=0)

# Set ground
# 2 values because one for each subsystem
electric_sys.set_ground(0, 3)

# Build system
electric_sys.build_system()

# Get and solve the system
sys, b = electric_sys.get_system()
sol = spsolve(sys.tocsr(), b)
frequencial_response = electric_sys.build_intensity_and_voltage_from_vector(sol)

# We see a tension appearing on the lonely coil (between node 3 and 4)
print(frequencial_response.potentials[3] - frequencial_response.potentials[4])
```

## Adding a Parallel Resistance

We want to add components in parallel with existing components, for instance inserting a resistor in parallel with the first inductance between nodes 0 and 2.

![Parallel system](../img/schema3.png)

In Python, simply add the resistance to the list of impedances in the first lines of the script:

```python
import numpy as np
from scipy.sparse.linalg import spsolve
from ElecSolver import FrequencySystemBuilder


# We add an additional resistance between 0 and 2
impedence_coords = np.array([[0, 0, 1, 3, 0], [1, 2, 2, 4, 2]], dtype=int)
impedence_data = np.array([1, 1j, 1, 1j, 1], dtype=complex)

# No need to change the couplings since indexes of the coils did not change
mutuals_coords = np.array([[1], [3]], dtype=int)
mutuals_data = np.array([2.0j], dtype=complex)
```
