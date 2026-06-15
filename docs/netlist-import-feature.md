# Netlist import feature

`NetlistParser` can import passive netlists and build a `TemporalSystemBuilder` instance. Solving the system then follows the same workflow as the other temporal examples.

```python
"""
*test netlist for python solver square.net
Iin 0 1 PWL(0 0 0.000000001 10)
L0 1 3 0.1
L1 2 0 0.1
R1 1 0 1
R2 2 3 1
c2 1 2 2
c3 0 3 2
.tran 0 4 0 0.08
.end
"""

import matplotlib.pyplot as plt
from numpy import arange
from scipy.sparse.linalg import spsolve
from ElecSolver import NetlistParser

parser = NetlistParser("square.net")
parser.map_netlist()

node_zero = parser.node_map["0"]
node_one = parser.node_map["1"]

elec_sys = parser.generate_temporal_system()
elec_sys.add_current_source(10, node_one, node_zero)
elec_sys.set_ground(node_zero)
elec_sys.build_system()

S_i, b = elec_sys.get_init_system()
sol = spsolve(S_i.tocsr(), b)
S1, S2, rhs = elec_sys.get_system()

dt = 0.08
steps = 50
vals_res1 = []
vals_res2 = []
vals_L1 = []
voltage_src = []

R1_index = list(parser.resistors.keys()).index("R1")
R2_index = list(parser.resistors.keys()).index("R2")
L1_index = list(parser.inductors.keys()).index("L1")

for i in range(steps):
    temporal_response = elec_sys.build_intensity_and_voltage_from_vector(sol)
    vals_res1.append(temporal_response.intensities_res[R1_index])
    vals_res2.append(temporal_response.intensities_res[R2_index])
    vals_L1.append(temporal_response.intensities_coil[L1_index])
    voltage_src.append(
        temporal_response.potentials[node_one] - temporal_response.potentials[node_zero]
    )
    sol = spsolve(S2 + dt * S1, b * dt + S2 @ sol)

plt.xlabel("Time")
plt.ylabel("Intensity")
plt.plot(arange(steps, dtype=float) * dt, vals_res1, label="intensity res 1")
plt.plot(arange(steps, dtype=float) * dt, vals_res2, label="intensity res 2")
plt.plot(arange(steps, dtype=float) * dt, vals_L1, label="intensity L1")
plt.plot(arange(steps, dtype=float) * dt, voltage_src, label="V(1-0)")

plt.legend()
plt.savefig("intensities_res.png")
```