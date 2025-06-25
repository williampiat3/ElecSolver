# ElectricSystemSolver
## Introduction
This repository provides a way to formalize Electric Systems as a linear problem, either for temporal studies or frequency studies.
This repository is not meant to be a general library for solving electric system, we limited our efforts to the construction of the linear system and let the user choose the solver for getting the solution. The aim is to be the bridge between the formalization of the electric system in terms of connectivity and the actual linear system to be solved

## Table of content

- [ElectricSystemSolver](#electricsystemsolver)
  - [Introduction](#introduction)
  - [Table of content](#table-of-content)
  - [FrequencySystemBuilder](#frequencysystembuilder)
    - [TemporalSystemBuilder](#temporalsystembuilder)
  - [Solver suggestions](#solver-suggestions)


This repository has two components:

1. FrequencySystemBuilder: A class for making frequency studies of linear electric systems. It takes as an input a graph of impedance in the form of a complex COO sparse matrix and another complex COO matrix for mutuals (i.e. remote interactions between impedances)
2. TemporalSystemBuilder: A class for making temporal evolutions of electric systems. It takes as an input lists and coordinates of coils, resistances and capacities along with a list of inductive and resistive mutuals.

The implementation is pure python but fully vectorized using numpy, as far as or testing goes, building a system with 3 million nodes took 10 seconds (most of the overhead taken by networkx analyzing the graph connectivity)

## FrequencySystemBuilder

The frequency system builder supports

* Tension and intensity sources
* Inductive and resistive mutual
* Auto detection and coupling of multiple subsystems
* Any complex impedence
* Any complex mutual



> [!TIP]
> Not all solvers support solving complex linear system, we therefore added utility function `cast_complex_system_in_real_system` in the `utils.py`  that casts the `n` dimensional linear system into a `2n` dimensional real linear system.

Exemple:

We would like to study the following system:
![Multiple system](img/schema.png)

this can simply be defined in the following manner (We took R=1, L=1 and M=2)
```python
import numpy as np
from scipy.sparse.linalg import spsolve
from FrequencySystemBuilder import FrequencySystemBuilder


# Sparse Python impedence matrix (notice coil impedence between points 0 and 2, and coil impedence between 3 and 4 )
impedence_coords = np.array([[0,0,1,3],[1,2,2,4]], dtype=int)
impedence_data = np.array([1, 1j, 1, 1j], dtype=complex)

# Mutual inductance or coupling
# The indexes here are the impedence indexes in impedence_data
# The coupling is inductive
mutuals_coords = np.array([[1],[3]], dtype=int)
mutuals_data = np.array([2.j], dtype=complex)

electric_sys = FrequencySystemBuilder(
    impedence_coords,
    impedence_data,
    mutuals_coords,
    mutuals_data
)

# Set node masses
electric_sys.set_mass(0, 3)
electric_sys.build_system()
electric_sys.build_second_member_intensity(intensity=10, input_node=2, output_node=0)

# Solve the system
sys, b = electric_sys.get_system()
## plotting the linear system
print(sys.todense())
print(b)

sol = spsolve(sys.tocsr(), b)
intensities, potentials = electric_sys.build_intensity_and_voltage_from_vector(sol)

## We see a tension appearing on the lonely coil (between node 3 and 4)
print(potentials[3]-potentials[4])
```

If you want to add dipoles in parallel, simply add them to the list of impedence, for instance if you want to add a resistance in parallel with the inductance to make the following scheme:
![Parallel system](img/schema3.png)
Simply change the very first lines of the script:
```python
import numpy as np
from scipy.sparse.linalg import spsolve
from FrequencySystemBuilder import FrequencySystemBuilder


# We add an additionnal resistance between 0 and 2
impedence_coords = np.array([[0,0,1,3,0],[1,2,2,4,2]], dtype=int)
impedence_data = np.array([1, 1j,1, 1j,1], dtype=complex)

# No need to change the couplings since indexes of the coils did not change
mutuals_coords = np.array([[1],[3]], dtype=int)
mutuals_data = np.array([2.j], dtype=complex)

```


### TemporalSystemBuilder

This class is meant for systems with a limited amount of features supported components

* Tension and intensity sources
* Inductive and resistive mutual
* Auto detection and coupling of multiple subsystems
* Resitances
* Capacities
* Coils
* Inductive mutuals
* Resistive mutuals

Exemple:

We would like to study the following system:
![Temporal system](img/schema2.png)

with R=1, L=0.1, C=2 this gives:
```python
import numpy as np
from scipy.sparse.linalg import spsolve
from TemporalSystemBuilder import TemporalSystemBuilder

## Defining resistances
res_coords  = np.array([[0,2],[1,3]],dtype=int)
res_data = np.array([1,1],dtype=float)
## Defining coils
coil_coords  = np.array([[1,0],[3,2]],dtype=int)
coil_data = np.array([0.1,0.1],dtype=float)
## Defining capacities
capa_coords = np.array([[1,3],[2,0]],dtype=int)
capa_data = np.array([2,2],dtype=float)

## Defining empty mutuals here
mutuals_coords=np.array([[],[]],dtype=int)
mutuals_data = np.array([],dtype=float)


res_mutuals_coords=np.array([[],[]],dtype=int)
res_mutuals_data = np.array([],dtype=float)

## initializing system
elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
## Seting mass at point 0
elec_sys.set_mass(0)
## Build second member
elec_sys.build_system()
elec_sys.build_second_member_intensity(10,1,0)
# getting initial condition system
S_i,b = elec_sys.get_init_system()
# initial condition
sol = spsolve(S_i.tocsr(),b)
# get system (S1 is real part, S2 derivative part)
S1,S2,rhs = elec_sys.get_system()
sol = spsolve(S_i,b)

## Solving using euler implicit scheme
dt=0.08
vals_res1 = []
vals_res2 = []
for i in range(50):
    currents_coil,currents_res,currents_capa,voltages,_ = elec_sys.build_intensity_and_voltage_from_vector(sol)
    vals_res1.append(currents_res[1])
    vals_res2.append(currents_res[0])
    sol = spsolve(S2+dt*S1,b*dt+S2@sol)
import matplotlib.pyplot as plt
plt.xlabel("Time")
plt.ylabel("Intensity")
plt.plot(vals_res1,label="intensity res 1")
plt.plot(vals_res2,label="intensity res 2")
plt.legend()
plt.savefig("intensities_res.png")
```

This gives the following graph that displays the intensity passing through the resistances
![Temporal system](img/intensities_res.png)


## Solver suggestions

For simple systems and use cases we advise using the sparse solver from scipy as presented in the exemples

However, for temporal resolutions we advise switching to MUMPS by using pyMUMPS as presented in the test function `tests.test_temporal_system.test_big_grid()`: MUMPS can alleviate the fact that only the second member changes through the temporal resolution making iterations fast and scalable after the analysis is performed once.