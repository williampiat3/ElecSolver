# ElectricSystemSolver

This repository provides a way to formalize Electric Systems as a linear problem, either for temporal studies or frequency studies.
This repository is not meant to be a general library for solving electric system, we limited our efforts to the construction of the linear system and let the user choose the solver for getting the solution. The aim is to be the bridge between the formalization of the electric system in terms of connectivity and the actual linear system to be solved

This repository has two components:

1. FrequencySystemBuilder: A class for making frequency studies of linear electric systems. It takes as an input a graph of impedance in the form of a complex COO sparse matrix and another complex COO matrix for mutuals (i.e. remote interactions between impedances)
2. TemporalSystemBuilder: A class for making temporal evolutions of electric systems. It takes as an input lists and coordinates of coils, resistances and capacities along with a list of inductive and resistive mutuals.

## FrequencySystemBuilder

The frequency system builder supports

* Tension and intensity sources
* Inductive and resistive mutual
* Auto detection and coupling of multiple subsystems
* Any complex impedence
* Any complex mutual



> [!TIP]
> Not all solvers support solving complex linear system, we therefore added utility function in the `utils.py`  that casts the `n` dimensional linear system into a `2n` dimensional linear system.

Exemple:

We would like to study the following system:
![Multiple system](img/schema.png)

this can simply be defined in the following manner (We took R=1, L=1 and M=2)
```python
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
print(elec_sys.build_intensity_and_voltage_from_vector(sol))
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
