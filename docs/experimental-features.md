# Experimental Features

## Robust init system

ElecSolver provides an init system feature that allows for users to obtain the inital condition of the electric system right when the simulation starts.

This feature is essential when assessing the capacity of the system's reponse to a threshold input signal.

Some circuits can't answer to a theshold input:

- A voltage source with capacitance connected to its terminals.
- A current source whose terminals need inductances to be in a connected graph

In these cases, the init system is not defined and the system can only answer to differentiable inputs.

When these conditions are not met the system can answer to a threshold input.

When it is working the init system provides the solution at time $t=0+$. However, this feature does not always work:
- It fails when having two capacitances in parallel in the system
- It fails when having an inductance that is connected to a dead end in the circuit.

To circonvent these limitation we provided a workaround by providing a robust init system that can handle these special cases more effectively. Should your init system fail under the standard conditions, the robust init system can be used as an alternative.

The robust init system does not provide a physical solution at $t=0+$ in the same sense as the standard init system. Instead, it offers a mathematically consistent approximation that allows the simulation to proceed even in cases where the standard init system would fail.

It is important to note that while the robust init system increases the range of circuits that can be simulated, a user should not rely on its solution as a substitute for the physically accurate initial conditions provided by the standard init system.

We are still trying to unify both solutions into a single coherent framework that can handle a wider range of initial conditions while maintaining physical accuracy.


## Example Usage
```python
## Simple tetrahedron
res_coords  = np.array([[1],[2]],dtype=int)
res_data = np.array([1],dtype=float)
## previous init system could not handle coils connecting a single node
coil_coords  = np.array([[],[]],dtype=int)
coil_data = np.array([],dtype=float)
## previous init system could not handle 2 capacitors connected to the same node
capa_coords = np.array([[0],[1]],dtype=int)
capa_data = np.array([1],dtype=float)

## total impedance
mutual_coords=np.array([[],[]],dtype=int)
mutual_data = np.array([],dtype=float)


res_mutual_coords=np.array([[],[]],dtype=int)
res_mutual_data = np.array([],dtype=float)

elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutual_coords,mutual_data,res_mutual_coords,res_mutual_data)
elec_sys.add_voltage_source(10,0,2)
elec_sys.set_ground(0)
elec_sys.build_system()
S_i_robust,b = elec_sys.get_robust_init_system()
sol = spsolve(S_i_robust.tocsr(),b)
S_i, rhs = elec_sys.get_init_system()
sol_init = spsolve(S_i.tocsr(),rhs)
print("The two init system solutions are different")
print("Robust init system solution:",sol)
print("Standard init system solution:",sol_init)
S1,S2,b = elec_sys.get_system()
dt = 0.1
print("Solving first timestep of the system with both init system solutions provide the same solution")
np.testing.assert_allclose(spsolve(S2+dt*S1,b*dt+S2@sol), spsolve(S2+dt*S1,b*dt+S2@sol_init))
```
