# ElecSolver

## Overview

**ElecSolver** formalizes electric systems as linear problems, suitable for both **temporal** and **frequency-domain** studies.
It focuses on constructing the linear system representation, leaving the actual numerical resolution to the user.

This repository is **not** a general-purpose electrical system solver. Instead, it acts as a **bridge** between:

- The graph-based description of an electric network
- The corresponding sparse linear system to solve

In a very simple way, ElecSolver takes as an input the Resistances $R$, Capacitances $C$, Inductances $L$, Mutuals $M$, Current sources $I$  and Voltage sources $V$ along with the connectivity graph $G$ of the system and outputs the linear system: matrix $S$ and vector $b$ to solve in order to get the solution of the electric problem.

$$f_{\text{ElecSolver}}(R,C,L,M,I,V,G) = (S,b)$$

It is then up to the user to choose the solver they want for solving the system $Sx=b$ and to manage the time iterations if needed for temporal problems.

The main goal of ElecSolver is to provide a friendly Python interface for simulating and optimizing analog electric systems. While suitable for small circuit simulations, its strength lies in its scalability: it is able to build linear systems with millions of nodes and components.


!!! tip

    ElecSolver has been designed with the following specifications in mind:

    - The time needed for building the linear system must be negligible compared to the time needed for solving it.
    - Handle natively inductive mutuals and resistive mutuals.
    - Handle as many coupled electric systems as needed.
    - Deal with lonely nodes and lonely edges in the electric graph when the problem is still well posed.
    - Allow backpropagation of gradients through the system for optimization purposes.

!!! warning

    Non-linear components are not supported. You must manage event detection and system updates yourself.

## How to install

ElecSolver is distributed on PyPI and can be installed with `pip` or with `conda`/`mamba`.

```bash
pip install ElecSolver
```

or

```bash
conda install elecsolver
```

For solving the linear systems, `python-mumps` is recommended when it is available on your platform.

```bash
conda install python-mumps
```

## Components

The documentation presents the two main components of ElecSolver:

- [FrequencySystemBuilder](components/frequency-system-builder.md)
- [TemporalSystemBuilder](components/temporal-system-builder.md)

Then presents extra features:

- [Solver suggestions](solver-suggestions.md)
- [Extra uses: Hydraulic or Thermal system modeling](extra-uses-hydraulic-or-thermal-system-modeling.md)
- [Netlist import feature](netlist-import-feature.md)
