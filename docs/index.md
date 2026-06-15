# ElecSolver

## Overview

**ElecSolver** formalizes electric systems as linear problems, suitable for both **temporal** and **frequency-domain** studies.
It focuses on constructing the linear system representation, leaving the actual numerical resolution to the user.

This repository is **not** a general-purpose electrical system solver. Instead, it acts as a **bridge** between:

- The graph-based description of an electric network
- The corresponding sparse linear system to solve

Its main goal is to provide a friendly Python interface for simulating analog electric systems. While suitable for small circuit simulations, its strength lies in its scalability: it is able to build linear systems with millions of nodes and components.

!!! warning

    ElecSolver has been designed with the following specifications in mind:

    - The time needed for building the linear system must be negligible compared to the time needed for solving it.
    - Handle natively inductive mutuals and resistive mutuals.
    - Handle as many coupled electric systems as needed.
    - Deal with lonely nodes and lonely edges in the electric graph when the problem is still well posed.

!!! note

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

The documentation is organized around the same component split as the README:

- [FrequencySystemBuilder](components/frequency-system-builder.md)
- [TemporalSystemBuilder](components/temporal-system-builder.md)

Then follow the same supporting sections:

- [Solver suggestions](solver-suggestions.md)
- [Extra uses: Hydraulic or Thermal system modeling](extra-uses-hydraulic-or-thermal-system-modeling.md)
- [Netlist import feature](netlist-import-feature.md)
