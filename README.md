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



> [!IMPORTANT]
> ElecSolver has been designed with the following specifications in mind:
> - The time needed for building the linear system must be negligible compared to the time needed for solving it.
> - Handle natively inductive mutuals and resistive mutuals (skin effects).
> - Handle as many coupled electric systems as one wants.
> - Deal with lonely nodes and lonely edges in the electric graph: the problem can be well-posed and thus solved.
> - Allow backpropagation of gradients through the system for optimization purposes.



> [!NOTE]
> Non-linear components are not supported. You must manage event detection and system updates yourself.


## Table of contents

- [ElecSolver](#elecsolver)
  - [Overview](#overview)
  - [Table of contents](#table-of-contents)
  - [How to install](#how-to-install)
    - [Using pip](#using-pip)
    - [Using conda/mamba](#using-condamamba)
    - [From source](#from-source)
  - [Documentation](#documentation)



## How to install
### Using pip
ElecSolver is distributed on PyPI and can be installed using pip or conda/mamba.
```
pip install ElecSolver
```
### Using conda/mamba
ElecSolver is also available on conda-forge and can be installed using conda or mamba:
```
conda install elecsolver
```
For solving linear systems, we advise using MUMPS through `python-mumps`, available for Linux, macOS, and Windows. It can be installed via conda:
```
conda install python-mumps
```
### From source
You can also install ElecSolver from source by cloning the repository and running:
```
pip install .
```
you will need to have `numpy`, `scipy` and `networkx` installed in your environment

## Documentation

The documentation is available in the `docs` folder and can be built using `zensical`. You can also view the documentation online at [https://williampiat3.github.io/ElecSolver/](https://williampiat3.github.io/ElecSolver/).