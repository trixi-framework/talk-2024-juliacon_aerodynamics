# JuliaCon 2024: Towards Aerodynamic Simulations in Julia with Trixi.jl

## Talk Information

### Abstract

In this presentation, we discuss recent developments in Trixi.jl, a numerical simulation framework designed for the accurate computation of the physics of flowing compressible media.

We focus on the technical aspect of importing existing meshes into Trixi.jl, exemplified by the utilization of the NACA4412 airfoil mesh from NASA's Turbulence Model Benchmarking Working Group.
This functionality provides users with the ability to integrate existing meshes into their simulations.

### Description

Mesh generation is a complex and critical process at the beginning of every grid-based simulation workflow.

We demonstrate how Trixi.jl facilitates this task by showcasing the usage of existing meshes, exemplified by the NACA4412 airfoil meshes from NASA's Turbulence Model Benchmarking Working Group.
This functionality represents a substantial leap forward, allowing users to seamlessly integrate well-established meshes into their simulations, a feature particularly significant given the intricacies and resource-intensive nature of mesh generation.
In more detail, we discuss the import of bilinear meshes into P4est.jl, a Julia wrapper around the p4est library which provides the datastructures for distributed memory adaptive meshes.

As we delve into the technicalities of Trixi.jl, we shine a spotlight on its pragmatic application in real-world aerodynamic scenarios. The presentation underscores the straightforward integration of advanced features, notably adaptive mesh refinement (AMR), automatic differentiation through ForwardDiff.jl and shock capturing, showcasing their utility in enhancing the precision, robustness and efficiency of simulations.

## About this Reproducibility Repository

This repository contains all necessary files to re-run the simulations involving the NACA4412 airfoil shown in the presentation.

### Meshes

The mesh-data (as downloaded from the [NASA TMR 4412 Grid page](https://turbmodels.larc.nasa.gov/naca4412sep_grids.html)) as well as the processed files are available in the `mesh_data` directory.
The scripts to make the downloaded mesh `Trixi.jl` ready are present in the `pre_proc_meshes` directory, which contains also a README describing the process in some more detail.
Please note that this step involves usage of the third-party [`p3d2gmsh.py`](https://github.com/mrklein/p3d2gmsh) tool.

### Simulation Setups

The simulation setups can be found in the `elixirs` directory.
Note that for both the inviscid and viscous case we provide two files.
One is the simulation run with the "vanilla" NASA TMR grid: `Inviscid_TMR.jl, Viscous_TMR.jl`.
The others are with adaptive mesh refinement (`Inviscid_TMR+AMR.jl, Viscous_TMR+AMR.jl`) using a warm-start from the restart files from the previous simulation without AMR.
These runs with AMR are then also used to monitor the pressure coefficient along the airfoil surface.
Note that we also provide the restart files for the AMR runs in the `restart_files` directory.