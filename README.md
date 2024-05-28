# JuliaCon 2024: Towards Aerodynamic Simulations in Julia with Trixi.jl

## Abstract

In this presentation, we discuss recent developments in Trixi.jl, a numerical simulation framework designed for the accurate computation of the physics of flowing compressible media.

We focus on the technical aspect of importing existing meshes into Trixi.jl, exemplified by the utilization of the NACA4412 airfoil mesh from NASA's Turbulence Model Benchmarking Working Group.
This functionality provides users with the ability to integrate existing meshes into their simulations.

## Description

Mesh generation is a complex and critical process at the beginning of every grid-based simulation workflow.

We demonstrate how Trixi.jl facilitates this task by showcasing the usage of existing meshes, exemplified by the NACA4412 airfoil meshes from NASA's Turbulence Model Benchmarking Working Group.
This functionality represents a substantial leap forward, allowing users to seamlessly integrate well-established meshes into their simulations, a feature particularly significant given the intricacies and resource-intensive nature of mesh generation.
In more detail, we discuss the import of bilinear meshes into P4est.jl, a Julia wrapper around the p4est library which provides the datastructures for distributed memory adaptive meshes.

As we delve into the technicalities of Trixi.jl, we shine a spotlight on its pragmatic application in real-world aerodynamic scenarios. The presentation underscores the straightforward integration of advanced features, notably adaptive mesh refinement (AMR), automatic differentiation through ForwardDiff.jl and shock capturing, showcasing their utility in enhancing the precision, robustness and efficiency of simulations.
