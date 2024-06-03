using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

prandtl_number() = 0.72
                                                     
AoA = 0.02181661564992912 # 1.25 degreee in radians

@inline function initial_condition_mach08_flow(x, t, equations)
    # set the freestream flow parameters
    rho_freestream = 1.4

    #=
    sin_AoA, cos_AoA = sincos(0.02181661564992912)
    v = 0.8

    v1 = cos_AoA * v
    v2 = sin_AoA * v
    =#

    v1 = 0.7998096216639273
    v2 = 0.017451908027648896

    p_freestream = 1.0

    prim = SVector(rho_freestream, v1, v2, p_freestream)
    return prim2cons(prim, equations)
end

initial_condition = initial_condition_mach08_flow

# Boundary conditions for free-stream testing
boundary_condition_free_stream = BoundaryConditionDirichlet(initial_condition)

polydeg = 1 # Low polynomial degree to keep memory footprint somewhat low
surface_flux = flux_hllc

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux)

cd(@__DIR__)
path = "../mesh_data/"
mesh = "NACA4412_2_2D_unique.inp"
mesh_file = path * mesh

boundary_symbols = [:b2_symmetry_y_strong, :b4_farfield_riem, :b5_farfield_riem, 
                    :b7_farfield_riem, :b6_viscous_solid, :b8_to_stitch_a]

mesh = P4estMesh{2}(mesh_file, polydeg = polydeg, boundary_symbols = boundary_symbols)

boundary_conditions = Dict(:b2_symmetry_y_strong => boundary_condition_free_stream,
                           :b4_farfield_riem => boundary_condition_free_stream,
                           :b5_farfield_riem => boundary_condition_free_stream,
                           :b7_farfield_riem => boundary_condition_free_stream,
                           :b6_viscous_solid => boundary_condition_slip_wall,
                           :b8_to_stitch_a => boundary_condition_free_stream)

semi = SemidiscretizationHyperbolic(mesh, equations,
                                    initial_condition, solver;
                                    boundary_conditions = boundary_conditions)

A = jacobian_ad_forward(semi)

using Plots, SparseArrays

N_plot = 800
A_sparse = sparse(A[1:N_plot, 1:N_plot])

spy(A_sparse, xlabel = "\$j\$", ylabel = "\$i\$", colorbar = false, legend = nothing)

using LinearAlgebra

Eigenvalues = eigvals(A)

EigValsReal = real(Eigenvalues)
EigValsImag = imag(Eigenvalues)

scatter(EigValsReal, EigValsImag, label = "Spectrum")