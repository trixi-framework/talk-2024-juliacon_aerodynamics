# See PR https://github.com/trixi-framework/Trixi.jl/pull/1920

using Trixi

mutable struct AnalysisCallback{Analyzer, AnalysisIntegrals, AnalysisPointwise,
                                InitialStateIntegrals, Cache}
    start_time::Float64
    start_time_last_analysis::Float64
    ncalls_rhs_last_analysis::Int
    start_gc_time::Float64
    interval::Int
    save_analysis::Bool
    output_directory::String
    analysis_filename::String
    analyzer::Analyzer
    analysis_errors::Vector{Symbol}
    analysis_integrals::AnalysisIntegrals
    analysis_pointwise::AnalysisPointwise
    initial_state_integrals::InitialStateIntegrals
    cache::Cache
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{<:Any, <:AnalysisCallback})
    @nospecialize cb # reduce precompilation time

    if get(io, :compact, false)
        show(io, cb)
    else
        analysis_callback = cb.affect!

        setup = Pair{String, Any}["interval" => analysis_callback.interval,
                                  "analyzer" => analysis_callback.analyzer]
        for (idx, error) in enumerate(analysis_callback.analysis_errors)
            push!(setup, "│ error " * string(idx) => error)
        end
        for (idx, integral) in enumerate(analysis_callback.analysis_integrals)
            push!(setup, "│ integral " * string(idx) => integral)
        end
        for (idx, quantity) in enumerate(analysis_callback.analysis_pointwise)
            push!(setup, "│ pointwise " * string(idx) => quantity)
        end
        push!(setup,
              "save analysis to file" => analysis_callback.save_analysis ? "yes" : "no")
        if analysis_callback.save_analysis
            push!(setup, "│ filename" => analysis_callback.analysis_filename)
            push!(setup,
                  "│ output directory" => abspath(normpath(analysis_callback.output_directory)))
        end
        Trixi.summary_box(io, "AnalysisCallback", setup)
    end
end

function Trixi.AnalysisCallback(mesh, equations::Trixi.AbstractEquations, solver, cache;
                                interval = 0,
                                save_analysis = false,
                                output_directory = "out",
                                analysis_filename = "analysis.dat",
                                extra_analysis_errors = Symbol[],
                                analysis_errors = union(default_analysis_errors(equations),
                                                        extra_analysis_errors),
                                extra_analysis_integrals = (),
                                analysis_integrals = union(default_analysis_integrals(equations),
                                                            extra_analysis_integrals),
                                analysis_pointwise = (),
                                RealT = real(solver),
                                uEltype = eltype(cache.elements),
                                kwargs...)
    # Decide when the callback is activated.
    # With error-based step size control, some steps can be rejected. Thus,
    #   `integrator.iter >= integrator.stats.naccept`
    #    (total #steps)       (#accepted steps)
    # We need to check the number of accepted steps since callbacks are not
    # activated after a rejected step.
    condition = (u, t, integrator) -> interval > 0 &&
        (integrator.stats.naccept % interval == 0 || Trixi.isfinished(integrator))

    analyzer = Trixi.SolutionAnalyzer(solver; kwargs...)
    cache_analysis = Trixi.create_cache_analysis(analyzer, mesh, equations, solver, cache,
                                           RealT, uEltype)

    analysis_callback = AnalysisCallback(0.0, 0.0, 0, 0.0,
                                         interval, save_analysis, output_directory,
                                         analysis_filename,
                                         analyzer,
                                         analysis_errors, Tuple(analysis_integrals),
                                         Tuple(analysis_pointwise),
                                         SVector(ntuple(_ -> zero(uEltype),
                                                        Val(nvariables(equations)))),
                                         cache_analysis)

    DiscreteCallback(condition, analysis_callback,
                     save_positions = (false, false),
                     initialize = initialize!)
end


function initialize!(cb::DiscreteCallback{Condition, Affect!}, u_ode, t,
                     integrator) where {Condition, Affect! <: AnalysisCallback}
    semi = integrator.p
    du_ode = first(get_tmp_cache(integrator))
    initialize!(cb, u_ode, du_ode, t, integrator, semi)
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u_ode, du_ode, t,
                     integrator, semi) where {Condition, Affect! <: AnalysisCallback}
    initial_state_integrals = integrate(u_ode, semi)
    _, equations, _, _ = Trixi.mesh_equations_solver_cache(semi)

    analysis_callback = cb.affect!
    analysis_callback.initial_state_integrals = initial_state_integrals
    @unpack save_analysis, output_directory, analysis_filename, analysis_errors, analysis_integrals, analysis_pointwise = analysis_callback

    if save_analysis && Trixi.mpi_isroot()
        mkpath(output_directory)

        # write header of output file
        open(joinpath(output_directory, analysis_filename), "w") do io
            Trixi.@printf(io, "#%-8s", "timestep")
            Trixi.@printf(io, "  %-14s", "time")
            Trixi.@printf(io, "  %-14s", "dt")
            if :l2_error in analysis_errors
                for v in varnames(cons2cons, equations)
                    Trixi.@printf(io, "   %-14s", "l2_"*v)
                end
            end
            if :linf_error in analysis_errors
                for v in varnames(cons2cons, equations)
                    Trixi.@printf(io, "   %-14s", "linf_"*v)
                end
            end
            if :conservation_error in analysis_errors
                for v in varnames(cons2cons, equations)
                    Trixi.@printf(io, "   %-14s", "cons_"*v)
                end
            end
            if :residual in analysis_errors
                for v in varnames(cons2cons, equations)
                    Trixi.@printf(io, "   %-14s", "res_"*v)
                end
            end
            if :l2_error_primitive in analysis_errors
                for v in varnames(cons2prim, equations)
                    Trixi.@printf(io, "   %-14s", "l2_"*v)
                end
            end
            if :linf_error_primitive in analysis_errors
                for v in varnames(cons2prim, equations)
                    Trixi.@printf(io, "   %-14s", "linf_"*v)
                end
            end

            for quantity in analysis_integrals
                Trixi.@printf(io, "   %-14s", pretty_form_ascii(quantity))
            end

            for quantity in analysis_pointwise
                Trixi.@printf(io, "   %-14s", pretty_form_ascii(quantity))
            end

            println(io)
        end
    end

    # Record current time using a high-resolution clock
    analysis_callback.start_time = time_ns()

    # Record current time for performance index computation
    analysis_callback.start_time_last_analysis = time_ns()

    # Record current number of `rhs!` calls for performance index computation
    analysis_callback.ncalls_rhs_last_analysis = Trixi.ncalls(semi.performance_counter)

    # Record total time spent in garbage collection so far using a high-resolution clock
    # Note: For details see the actual callback function below
    analysis_callback.start_gc_time = Base.gc_time_ns()

    analysis_callback(u_ode, du_ode, integrator, semi)
    return nothing
end


function (analysis_callback::AnalysisCallback)(u_ode, du_ode, integrator, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    @unpack dt, t = integrator
    iter = integrator.stats.naccept

    # Compute the percentage of the simulation that is done
    t = integrator.t
    t_initial = first(integrator.sol.prob.tspan)
    t_final = last(integrator.sol.prob.tspan)
    sim_time_percentage = (t - t_initial) / (t_final - t_initial) * 100

    # Record performance measurements and compute performance index (PID)
    runtime_since_last_analysis = 1.0e-9 * (time_ns() -
                                   analysis_callback.start_time_last_analysis)
    # PID is an MPI-aware measure of how much time per global degree of freedom (i.e., over all ranks)
    # and per `rhs!` evaluation is required. MPI-aware means that it essentially adds up the time
    # spent on each MPI rank. Thus, in an ideally parallelized program, the PID should be constant
    # independent of the number of MPI ranks used, since, e.g., using 4x the number of ranks should
    # divide the runtime on each rank by 4. See also the Trixi.jl docs ("Performance" section) for
    # more information.
    Trixi.ncalls_rhs_since_last_analysis = (Trixi.ncalls(semi.performance_counter)
                                      -
                                      analysis_callback.ncalls_rhs_last_analysis)
    performance_index = runtime_since_last_analysis * Trixi.mpi_nranks() /
                        (Trixi.Trixi.ndofsglobal(mesh, solver, cache)
                         *
                         Trixi.ncalls_rhs_since_last_analysis)

    # Compute the total runtime since the analysis callback has been initialized, in seconds
    runtime_absolute = 1.0e-9 * (time_ns() - analysis_callback.start_time)

    # Compute the relative runtime as time spent in `rhs!` divided by the number of calls to `rhs!`
    # and the number of local degrees of freedom
    # OBS! This computation must happen *after* the PID computation above, since `take!(...)`
    #      will reset the number of calls to `rhs!`
    runtime_relative = 1.0e-9 * take!(semi.performance_counter) / Trixi.ndofs(semi)

    # Compute the total time spent in garbage collection since the analysis callback has been
    # initialized, in seconds
    # Note: `Base.gc_time_ns()` is not part of the public Julia API but has been available at least
    #        since Julia 1.6. Should this function be removed without replacement in a future Julia
    #        release, just delete this analysis quantity from the callback.
    # Source: https://github.com/JuliaLang/julia/blob/b540315cb4bd91e6f3a3e4ab8129a58556947628/base/timing.jl#L83-L84
    gc_time_absolute = 1.0e-9 * (Base.gc_time_ns() - analysis_callback.start_gc_time)

    # Compute the percentage of total time that was spent in garbage collection
    gc_time_percentage = gc_time_absolute / runtime_absolute * 100

    # Obtain the current memory usage of the Julia garbage collector, in MiB, i.e., the total size of
    # objects in memory that have been allocated by the JIT compiler or the user code.
    # Note: `Base.gc_live_bytes()` is not part of the public Julia API but has been available at least
    #        since Julia 1.6. Should this function be removed without replacement in a future Julia
    #        release, just delete this analysis quantity from the callback.
    # Source: https://github.com/JuliaLang/julia/blob/b540315cb4bd91e6f3a3e4ab8129a58556947628/base/timing.jl#L86-L97
    memory_use = Base.gc_live_bytes() / 2^20 # bytes -> MiB

    Trixi.@trixi_timeit Trixi.timer() "analyze solution" begin
        # General information
        Trixi.mpi_println()
        Trixi.mpi_println("─"^100)
        Trixi.mpi_println(" Simulation running '", Trixi.get_name(equations), "' with ",
                    summary(solver))
        Trixi.mpi_println("─"^100)
        Trixi.mpi_println(" #timesteps:     " * Trixi.@sprintf("% 14d", iter) *
                    "               " *
                    " run time:       " * Trixi.@sprintf("%10.8e s", runtime_absolute))
        Trixi.mpi_println(" Δt:             " * Trixi.@sprintf("%10.8e", dt) *
                    "               " *
                    " └── GC time:    " *
                    Trixi.@sprintf("%10.8e s (%5.3f%%)", gc_time_absolute, gc_time_percentage))
        Trixi.mpi_println(rpad(" sim. time:      " *
                         Trixi.@sprintf("%10.8e (%5.3f%%)", t, sim_time_percentage), 46) *
                    " time/DOF/rhs!:  " * Trixi.@sprintf("%10.8e s", runtime_relative))
        Trixi.mpi_println("                 " * "              " *
                    "               " *
                    " PID:            " * Trixi.@sprintf("%10.8e s", performance_index))
        Trixi.mpi_println(" #DOFs per field:" * Trixi.@sprintf("% 14d", Trixi.ndofsglobal(semi)) *
                    "               " *
                    " alloc'd memory: " * Trixi.@sprintf("%14.3f MiB", memory_use))
        Trixi.mpi_println(" #elements:      " *
                    Trixi.@sprintf("% 14d", Trixi.nelementsglobal(mesh, solver, cache)))

        # Level information (only show for AMR)
        Trixi.print_amr_information(integrator.opts.callback, mesh, solver, cache)
        Trixi.mpi_println()

        # Open file for appending and store time step and time information
        if Trixi.mpi_isroot() && analysis_callback.save_analysis
            io = open(joinpath(analysis_callback.output_directory,
                               analysis_callback.analysis_filename), "a")
            Trixi.@printf(io, "% 9d", iter)
            Trixi.@printf(io, "  %10.8e", t)
            Trixi.@printf(io, "  %10.8e", dt)
        else
            io = devnull
        end

        # Calculate current time derivative (needed for semidiscrete entropy time derivative, residual, etc.)
        # `integrator.f` is usually just a call to `rhs!`
        # However, we want to allow users to modify the ODE RHS outside of Trixi.jl
        # and allow us to pass a combined ODE RHS to OrdinaryDiffEq, e.g., for
        # hyperbolic-parabolic systems.
        Trixi.@notimeit Trixi.timer() integrator.f(du_ode, u_ode, semi, t)
        u = Trixi.wrap_array(u_ode, mesh, equations, solver, cache)
        du = Trixi.wrap_array(du_ode, mesh, equations, solver, cache)
        # Compute l2_error, linf_error among other quantities
        analysis_callback(io, du, u, u_ode, t, semi, iter)

        Trixi.mpi_println("─"^100)
        Trixi.mpi_println()

        # Add line break and close analysis file if it was opened
        if Trixi.mpi_isroot() && analysis_callback.save_analysis
            # This resolves a possible type instability introduced above, since `io`
            # can either be an `IOStream` or `devnull`, but we know that it must be
            # an `IOStream here`.
            println(io::IOStream)
            close(io::IOStream)
        end
    end

    # avoid re-evaluating possible FSAL stages
    u_modified!(integrator, false)

    # Reset performance measurements
    analysis_callback.start_time_last_analysis = time_ns()
    analysis_callback.ncalls_rhs_last_analysis = Trixi.ncalls(semi.performance_counter)

    return nothing
end

function (analysis_callback::AnalysisCallback)(io, du, u, u_ode, t, semi, iter)
    @unpack analyzer, analysis_errors, analysis_integrals, analysis_pointwise = analysis_callback
    cache_analysis = analysis_callback.cache
    _, equations, _, _ = Trixi.mesh_equations_solver_cache(semi)

    # Calculate and print derived quantities (error norms, entropy etc.)
    # Variable names required for L2 error, Linf error, and conservation error
    if any(q in analysis_errors
           for q in (:l2_error, :linf_error, :conservation_error, :residual)) &&
       Trixi.mpi_isroot()
        print(" Variable:    ")
        for v in eachvariable(equations)
            Trixi.@printf("   %-14s", varnames(cons2cons, equations)[v])
        end
        println()
    end

    if :l2_error in analysis_errors || :linf_error in analysis_errors
        # Calculate L2/Linf errors
        l2_error, linf_error = calc_error_norms(u_ode, t, analyzer, semi,
                                                cache_analysis)

        if Trixi.mpi_isroot()
            # L2 error
            if :l2_error in analysis_errors
                print(" L2 error:    ")
                for v in eachvariable(equations)
                    Trixi.@printf("  % 10.8e", l2_error[v])
                    Trixi.@printf(io, "  % 10.8e", l2_error[v])
                end
                println()
            end

            # Linf error
            if :linf_error in analysis_errors
                print(" Linf error:  ")
                for v in eachvariable(equations)
                    Trixi.@printf("  % 10.8e", linf_error[v])
                    Trixi.@printf(io, "  % 10.8e", linf_error[v])
                end
                println()
            end
        end
    end

    # Conservation error
    if :conservation_error in analysis_errors
        @unpack initial_state_integrals = analysis_callback
        state_integrals = integrate(u_ode, semi)

        if Trixi.mpi_isroot()
            print(" |∑U - ∑U₀|:  ")
            for v in eachvariable(equations)
                err = abs(state_integrals[v] - initial_state_integrals[v])
                Trixi.@printf("  % 10.8e", err)
                Trixi.@printf(io, "  % 10.8e", err)
            end
            println()
        end
    end

    # Residual (defined here as the vector maximum of the absolute values of the time derivatives)
    if :residual in analysis_errors
        mpi_print(" max(|Uₜ|):   ")
        for v in eachvariable(equations)
            # Calculate maximum absolute value of Uₜ
            res = maximum(abs, view(du, v, ..))
            if mpi_isparallel()
                # TODO: Debugging, here is a type instability
                global_res = MPI.Reduce!(Ref(res), max, mpi_root(), mpi_comm())
                if Trixi.mpi_isroot()
                    res::eltype(du) = global_res[]
                end
            end
            if Trixi.mpi_isroot()
                Trixi.@printf("  % 10.8e", res)
                Trixi.@printf(io, "  % 10.8e", res)
            end
        end
        Trixi.mpi_println()
    end

    # L2/L∞ errors of the primitive variables
    if :l2_error_primitive in analysis_errors ||
       :linf_error_primitive in analysis_errors
        l2_error_prim, linf_error_prim = calc_error_norms(cons2prim, u_ode, t, analyzer,
                                                          semi, cache_analysis)

        if Trixi.mpi_isroot()
            print(" Variable:    ")
            for v in eachvariable(equations)
                Trixi.@printf("   %-14s", varnames(cons2prim, equations)[v])
            end
            println()

            # L2 error
            if :l2_error_primitive in analysis_errors
                print(" L2 error prim.: ")
                for v in eachvariable(equations)
                    Trixi.@printf("%10.8e   ", l2_error_prim[v])
                    Trixi.@printf(io, "  % 10.8e", l2_error_prim[v])
                end
                println()
            end

            # L∞ error
            if :linf_error_primitive in analysis_errors
                print(" Linf error pri.:")
                for v in eachvariable(equations)
                    Trixi.@printf("%10.8e   ", linf_error_prim[v])
                    Trixi.@printf(io, "  % 10.8e", linf_error_prim[v])
                end
                println()
            end
        end
    end

    # additional integrals
    Trixi.analyze_integrals(analysis_integrals, io, du, u, t, semi)
    # additional pointwise quantities
    analyze_pointwise(analysis_pointwise, du, u, t, semi, iter)

    return nothing
end

# Iterate over tuples of pointwise analysis quantities in a type-stable way using "lispy tuple programming".
function analyze_pointwise(analysis_quantities::NTuple{N, Any}, du, u, t,
                           semi, iter) where {N}

    # Extract the first pointwise analysis quantity and process it; keep the remaining to be processed later
    quantity = first(analysis_quantities)
    remaining_quantities = Base.tail(analysis_quantities)

    Trixi.analyze(quantity, du, u, t, semi, iter)

    # Recursively call this method with the unprocessed pointwise analysis quantities
    analyze_pointwise(remaining_quantities, du, u, t, semi, iter)
    return nothing
end

# terminate the type-stable iteration over tuples
function analyze_pointwise(analysis_quantities::Tuple{}, du, u, t, semi, iter)
    nothing
end

function Trixi.analyze(quantity, du, u, t, semi::Trixi.AbstractSemidiscretization, iter)
                 mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
                 Trixi.analyze(quantity, du, u, t, mesh, equations, solver, cache, iter)
end

function Trixi.analyze(quantity::AnalysisSurfaceIntegral{Variable}, du, u, t,
                 semi::Trixi.AbstractSemidiscretization) where {Variable}
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    Trixi.analyze(quantity, du, u, t, mesh, equations, solver, cache, semi)
end

"""
    AnalysisSurfacePointwise{Variable, NBoundaries}(boundary_symbol_or_boundary_symbols,
                                                    variable, output_directory = "out")

This struct is used to compute pointwise surface values of a quantity of interest `variable`
alongside the boundary/boundaries associated with particular name(s) given in `boundary_symbol`
or `boundary_symbols`.
For instance, this can be used to compute the surface pressure coefficient [`SurfacePressureCoefficient`](@ref) or
surface friction coefficient [`SurfaceFrictionCoefficient`](@ref) of e.g. an airfoil with the boundary
symbol `:Airfoil` in 2D.

- `boundary_symbols::NTuple{NBoundaries, Symbol}`: Name(s) of the boundary/boundaries
  where the quantity of interest is computed
- `variable::Variable`: Quantity of interest, like lift or drag
- `output_directory = "out"`: Directory where the pointwise value files are stored.
"""
struct AnalysisSurfacePointwise{Variable, NBoundaries}
    variable::Variable # Quantity of interest, like lift or drag
    boundary_symbols::NTuple{NBoundaries, Symbol} # Name(s) of the boundary/boundaries
    output_directory::String

    function AnalysisSurfacePointwise(boundary_symbols::NTuple{NBoundaries, Symbol},
                                      variable,
                                      output_directory = "out") where {NBoundaries}
        return new{typeof(variable), NBoundaries}(variable, boundary_symbols,
                                                  output_directory)
    end
end

struct FlowState{RealT <: Real}
    rhoinf::RealT
    uinf::RealT
    linf::RealT
end

struct SurfacePressureCoefficient{RealT <: Real}
    pinf::RealT # Free stream pressure
    flow_state::FlowState{RealT}
end

struct SurfaceFrictionCoefficient{RealT <: Real} <: Trixi.VariableViscous
    flow_state::FlowState{RealT}
end

"""
    SurfacePressureCoefficient(pinf, rhoinf, uinf, linf)

Compute the surface pressure coefficient
```math
C_p \\coloneqq \\frac{p - p_{p_\\infty}}
                     {0.5 \\rho_{\\infty} U_{\\infty}^2 L_{\\infty}}
```
based on the pressure distribution along a boundary.
Supposed to be used in conjunction with [`AnalysisSurfacePointwise`](@ref)
which stores the boundary information and semidiscretization.

- `pinf::Real`: Free-stream pressure
- `rhoinf::Real`: Free-stream density
- `uinf::Real`: Free-stream velocity
- `linf::Real`: Reference length of geometry (e.g. airfoil chord length)
"""
function SurfacePressureCoefficient(pinf, rhoinf, uinf, linf)
    return SurfacePressureCoefficient(pinf, FlowState(rhoinf, uinf, linf))
end

"""
SurfaceFrictionCoefficient(rhoinf, uinf, linf)

Compute the surface skin friction coefficient
```math
C_f \\coloneqq \\frac{\\boldsymbol \\tau_w  \\boldsymbol n^\\perp}
                     {0.5 \\rho_{\\infty} U_{\\infty}^2 L_{\\infty}}
```
based on the wall shear stress vector ``\\tau_w`` along a boundary.
Supposed to be used in conjunction with [`AnalysisSurfacePointwise`](@ref)
which stores the boundary information and semidiscretization.

- `rhoinf::Real`: Free-stream density
- `uinf::Real`: Free-stream velocity
- `linf::Real`: Reference length of geometry (e.g. airfoil chord length)
"""
function SurfaceFrictionCoefficient(rhoinf, uinf, linf)
    return SurfaceFrictionCoefficient(FlowState(rhoinf, uinf, linf))
end

function (pressure_coefficient::SurfacePressureCoefficient)(u, equations)
    p = pressure(u, equations)
    @unpack pinf = pressure_coefficient
    @unpack rhoinf, uinf, linf = pressure_coefficient.flow_state
    return (p - pinf) / (0.5 * rhoinf * uinf^2 * linf)
end

function (surface_friction::SurfaceFrictionCoefficient)(u, normal_direction, x, t,
                                                        equations_parabolic,
                                                        gradients_1, gradients_2)
    viscous_stress_vector_ = viscous_stress_vector(u, normal_direction,
                                                   equations_parabolic,
                                                   gradients_1, gradients_2)
    @unpack rhoinf, uinf, linf = surface_friction.flow_state

    # Normalize as `normal_direction` is not necessarily a unit vector
    n = normal_direction / norm(normal_direction)
    n_perp = (-n[2], n[1])
    return (viscous_stress_vector_[1] * n_perp[1] +
            viscous_stress_vector_[2] * n_perp[2]) /
           (0.5 * rhoinf * uinf^2 * linf)
end

function Trixi.analyze(surface_variable::AnalysisSurfacePointwise, du, u, t,
                 mesh::P4estMesh{2},
                 equations, dg::DGSEM, cache, semi, iter)
    @unpack boundaries = cache
    @unpack surface_flux_values, node_coordinates, contravariant_vectors = cache.elements
    @unpack weights = dg.basis

    @unpack variable, boundary_symbols = surface_variable
    @unpack boundary_symbol_indices = semi.boundary_conditions
    indices = Trixi.get_boundary_indices(boundary_symbols, boundary_symbol_indices)

    dim = ndims(mesh)
    n_nodes = nnodes(dg)
    n_elements = length(indices)

    coordinates = Matrix{real(dg)}(undef, n_elements * n_nodes, dim) # physical coordinates of indices
    values = Vector{real(dg)}(undef, n_elements * n_nodes) # variable values at indices

    index_range = eachnode(dg)

    global_node_index = 1 # Keeps track of solution point number on the surface
    for boundary in indices
        element = boundaries.neighbor_ids[boundary]
        node_indices = boundaries.node_indices[boundary]

        i_node_start, i_node_step = Trixi.index_to_start_step_2d(node_indices[1], index_range)
        j_node_start, j_node_step = Trixi.index_to_start_step_2d(node_indices[2], index_range)

        i_node = i_node_start
        j_node = j_node_start
        for node_index in index_range
            u_node = Trixi.get_node_vars(cache.boundaries.u, equations, dg, node_index,
                                         boundary)

            x = Trixi.get_node_coords(node_coordinates, equations, dg, i_node, j_node,
                                element)
            value = variable(u_node, equations)

            coordinates[global_node_index, 1] = x[1]
            coordinates[global_node_index, 2] = x[2]
            values[global_node_index] = value

            i_node += i_node_step
            j_node += j_node_step
            global_node_index += 1
        end
    end

    save_pointwise_file(surface_variable.output_directory, varname(variable),
                        coordinates,
                        values, t, iter)
end

function Trixi.analyze(surface_variable::AnalysisSurfacePointwise{Variable},
                 du, u, t, mesh::P4estMesh{2},
                 equations, equations_parabolic,
                 dg::DGSEM, cache, semi,
                 cache_parabolic, iter) where {Variable <: Trixi.VariableViscous}
    @unpack boundaries = cache
    @unpack surface_flux_values, node_coordinates, contravariant_vectors = cache.elements
    @unpack weights = dg.basis

    @unpack variable, boundary_symbols = surface_variable
    @unpack boundary_symbol_indices = semi.boundary_conditions
    indices = Trixi.get_boundary_indices(boundary_symbols, boundary_symbol_indices)

    dim = ndims(mesh)
    n_nodes = nnodes(dg)
    n_elements = length(indices)

    coordinates = Matrix{real(dg)}(undef, n_elements * n_nodes, dim) # physical coordinates of indices
    values = Vector{real(dg)}(undef, n_elements * n_nodes) # variable values at indices

    # Additions for parabolic
    @unpack viscous_container = cache_parabolic
    @unpack gradients = viscous_container

    gradients_x, gradients_y = gradients

    index_range = eachnode(dg)
    global_node_index = 1 # Keeps track of solution point number on the surface
    for boundary in indices
        element = boundaries.neighbor_ids[boundary]
        node_indices = boundaries.node_indices[boundary]
        direction = indices2direction(node_indices)

        i_node_start, i_node_step = Trixi.index_to_start_step_2d(node_indices[1], index_range)
        j_node_start, j_node_step = Trixi.index_to_start_step_2d(node_indices[2], index_range)

        i_node = i_node_start
        j_node = j_node_start
        for node_index in index_range
            u_node = Trixi.get_node_vars(cache.boundaries.u, equations, dg, node_index,
                                         boundary)

            x = Trixi.get_node_coords(node_coordinates, equations, dg, i_node, j_node,
                                element)
            # Extract normal direction at nodes which points from the elements outwards,
            # i.e., *into* the structure.
            normal_direction = get_normal_direction(direction, contravariant_vectors,
                                                    i_node, j_node,
                                                    element)

            gradients_1 = get_node_vars(gradients_x, equations_parabolic, dg, i_node,
                                        j_node, element)
            gradients_2 = get_node_vars(gradients_y, equations_parabolic, dg, i_node,
                                        j_node, element)

            # Integral over whole boundary surface
            value = variable(u_node, normal_direction, x, t, equations_parabolic,
                             gradients_1, gradients_2)

            coordinates[global_node_index, 1] = x[1]
            coordinates[global_node_index, 2] = x[2]
            values[global_node_index] = value

            i_node += i_node_step
            j_node += j_node_step
            global_node_index += 1
        end
    end

    save_pointwise_file(surface_variable.output_directory, varname(variable),
                        coordinates,
                        values, t, iter)
end

varname(::Any) = @assert false "Surface variable name not assigned" # This makes sure default behaviour is not overwriting
varname(pressure_coefficient::SurfacePressureCoefficient) = "CP_x"
varname(friction_coefficient::SurfaceFrictionCoefficient) = "CF_x"

function save_pointwise_file(output_directory, varname, coordinates, values, t, iter)
    n_points = length(values)

    filename = joinpath(output_directory, varname) * Trixi.@sprintf("_%06d.h5", iter)

    Trixi.h5open(filename, "w") do file
        # Add context information as Trixi.attributes
        Trixi.attributes(file)["n_points"] = n_points
        Trixi.attributes(file)["variable_name"] = varname

        file["time"] = t
        file["timestep"] = iter
        file["point_coordinates"] = coordinates
        file["point_data"] = values
    end
end

function pretty_form_ascii(::AnalysisSurfacePointwise{<:SurfacePressureCoefficient{<:Any}})
    "CP(x)"
end
function pretty_form_utf(::AnalysisSurfacePointwise{<:SurfacePressureCoefficient{<:Any}})
    "CP(x)"
end

function pretty_form_ascii(::AnalysisSurfacePointwise{<:SurfaceFrictionCoefficient{<:Any}})
    "CF(x)"
end
function pretty_form_utf(::AnalysisSurfacePointwise{<:SurfaceFrictionCoefficient{<:Any}})
    "CF(x)"
end

function Trixi.analyze(quantity::AnalysisSurfacePointwise{Variable},
                 du, u, t,
                 semi::Trixi.AbstractSemidiscretization,
                 iter) where {Variable}
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    Trixi.analyze(quantity, du, u, t, mesh, equations, solver, cache, semi, iter)
end
function Trixi.analyze(quantity::AnalysisSurfacePointwise{Variable},
                 du, u, t,
                 semi::SemidiscretizationHyperbolicParabolic,
                 iter) where {
                              Variable <:
                              Trixi.VariableViscous}
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    equations_parabolic = semi.equations_parabolic
    cache_parabolic = semi.cache_parabolic
    analyze(quantity, du, u, t, mesh, equations, equations_parabolic, solver, cache, semi,
            cache_parabolic, iter)
end