#-------------------------------------------------------------------------------------------
#=
    Purpose:    Solver interface for solving LPs
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       September 2024
=#
#-------------------------------------------------------------------------------------------

module Solve

export SolverConfig, changeSparseQFCASolver

using JuMP, HiGHS, Clp, Cbc, GLPK, ECOS, SCS

"""
    SolverConfig(name, handle)

Definition of a common solver type, which inclues the name of the solver and other parameters

- `name`:           Name of the solver
- `handle`:         Solver handle used to refer to the solver

"""

mutable struct SolverConfig
    name      ::String
    handle
end

#-------------------------------------------------------------------------------------------

"""
    changeSparseQFCASolver(name, printLevel)

Function used to change the solver and include the respective solver interfaces

# INPUT

- `name`:           Name of the solver

# OPTIONAL INPUT

- `printLevel`:     Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUT

- `model`:          A JuMP model for solving LPs.
- `solver`:         Solver object with a `handle` field

# EXAMPLES

```julia
julia> SolverName = :GLPK
julia> model, solver = changeCobraSolver(SolverName)
```

"""
function changeSparseQFCASolver(name, printLevel::Int=1)

    # Convert the input 'name' to a string if it's not already a string:
    if typeof(name) != :String
        name = string(name)
    end

    # Define an empty solver object with the provided name and initial handle set to 0:
    solver = SolverConfig(name, 0)

    ## Define the solver handle based on the solver name provided
    # If the solver is "HiGHS":
    if name == "HiGHS"
        try
            # Set the solver handle to HiGHS.Optimizer:
            solver.handle = HiGHS.Optimizer
            # Create a model using the HiGHS solver:
            model = Model(solver.handle)
            # Set specific attributes for the HiGHS solver:
            set_attribute(model, "presolve", "on")
            set_attribute(model, "time_limit", 60.0)
            # Set verbose attribute to false (disable verbose output):
            set_optimizer_attribute(model, "output_flag", false)
            # Set tolerances for increased accuracy
            set_optimizer_attribute(model, "primal_feasibility_tolerance", 1e-9)
            set_optimizer_attribute(model, "dual_feasibility_tolerance", 1e-9)
            # Return the created model and solver configuration:
            return model, solver
        catch
            # Handle the error if HiGHS cannot be set:
            error("The solver `HiGHS` cannot be set using `changeSparseQFCASolver()`.")
        end
#=
    # If the solver is "CPLEX":
    elseif name == "CPLEX"
        try
            # Set the solver handle to CPLEX.Optimizer:
            solver.handle = CPLEX.Optimizer
            # Create a model using the CPLEX solver:
            model = Model(solver.handle)
            # Set specific attributes for the CPLEX solver:
            set_attribute(model, "CPX_PARAM_EPINT", 1e-8)
            # Set verbose attribute to false (disable verbose output):
            set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
            # Return the created model and solver configuration:
            return model, solver
        catch
            # Handle the error if CPLEX cannot be set:
            error("The solver `CPLEX` cannot be set using `changeSparseQFCASolver()`.")
        end
=#
    # If the solver is "Clp":
    elseif name == "Clp"
        try
            # Set the solver handle to Clp.Optimizer:
            solver.handle = Clp.Optimizer
            # Create a model using the Clp solver:
            model = Model(solver.handle)
            # Set specific attributes for the Clp solver:
            set_attribute(model, "LogLevel", 1)
            set_attribute(model, "Algorithm", 4)
            # Return the created model and solver configuration:
            return model, solver
        catch
            # Handle the error if Clp cannot be set:
            error("The solver `Clp` cannot be set using `changeSparseQFCASolver()`.")
        end

    # If the solver is "Cbc":
    elseif name == "Cbc"
        try
            # Set the solver handle to Cbc.Optimizer:
            solver.handle = Cbc.Optimizer
            # Create a model using the Cbc solver:
            model = Model(solver.handle)
            # Set specific attributes for the Cbc solver:
            set_attribute(model, "logLevel", 1)
            # Return the created model and solver configuration:
            return model, solver
        catch
            # Handle the error if Cbc cannot be set:
            error("The solver `Cbc` cannot be set using `changeSparseQFCASolver()`.")
        end

    # If the solver is "GLPK":
    elseif name == "GLPK"
        try
            # Set the solver handle to GLPK.Optimizer:
            solver.handle = GLPK.Optimizer
            # Create a model using the GLPK solver:
            model = Model(solver.handle)
            # Set specific attributes for the GLPK solver:
            set_attribute(model, "tm_lim", 60 * 1_000)
            set_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
            # Return the created model and solver configuration:
            return model, solver

        catch
            # Handle the error if GLPK cannot be set:
            error("The solver `GLPK` cannot be set using `changeSparseQFCASolver()`.")
        end

    # If the solver is "ECOS":
    elseif name == "ECOS"
        try
            # Set the solver handle to ECOS.Optimizer:
            solver.handle = ECOS.Optimizer
            # Create a model using the ECOS solver:
            model = Model(solver.handle)
            # Set specific attributes for the ECOS solver:
            set_attribute(model, "maxit", 100)
            # Return the created model and solver configuration:
            return model, solver
        catch
            # Handle the error if ECOS cannot be set:
            error("The solver `ECOS` cannot be set using `changeSparseQFCASolver()`.")
        end

    # If the solver is "SCS":
    elseif name == "SCS"
        try
            # Set the solver handle to SCS.Optimizer:
            solver.handle = SCS.Optimizer
            # Create a model using the SCS solver:
            model = Model(solver.handle)
            # Return the created model and solver configuration:
            return model, solver
        catch
            # Handle the error if SCS cannot be set:
            error("The solver `SCS` cannot be set using `changeSparseQFCASolver()`.")
        end
#=
    # If the solver is "Mosek":
    elseif name == "Mosek"
        try
            # Set the solver handle to Mosek.Optimizer:
            solver.handle = Mosek.Optimizer
            # Create a model using the Mosek solver:
            model = Model(solver.handle)
            # Set specific attributes for the Mosek solver:
            set_attribute(model, "QUIET", true)
            set_attribute(model, "INTPNT_CO_TOL_DFEAS", 1e-7)
            # Return the created model and solver configuration:
            return model, solver
        catch
            # Handle the error if Mosek cannot be set:
            error("The solver `Mosek` cannot be set using `changeSparseQFCASolver()`.")
        end
=#
    # If the solver is not recognized or supported:
    else
        solver.handle = -1
        # Display an error message if the solver name is invalid:
        error("The solver is not supported. Please set the solver name to one the supported solvers.")
    end

end

end

#-----------------------------------------------------------------------------------------------------------------------------------------------------
