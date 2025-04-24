# This code is expected to be run from an sbatch script after a module load julia command has been run.
# It starts the remote processes with srun within an allocation specified in the sbatch script.

using Pkg

Pkg.activate("..")
Pkg.instantiate()
Pkg.precompile()

using DrWatson
using Distributed
using ClusterManagers

projectdir_static = dirname(Base.active_project())

cluster_calc = true

if cluster_calc
    n_tasks = parse(Int, ENV["SLURM_NTASKS"])
    addprocs(SlurmManager(n_tasks))
    @everywhere using Pkg
    @everywhere Pkg.activate("..")
end

@everywhere begin
    using DrWatson
    using JLD2
    using Printf
    using Base.Threads
    using Base.Threads: @spawn
    using ParallelDataTransfer

    using DifferentialEquations
    using StatsBase
    using Distributions
    using XLSX
    using DataFrames
    using Optimization,OptimizationNOMAD
end

@everywhere projectdirx(args...) = joinpath($projectdir_static, args...)

for dir_type ∈ ("data", "src", "plots", "scripts", "papers")
    function_name = Symbol(dir_type * "dirx")
    @everywhere @eval begin
        $function_name(args...) = projectdirx($dir_type, args...)
    end
end

@everywhere include(srcdirx("NodalLefty_E.jl"))
@everywhere include(srcdirx("FittingFunctions.jl"))

all_experiments = ["NodalLefty_DiffusionDominated"]

@everywhere include(scriptsdirx("LoadData.jl"))

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/" * $exp_name * ".jl"))

    p,p_cp,p_lm = get_params(pv_orig)

    tspan = (0,Inf)

    u0 = ones(Nc,4)
    
    u0[:,1] .= 1e-10
    u0[:,2] .= 1e-10
    u0[:,3] .= 0.
    u0[:,4] .= α0

    prob = ODEProblem(nodal_lefty_spatial_diff!,u0,tspan,p)

    lb = copy(pv_orig)
    ub = copy(pv_orig)

    lb[1:2] = 0.9 .* lb[1:2]
    ub[1:2] = 1.1 .* ub[1:2]

    lb[3:end] = (1 - γ) .* lb[3:end]
    ub[3:end] = (1 + γ) .* ub[3:end];

    lb_pow = copy(pv_orig)
    ub_pow = copy(pv_orig)

    lb_pow[1:2] = 0.9 .* lb_pow[1:2]
    ub_pow[1:2] = 1.1 .* ub_pow[1:2]

    lb_pow[3:end] = 10^(-j) .* lb_pow[3:end]
    ub_pow[3:end] = 10^(j) .* ub_pow[3:end];

    cp_list = [0.05,0.1,0.2,0.3]

    max_iter = 1000

    results = pmap(cp->optimize_params(prob,cp,pv_orig,lb,ub,max_iter),cp_list)
    results_pow = pmap(cp->optimize_params(prob,cp,pv_orig,lb_pow,ub_pow,max_iter),cp_list)

    summaryd = Dict{String, Any}()

    summaryd["cplist"] = cp_list 
    summaryd["OptimalParam"] = first.(results)
    summaryd["OptimalParam_obj"] = last.(results)
    summaryd["OptimalParam_pow"] = first.(results_pow)
    summaryd["OptimalParam_pow_obj"] = last.(results_pow)

    @tag!(summaryd)

    safesave(datadirx("exp_raw",exp_name * "_OptParmas.jld2"), summaryd)

    ########################################

end
