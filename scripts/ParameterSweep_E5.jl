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
    # using Optimization, OptimizationPolyalgorithms, SciMLSensitivity,OptimizationOptimJL,OptimizationBBO,OptimizationNOMAD
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

    # pv_orig = [DN0,DL0,kN0,kL0,kE,kNL,σN0,σL0,Na,NL,NE,LN,s0]

    lb = copy(pv_orig)
    ub = copy(pv_orig)

    lb = 0.1 .* lb
    ub = 10. .* ub

    var_id = [3,4,5,6,8,10,11,12]

    order_restr = [(11,10)]
    
    p_set =  generate_param_set(lb,ub,var_id,pv_orig,order_restr,N_sim)

    cp_set = [0.05,0.1,0.2]

    # sim = pmap(pv-> get_summary_metrics_safe(pv,prob,data,alpha_data,cp),p_set)

    sim = pmap(pv->get_summary_metrics_cpset(pv,prob,data,alpha_data,cp_set),p_set)

    summaryd = Dict{String, Any}()

    summaryd["Parameters"] = p_set
    summaryd["cp_set"] = cp_set
    summaryd["Results"] =  sim 

    @tag!(summaryd)

    safesave(datadirx("exp_raw",exp_name * "_CPSet_Sweep_RestrID_10fold.jld2"), summaryd)

    ########################################

end
