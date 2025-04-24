using Distributions

const N_samp = 5000
const t_grid_N = 100
const t_plot_N = 1000

function get_params(p_vector::Vector{Float64})

    p = (DN0 = p_vector[1],DL0 = p_vector[2],kN0 = p_vector[3],kL0 = p_vector[4],kE = p_vector[5],kNL = p_vector[6],σN0 = p_vector[7],σL0 = p_vector[8],Na = p_vector[9],NL = p_vector[10],NE = p_vector[11],mN = default_mN,mL = default_mL,mNL = default_mNL,LN = p_vector[12],s0 = p_vector[13])

    p_cp =  (DN0 = p_vector[1],DL0 = p_vector[2],kN0 = p_vector[3],kL0 = p_vector[4],kE = p_vector[5],kNL = p_vector[6],σN0 = p_vector[7],σL0 = p_vector[8],Na = p_vector[9],NL = p_vector[10],NE = 1e8,mN = default_mN,mL = default_mL,mNL = default_mNL,LN = p_vector[12],s0 = p_vector[13])

    p_lm =  (DN0 = p_vector[1],DL0 = p_vector[2],kN0 = p_vector[3],kL0 = p_vector[4],kE = p_vector[5],kNL = 0.,σN0 = p_vector[7],σL0 = p_vector[8],Na = p_vector[9],NL = p_vector[10],NE = p_vector[11],mN = default_mN,mL = default_mL,mNL = default_mNL,LN = p_vector[12],s0 = p_vector[13])

    p,p_cp,p_lm
end 

function get_params(p_orig::NamedTuple)

    p_cp =  (DN0 = p_orig[:DN0],DL0 = p[:DL0],kN0 = p[:kN0],kL0 = p[:kL0],kE = p[:kE],kNL = p[:kNL],σN0 = p[:σN0],σL0 = p[:σL0],Na = p[:Na],NL = p[:NL],NE = 1e8,mN = p[:mN],mL = p[:mL],mNL = p[:mNL],LN = p[:LN],s0 = p[:s0])

    p_lm =  (DN0 = p_orig[:DN0],DL0 = p[:DL0],kN0 = p[:kN0],kL0 = p[:kL0],kE = p[:kE],kNL = 0.,σN0 = p[:σN0],σL0 = p[:σL0],Na = p[:Na],NL = p[:NL],NE = p[:NE],mN = p[:mN],mL = p[:mL],mNL = p[:mNL],LN = p[:LN],s0 = p[:s0])

    p_orig,p_cp,p_lm
end 

function check_inbounds(pv,lb,ub)

    ib = all(pv.>=lb) & all(pv .<= ub)

    if ib
        return (true,[])
    else
        return(false,findall(pv.<lb),findall(pv .> ub))
    end
end

function get_lambda_half(sol,t_range)
    c0_t = [maximum(sol(t)[:,1]) for t in t_range] .- 1e-10

    λhalf_id = [findall(sol(t)[:,1] .- 1e-10 .> 0.5*c0) for (t,c0) in zip(t_range,c0_t)]
    λhalf_x = [length(id_list) != 0 ? tissue[maximum(id_list)] : 0. for id_list in λhalf_id];

    λhalf_x, t_range[argmax(λhalf_x)]
end

function mse_xmax_profiles(sol,sol_cp,wt_t0,c_level,xmax_profile_wt,xmax_profile_cp)

    level_x_wt_mse = get_level_x(sol,c_level,exp_times_times_norm .* wt_t0)  
    level_x_cp_mse = get_level_x(sol_cp,c_level,exp_times_times_norm .* wt_t0) 

    error = mean(((xmax_profile_wt  .-  level_x_wt_mse) ./ max_exp_wt).^2)
    error_cp = mean(((xmax_profile_cp  .-  level_x_cp_mse) ./ max_exp_wt).^2)

    return error,error_cp
end


function mse_alpha_profile(sol,t_grid_alpha,alpha_profiles)

    dyn_alpha = [sol(t)[alpha_x,4] for t in t_grid_alpha]

    error = sum([mean((dyn .- ap).^2) for (dyn,ap) in zip(dyn_alpha,alpha_profiles)])

    return error
end

function get_level_x(sol,c_level)
    c_level_id = [findall(sol[:,1] .>= c_level) for sol in sol.u]
    c_level_x = [length(id_list) != 0 ? tissue[maximum(id_list)] : 0. for id_list in c_level_id];
    return  c_level_x
end

function get_level_x(sol,c_level,t_grid)
    c_level_id = [findall(sol(t)[:,1] .>= c_level) for t in t_grid]
    c_level_x = [length(id_list) != 0 ? tissue[maximum(id_list)] : 0. for id_list in c_level_id];
    return  c_level_x
end


function get_summary_metrics(p_vector,prob,xmax_data,alpha_data,cp)

    p,p_cp,p_lm = get_params(p_vector)

    sol = solve(prob, p = p, FBDF(),abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_cp = solve(prob, p = p_cp, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_lm = solve(prob, p = p_lm, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));

    λ_trange = LinRange(0.,sol.t[end],N_samp)
    
    λhalf,λhalf_max_t = get_lambda_half(sol,λ_trange)
    
    c_max_wt = maximum(sol(λhalf_max_t)[:,1])

    c_level = cp*c_max_wt

    t_grid = LinRange(0,sol.t[end],t_grid_N)

    level_x_wt = get_level_x(sol,c_level,t_grid);
    level_x_cp = get_level_x(sol_cp,c_level,t_grid)
    level_x_lm = get_level_x(sol_lm,c_level,t_grid);

    wt_t0 = t_grid[argmax(level_x_wt)];
    cp_t0 = t_grid[argmax(level_x_cp)];

    wt_xMax = maximum(level_x_wt)
    cp_xMax = maximum(level_x_cp)
    lm_xMax = maximum(level_x_lm)

    wt_d0 = level_x_wt[end] ./ wt_xMax
    cp_d0 = level_x_cp[end] ./ cp_xMax
    lm_d0 = level_x_lm[end] ./ lm_xMax

    xmax_peak_ratio = cp_t0 / wt_t0 

    xmax_mse = mse_xmax_profiles(sol,sol_cp,wt_t0,c_level,xmax_data[:,"WT"],xmax_data[:,"SLB"])

    t_grid_alpha = alpha_data_times_norm .* wt_t0;

    alpha_mse = mse_alpha_profile(sol,t_grid_alpha,eachcol(alpha_data[:,2:end]))

    cp_lprod_t0,wt_lprod_t0 = get_integrated_lefty_prod(sol,sol_cp,t_grid)

    return (wt_t0 = wt_t0,cp_t0 = cp_t0,wt_xMax = wt_xMax,cp_xMax = cp_xMax,lm_xMax = lm_xMax,wt_d0 = wt_d0,cp_d0 = cp_d0,lm_d0 = lm_d0,xmax_peak_ratio = xmax_peak_ratio,xmax_mse = xmax_mse,alpha_mse = alpha_mse,cp_lprod_t0 = cp_lprod_t0,wt_lprod_t0 = wt_lprod_t0,retcodes = (sol.retcode,sol_cp.retcode,sol_lm.retcode))
end

function get_summary_metrics_cpset(p_vector,prob,xmax_data,alpha_data,cp_set)

    p,p_cp,p_lm = get_params(p_vector)

    sol = solve(prob, p = p, FBDF(),abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_cp = solve(prob, p = p_cp, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_lm = solve(prob, p = p_lm, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));

    λ_trange = LinRange(0.,sol.t[end],N_samp)
    
    λhalf,λhalf_max_t = get_lambda_half(sol,λ_trange)
    
    c_max_wt = maximum(sol(λhalf_max_t)[:,1])

    all_results = []

    for cp in cp_set

        c_level = cp*c_max_wt

        t_grid = LinRange(0,sol.t[end],t_grid_N)

        level_x_wt = get_level_x(sol,c_level,t_grid);
        level_x_cp = get_level_x(sol_cp,c_level,t_grid)
        level_x_lm = get_level_x(sol_lm,c_level,t_grid);

        wt_t0 = t_grid[argmax(level_x_wt)];
        cp_t0 = t_grid[argmax(level_x_cp)];

        wt_xMax = maximum(level_x_wt)
        cp_xMax = maximum(level_x_cp)
        lm_xMax = maximum(level_x_lm)

        wt_d0 = level_x_wt[end] ./ wt_xMax
        cp_d0 = level_x_cp[end] ./ cp_xMax
        lm_d0 = level_x_lm[end] ./ lm_xMax

        xmax_peak_ratio = cp_t0 / wt_t0 

        xmax_mse = mse_xmax_profiles(sol,sol_cp,wt_t0,c_level,xmax_data[:,"WT"],xmax_data[:,"SLB"])

        t_grid_alpha = alpha_data_times_norm .* wt_t0;

        alpha_mse = mse_alpha_profile(sol,t_grid_alpha,eachcol(alpha_data[:,2:end]))

        cp_lprod_t0,wt_lprod_t0 = get_integrated_lefty_prod(sol,sol_cp,t_grid)

        push!(all_results,(wt_t0 = wt_t0,cp_t0 = cp_t0,wt_xMax = wt_xMax,cp_xMax = cp_xMax,lm_xMax = lm_xMax,wt_d0 = wt_d0,cp_d0 = cp_d0,lm_d0 = lm_d0,xmax_peak_ratio = xmax_peak_ratio,xmax_mse = xmax_mse,alpha_mse = alpha_mse,cp_lprod_t0 = cp_lprod_t0,wt_lprod_t0 = wt_lprod_t0,retcodes = (sol.retcode,sol_cp.retcode,sol_lm.retcode)))
    end

    return all_results
end

function get_integrated_lefty_prod(sol,sol_cp,t_grid)
    cNt_cp = [sol_cp(t)[:,1] for t in t_grid]
    αt_cp = [sol_cp(t)[:,4] for t in t_grid]

    cNt = [sol(t)[:,1] for t in t_grid]
    αt = [sol(t)[:,4] for t in t_grid]

    νN_int_cp = [trapezoid_rule(ν.(cN,σ.(σL0,ϕ0,ϕ.(α)),NL,mL),1) for (cN,α) in zip(cNt_cp,αt_cp)];
    νN_int = [trapezoid_rule(ν.(cN,σ.(σL0,ϕ0,ϕ.(α)),NL,mL),1) for (cN,α) in zip(cNt,αt)];

    return (t_grid[argmax(νN_int_cp)] ,t_grid[argmax(νN_int)])
end

function get_summary_metrics_safe(p_vector,prob,xmax_data,alpha_data,cp)

    p,p_cp,p_lm = get_params(p_vector)

    try
        get_summary_metrics(p_vector,prob,xmax_data,alpha_data,cp)
    catch 
        return p
    end
end

function get_summary_metrics_cpset_safe(p_vector,prob,xmax_data,alpha_data,cp_set)

    p,p_cp,p_lm = get_params(p_vector)

    try
        get_summary_metrics_cpset(p_vector,prob,xmax_data,alpha_data,cp_set)
    catch 
        return p
    end
end


function get_alpha_xmax_lambda(p_vector,prob,cp)

    p,p_cp,p_lm = get_params(p_vector)

    sol = solve(prob, p = p, FBDF(),abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_cp = solve(prob, p = p_cp, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_lm = solve(prob, p = p_lm, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));

    λ_trange = LinRange(0.,sol.t[end],N_samp)
    
    λhalf,λhalf_max_t = get_lambda_half(sol,λ_trange)
    
    c_max_wt = maximum(sol(λhalf_max_t)[:,1])

    c_level = cp*c_max_wt

    level_x_wt = get_level_x(sol,c_level,λ_trange);

    wt_t0 = λ_trange[argmax(level_x_wt)];

    t_plot = LinRange(0,exp_times_times_norm[end],t_plot_N)

    level_x_wt_rescaled = get_level_x(sol,c_level,t_plot .* wt_t0)  
    level_x_cp_rescaled = get_level_x(sol_cp,c_level,t_plot .* wt_t0)
    level_x_lm_rescaled  = get_level_x(sol_lm,c_level,t_plot .* wt_t0)

    t_grid_alpha = alpha_data_times_norm .* wt_t0;

    dyn_alpha = [sol(t)[alpha_x,4] for t in t_grid_alpha]

    porosity_dyn = [mean(ϕ.(sol(t)[1:50,4])) for t in t_plot .* wt_t0];

    porosity_dyn_cp =   [mean(ϕ.(sol_cp(t)[1:50,4])) for t in t_plot .* wt_t0];

    return (t_grid_alpha,dyn_alpha),(t_plot,(level_x_wt_rescaled,level_x_cp_rescaled,level_x_lm_rescaled )),(porosity_dyn,porosity_dyn_cp),c_level,(sol,sol_cp,sol_lm)

end

function loss(p_vector,prob,xmax_data,alpha_data,cp)

    p,p_cp,p_lm = get_params(p_vector)

    sol = solve(prob, p = p, FBDF(),abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_cp = solve(prob, p = p_cp, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));

    λ_trange = LinRange(0.,sol.t[end],N_samp)
    
    λhalf,λhalf_max_t = get_lambda_half(sol,λ_trange)
    
    c_max_wt = maximum(sol(λhalf_max_t)[:,1])

    c_level = cp*c_max_wt

    t_grid = LinRange(0,sol.t[end],t_grid_N)

    level_x_wt = get_level_x(sol,c_level,t_grid);

    wt_t0 = t_grid[argmax(level_x_wt)];

    xmax_mse = mse_xmax_profiles(sol,sol_cp,wt_t0,c_level,xmax_data[:,"WT"],xmax_data[:,"SLB"])

    t_grid_alpha = alpha_data_times_norm .* wt_t0;

    alpha_mse = mse_alpha_profile(sol,t_grid_alpha,[c for c in eachcol(alpha_data[:,2:end])])

    return sum(xmax_mse) + alpha_mse
end

function metric_loss(p_vector,prob,data_metrics,alpha_data,cp)

    p,p_cp,p_lm = get_params(p_vector)

    sol = solve(prob, p = p, FBDF(),abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_cp = solve(prob, p = p_cp, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));
    sol_lm = solve(prob, p = p_lm, FBDF(),abstol = de_abstol,reltol = de_reltol,maxiters = 1e6,callback = TerminateSteadyState(1e-6,1e-4),isoutofdomain = (u,p,t) -> any(x->x<0, u));

    λ_trange = LinRange(0.,sol.t[end],N_samp)
    
    λhalf,λhalf_max_t = get_lambda_half(sol,λ_trange)
    
    c_max_wt = maximum(sol(λhalf_max_t)[:,1])

    c_level = cp*c_max_wt

    t_grid = LinRange(0,sol.t[end],t_grid_N)

    level_x_wt = get_level_x(sol,c_level,t_grid);
    level_x_cp = get_level_x(sol_cp,c_level,t_grid)
    level_x_lm = get_level_x(sol_lm,c_level,t_grid);

    wt_t0 = t_grid[argmax(level_x_wt)];
    cp_t0 = t_grid[argmax(level_x_cp)];

    wt_xMax = maximum(level_x_wt)
    cp_xMax = maximum(level_x_cp)
    lm_xMax = maximum(level_x_lm)

    wt_d0 = level_x_wt[end] ./ wt_xMax
    cp_d0 = level_x_cp[end] ./ cp_xMax
    lm_d0 = level_x_lm[end] ./ lm_xMax

    xmax_peak_ratio = cp_t0 / wt_t0 

    t_grid_alpha = alpha_data_times_norm .* wt_t0;

    alpha_mse = mse_alpha_profile(sol,t_grid_alpha,eachcol(alpha_data[:,2:end]))

    # cp_lprod_t0,wt_lprod_t0 = get_integrated_lefty_prod(sol,sol_cp,t_grid)

    return (wt_d0 - data_metrics[:wt_d0])^2 + (cp_d0 - data_metrics[:cp_d0])^2 + (xmax_peak_ratio - data_metrics[:xmax_peak_ratio])^2 + alpha_mse  + ((wt_xMax/data_metrics[:wt_xMax]) - 1)^2 + ((cp_xMax/data_metrics[:cp_xMax]) - 1)^2 
end

function metric_loss_safe(p_vector,prob,data_metrics,alpha_data,cp)
    try 
        metric_loss(p_vector,prob,data_metrics,alpha_data,cp)
    catch
        1e8
    end
end

function loss_safe(p_vector,prob,xmax_data,alpha_data,cp)
    try 
        loss(p_vector,prob,xmax_data,alpha_data,cp)
    catch
        1e8
    end
end

function generate_param_set(lb,ub,N)
    all_p = zeros(N,length(lb))
    for (n,(lbi,ubi)) in enumerate(zip(lb,ub))
        all_p[:,n] .= rand(Uniform(lbi,ubi),N)
    end

    return [collect(p) for p in eachrow(all_p)]
end

function optimize_params(prob,cp,pv_orig,lb,ub,max_iter)

    optf = Optimization.OptimizationFunction((x, p) -> loss_safe(x,prob,data,alpha_data,cp))
        
    optprob = Optimization.OptimizationProblem(optf,pv_orig,lb = lb, ub = ub);

    result_opt = Optimization.solve(optprob,NOMADOpt(),maxiters = max_iter);

    return result_opt.u,result_opt.objective
end
