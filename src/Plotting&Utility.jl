function trapezoid_rule(f_sampled,dx)
    (dx/2)*sum(f_sampled[1:end-1] .+ f_sampled[2:end])
end

function construct_named_p(p_vector)
    p = (DN0 = p_vector[1],DL0 = p_vector[2],kN0 = p_vector[3],kL0 = p_vector[4],kE = p_vector[5],kNL = p_vector[6],σN0 = p_vector[7],σL0 = p_vector[8],Na = p_vector[9],NL = p_vector[10],NE = p_vector[11],mN = p_vector[12],mL = p_vector[13],mNL = p_vector[14],LN = p_vector[15],s0 = p_vector[16])
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


function plot_steady_state_profiles!(ax,ax_alpha,sol,linestyle)
    lines!(ax,sol.u[end][:,1],color = :blue,linestyle = linestyle)
    lines!(ax,sol.u[end][:,2],color = :orange,linestyle = linestyle)
    lines!(ax_alpha,sol.u[end][:,4],color = :red,linestyle = linestyle)
end

function save_dynamics_movie(sol)
    conc_N = Observable(sol.u[1][:,1])
    conc_L = Observable(sol.u[1][:,2])
    conc_α = Observable(sol.u[1][:,3])

    max_v = maximum(reduce(hcat,[sol for sol in sol.u]))

    fig = Figure(size = (800,800), fontsize = 20.)

    ax = Axis(fig[1,1], xlabel = L"\text{Tissue, } x", ylabel = L"\text{Signal (a.u.)}")
    ax1 = Axis(fig[1,1], yaxisposition = :right,ygridvisible = false,xgridvisible = false,yticklabelcolor = :red,ylabel = L"\alpha")

    hidespines!(ax1)
    hidexdecorations!(ax1)

    lines!(ax,tissue,conc_N)
    lines!(ax,tissue,conc_L)
    lines!(ax1,tissue,conc_α,color = :red)

    ylims!(ax,0.,max_v)

    record(fig, "MovieDyn.mp4", 1:length(sol.u)) do i
        conc_N[] = sol.u[i][:,1]
        conc_L[] = sol.u[i][:,2]
        conc_α[] = sol.u[i][:,3]
        sleep(0.05);
    end
end

function calculate_entropy(ptest,n_bin,lb,ub)

    # bin_boundaries = LinRange(lb[p_names_id[p_name]],ub[p_names_id[p_name]],n_bin+1)
    bin_boundaries = LinRange(lb,ub,n_bin+1)

    n_ptest = length(ptest)

    uni_prob = fit(Histogram, ptest, bin_boundaries; closed = :left) 

    ptest_bins = map(f->StatsBase.binindex(uni_prob, f),ptest);

    prob = [count(x->x==i,ptest_bins) / n_ptest for i in 1:n_bin]

    entropy(prob,n_bin)
end


metric_names = [:wt_t0,:cp_t0,:wt_xMax,:cp_xMax,:lm_xMax,:wt_d0,:cp_d0,:lm_d0,:xmax_peak_ratio,:xmax_mse,:alpha_mse,:cp_lprod_t0,:wt_lprod_t0,:retcodes]
metric_names_string = Dict(:wt_t0 => "wt_t0",:cp_t0=>"cp_t0",:wt_xMax=>"wt_xMax",:cp_xMax=>"cp_xMax",:lm_xMax=>"lm_xMax",:wt_d0=>"wt_d0",:cp_d0=>"cp_d0",:lm_d0=>"lm_d0",:xmax_peak_ratio=>"xmax_peak_ratio",:xmax_mse=>"xmax_mse",:alpha_mse=>"alpha_mse",:cp_lprod_t0=>"cp_lprod_t0",:wt_lprod_t0=>"wt_lprod_t0")

metric_names_latex = Dict(:wt_t0 => L"t_{\text{WT}}",:cp_t0 => L"t_{\text{wnt11}}",:wt_xMax=>L"X_{\text{max}}^{\text{WT}}",:cp_xMax=>L"X_{\text{max}}^{\text{wnt11}}",
:lm_xMax=>L"X_{\text{max}}^{\text{lefty}}",:wt_d0=>L"X_{\text{end}}^{\text{WT}} / X_{\text{max}}^{\text{WT}}",:cp_d0=>L"X_{\text{end}}^{\text{wnt11}} / X_{\text{max}}^{\text{wnt11}}",
:lm_d0=>L"X_{\text{end}}^{\text{lefty}} / X_{\text{max}}^{\text{lefty}}",:xmax_peak_ratio=>L"t_{\text{wnt11}} / t_{\text{WT}}",:xmax_mse => L"X_{\text{max}}^{MSE}",:alpha_mse=>L"\alpha_{\text{max}}^{MSE}",:cp_lprod_t0=>L"\text{Peak time L production - wnt11}",:wt_lprod_t0=>L"\text{Peak time L production - WT}")

transformations = Dict(:wt_t0 => t -> t / 60,:cp_t0 => t -> t / 60,:wt_xMax => t -> t ,:cp_xMax => t -> t ,:lm_xMax => t -> t ,:wt_d0 => t -> t ,:cp_d0 => t -> t ,:lm_d0 => t -> t ,:xmax_peak_ratio => t -> t ,:xmax_mse => t -> t,:alpha_mse=> t -> t,:wt_lprod_t0=>t->t / 60,:cp_lprod_t0=>t->t/60);

p_names = [:DN0,:DL0,:kN0,:kL0,:kE,:kNL,:σN0,:σL0,:Na,:NL,:NE,:LN,:s0]
p_names_id = Dict(:DN0=>1,:DL0=>2,:kN0=>3,:kL0=>4,:kE=>5,:kNL=>6,:σN0=>7,:σL0=>8,:Na=>9,:NL=>10,:NE=>11,:LN=>12,:s0=>13)

p_names_latex = Dict(:DN0 => L"D_N^0",:DL0 =>L"D_L^0",:kN0=>L"K_N^0",:kL0=>L"K_L^0",:kE => L"k_E",:kNL => L"k_{NL}",:σN0 => L"\sigma_{N}^{0}",:σL0 => L"\sigma_{L}^{0}",:Na => L"N_a",:NL => L"N_L",:NE => L"N_E",:LN => L"L_N",:s0 => L"s_0");
p_names_string = Dict(:DN0 => "D_N0",:DL0 => "D_L0",:kN0=>"K_N0",:kL0 =>"K_L0",:kE => "k_E",:kNL => "k_NL",:σN0 => "σN0",:σL0 => "σL0",:Na => "N_a",:NL => "N_L",:NE => "N_E",:LN => "L_N",:s0 => "s_0");