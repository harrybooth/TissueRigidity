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

# plotting function to compare a single parameter set 

function plot_summary_newtimes!(fig,pv,prob)

    # try

        p_tuple,p_cp_tuple,p_lm_tuple,p_ro_tuple = get_params_ro(pv)

        p_orig,_,_ = get_params(pv_orig)

        ax = Axis(fig[1,1], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Position of 20% max } c_N \text{ } (μm) ",ygridvisible = false,xgridvisible = false,title = "w/ cmax 0.2")
        ax_por = Axis(fig[1,1], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Average porosity, } \phi", yaxisposition = :right,ylabelcolor = :red,yticklabelcolor = :red,ygridvisible = false,xgridvisible = false,xticksvisible = false)

        (t_grid_alpha,dyn_alpha),(t_plot,(level_x_wt_rescaled,level_x_cp_rescaled,level_x_lm_rescaled )),(porosity_dyn,porosity_dyn_cp),c_level,(sol,sol_cp,sol_lm),(λhalf_t,λhalf_ro_t),(λhalf,λhalf_ro),λ_trange = get_alpha_xmax_lambda(pv,prob,0.2);

        hidexdecorations!(ax_por)

        lines!(ax,t_plot,level_x_wt_rescaled,color = :black,label = L"\text{WT}")
        lines!(ax,t_plot,level_x_cp_rescaled,color = :orange,label = L"\text{Wnt11}")
        lines!(ax,t_plot,level_x_lm_rescaled,color = :pink,label = L"\text{Lefty mutant}")

        lines!(ax_por,t_plot,porosity_dyn ,linestyle = :dash,color = :black,label = L"\text{ϕ}")
        lines!(ax_por,t_plot,porosity_dyn_cp ,linestyle = :dash,color = :orange,label = L"\text{ϕ}")

        axislegend(ax,position = :lt)

        ylims!(ax_por,0.,0.14)

        ax.xticks = (0:0.5:3.5,string.(0:0.5:3.5))

        orig_metrics = get_summary_metrics(pv,prob,data,alpha_data,0.2)

        t_plot_int = LinRange(0,3*orig_metrics[:wt_t0],t_plot_N)

        νN_int_cp,νN_int = get_integrated_lefty_prod_values(sol,sol_cp,t_plot_int,p_tuple,p_cp_tuple)

        # nodal_prod_cp,nodal_prod = get_integrated_nodal_prod_values(sol,sol_cp,t_plot_int)

        ax = Axis(fig[1,2], xlabel = L"\text{Rescaled time} t^{'}", ylabel= L"\int_0^L ν_L(N(t,x)) dx" )

        lines!(ax,LinRange(0,3,1000),νN_int_cp,linestyle = :dash,color = :grey, label = L"\text{Wnt11}")
        lines!(ax,LinRange(0,3,1000),νN_int,color = :grey, label = L"\text{WT}")
        
        axislegend(ax,position = :rt)

        ax = Axis(fig[1,3], xlabel = L"\text{Position}", ylabel= L"α(t)", title = "Beta = " * string(β) )

        for (n,d) in enumerate(dyn_alpha[1:end-1])
            lines!(ax,alpha_x,d, label = string(alpha_data_times_norm[n])* "* t_wt")
            scatter!(ax,alpha_x,alpha_data[:,n+1])
        end

        axislegend(ax,position = :rb)

        axr = Axis(fig[6,3], xlabel = L"\text{Position}", ylabel= L"ϕ(α(t,x))", title = "Beta = " * string(β) )

        for (n,d) in enumerate(dyn_alpha[1:end-1])
            lines!(axr,alpha_x,ϕ.(d), label = string(alpha_data_times_norm[n])* "* t_wt")
        end

        axislegend(ax,position = :rb)

        alpha_data_times_norm_v = [0.,0.3,0.6,0.9,1.2,1.5,1.8]

        t_check =  alpha_data_times_norm_v[2:end] .* orig_metrics[:wt_t0]

        prob_finite = remake(prob,tspan = (0, alpha_data_times_norm_v[end] * orig_metrics[:wt_t0]))

        sol_profiles = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);
        sol_profiles_cp = solve(prob_finite, p = p_cp_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);
        sol_profiles_lm = solve(prob_finite, p = p_lm_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);

        dyn_N = [sol[:,1] for sol in sol_profiles.u[1:length(t_check)]]
        dyn_L = [sol[:,2] for sol in sol_profiles.u[1:length(t_check)]];
        dyn_α = [sol[:,4] for sol in sol_profiles.u[1:length(t_check)]];

        dyn_N_cp = [sol[:,1] for sol in sol_profiles_cp.u[1:length(t_check)]];
        dyn_L_cp = [sol[:,2] for sol in sol_profiles_cp.u[1:length(t_check)]];
        dyn_α_cp = [sol[:,4] for sol in sol_profiles_cp.u[1:length(t_check)]];

        dyn_N_lm = [sol[:,1] for sol in sol_profiles_lm.u[1:length(t_check)]];

        lefty_prod_profiles_cp = [ν.(cN,σ.(p_cp_tuple[:σL0],ϕ0,ϕ.(α)),p_cp_tuple[:NL],p_cp_tuple[:mL]) for (cN,α) in zip(dyn_N_cp,dyn_α_cp)];
        lefty_prod_profiles = [ν.(cN,σ.(p_tuple[:σL0],ϕ0,ϕ.(α)),p_tuple[:NL],p_tuple[:mL]) for (cN,α) in zip(dyn_N,dyn_α)];

        nodal_prod_profiles_cp = [ν.(cN,σ.(p_cp_tuple[:σN0],ϕ0,ϕ.(α)),p_cp_tuple[:Na],p_cp_tuple[:mN]) for (cN,α) in zip(dyn_N_cp,dyn_α_cp)];
        nodal_prod_profiles = [ν.(cN,σ.(p_tuple[:σN0],ϕ0,ϕ.(α)),p_tuple[:Na],p_tuple[:mN]) for (cN,α) in zip(dyn_N,dyn_α)];

        # dyn_N = [sol[:,1] for sol in sol_profiles.u]
        # dyn_L = [sol[:,2] for sol in sol_profiles.u];

        # dyn_N_cp = [sol[:,1] for sol in sol_profiles_cp.u];
        # dyn_L_cp = [sol[:,2] for sol in sol_profiles_cp.u];

        ax1 = Axis(fig[2,1], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = "WT")

        for (i,N) in enumerate(dyn_N)
            lines!(ax1,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        ax2 = Axis(fig[2,2], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = "Wnt11")

        for (i,N) in enumerate(dyn_N_cp)
            lines!(ax2,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax2,position = :rt)

        ax3 = Axis(fig[2,3], xlabel = L"\text{Position}", ylabel= L"\text{Lefty Concentration}", title = "WT")

        for (i,N) in enumerate(dyn_L)
            lines!(ax3,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax3,position = :rt)

        ax4 = Axis(fig[3,1], xlabel = L"\text{Position}", ylabel= L"\text{Lefty Concentration}", title = "Wnt11")

        for (i,N) in enumerate(dyn_L_cp)
            lines!(ax4,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax4,position = :rt)

        ax5 = Axis(fig[6,2], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = "Lefty mutant")

        for (i,N) in enumerate(dyn_N_lm)
            lines!(ax5,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax5,position = :rt)

        axt1 = Axis(fig[3,2])
        axt2 = Axis(fig[3,3])

        hidedecorations!(axt1)
        hidedecorations!(axt2)

        text_pos = [Point(-0.5,i) for i in LinRange(-1,1,6)]

        for (n,p) in enumerate(p_names[1:6])

            if p_tuple[p] != p_orig[p]
                col = :red
            else
                col = :black
            end

            if p == :s0
                text!(axt1,text_pos[n], text = p_names_string[p] * " = " * string(round(2*p_tuple[p],digits = 15)),color = col)
            else
                text!(axt1,text_pos[n], text = p_names_string[p] * " = " * string(round(p_tuple[p],digits = 15)),color = col)
            end
        end

        ylims!(axt1,-1.8,1.8)
        ylims!(axt2,-1.8,1.8)
        xlims!(axt1,-1.5,1.5)
        xlims!(axt2,-1.5,1.5)

        text_pos = [Point(-0.5,i) for i in LinRange(-1,1,7)]

        for (n,p) in enumerate(p_names[7:end])

            if p_tuple[p] != p_orig[p]
                col = :red
            else
                col = :black
            end

            if p == :s0
                text!(axt2,text_pos[n], text = p_names_string[p] * " = " * string(round(2*p_tuple[p],digits = 15)),color = col)
            else
                text!(axt2,text_pos[n], text = p_names_string[p] * " = " * string(round(p_tuple[p],digits = 15)),color = col)
            end
        end

        #######

        ax1 = Axis(fig[4,1], xlabel = L"\text{Position}", ylabel= L"\text{Nodal production}", title = "WT")

        for (i,N) in enumerate(nodal_prod_profiles)
            lines!(ax1,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        ax2 = Axis(fig[4,2], xlabel = L"\text{Position}", ylabel= L"\text{Nodal production}", title = "Wnt11")

        for (i,N) in enumerate(nodal_prod_profiles_cp)
            lines!(ax2,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax2,position = :rt)

        ax3 = Axis(fig[4,3], xlabel = L"\text{Position}", ylabel= L"\text{Lefty production}", title = "WT")

        for (i,N) in enumerate(lefty_prod_profiles)
            lines!(ax3,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax3,position = :rt)

        ax4 = Axis(fig[5,1], xlabel = L"\text{Position}", ylabel= L"\text{Lefty production}", title = "Wnt11")

        for (i,N) in enumerate(lefty_prod_profiles_cp)
            lines!(ax4,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        nodal_degr = [p_tuple[:kNL] .*(1 ./ (1 .+ (p_tuple[:LN] ./ L) .^ p_tuple[:mNL])) for L in dyn_L]
        nodal_degr_cp = [p_cp_tuple[:kNL].*(1 ./ (1 .+ (p_cp_tuple[:LN] ./ L) .^ p_cp_tuple[:mNL])) for L in dyn_L_cp]

        ax1 = Axis(fig[5,2], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = "WT")

        for (i,N) in enumerate(nodal_degr)
            lines!(ax1,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        ax2 = Axis(fig[5,3], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = "Wnt11")

        for (i,N) in enumerate(nodal_degr_cp)
            lines!(ax2,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        ts = max(λhalf_t,λhalf_ro_t)
    
        prob_finite = remake(prob,tspan = (0, 3*ts))

        λ_trange = LinRange(0.,3*ts,10000)

        sol1  = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);
        sol_ro1  = solve(prob_finite, p = p_ro_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);

        N0t = [sol(t)[1,1] for t in λ_trange]     

        λhalf,λhalf_max_t = get_lambda_half(sol1,λ_trange)
        λhalf_ro,λhalf_max_t_ro = get_lambda_half(sol_ro1,λ_trange)

        text!(axt1,-0.5,-1.5,text = "normal : " * string(round(λhalf_max_t,digits = 2)),color = :blue)
        text!(axt2,-0.5,-1.5,text = "ro : " * string(round(λhalf_max_t_ro,digits = 2)),color = :green)

        text!(axt2,-0.5,-1.7,text = "ratio ro/normal : " * string(round(λhalf_max_t_ro / λhalf_max_t,digits = 2)),color = :orange)

        ax = Axis(fig[6,1], xlabel = L"\text{time (mins)}", ylabel= L"\lambda_{1/2} \quad (\mu m)",ygridvisible = false,xgridvisible = false)

        ax_N = Axis(fig[6,1], xlabel = L"\text{Rescaled time, t}", ylabel= L"N(x=0,t)", yaxisposition = :right,ylabelcolor = :orange,yticklabelcolor = :orange,ygridvisible = false,xgridvisible = false,xticksvisible = false)

        hidexdecorations!(ax_N)

        vlines!(ax,λhalf_max_t,color = :blue,linestyle = :dash)
        vlines!(ax,λhalf_max_t_ro,color = :green,linestyle = :dash)

        vlines!(ax_N,λhalf_t,color = :orange,linestyle = :dash)

        lines!(ax,λ_trange,λhalf,color = :blue, label = "WT")

        lines!(ax_N,λ_trange,N0t,color = :orange, label = "WT")

        lines!(ax,λ_trange,λhalf_ro,color = :green, label = "relay off (WT)")

        return fig

    # catch DomainError
    #     return "Numerical instability. Increase abstol/reltol"
    # end

end

# plotting functions to compare a list of parameters in pv_list (e.g. for relay to diffusion SM analysis)

function plot_summary_newtimes_compare!(fig,pv_list,prob)

    legend_labels = [L"0.3t_{\text{WT}}",L"0.6t_{\text{WT}}",L"0.9t_{\text{WT}}",L"1.2t_{\text{WT}}",L"1.5t_{\text{WT}}",L"1.8t_{\text{WT}}"]

    ax_smad = []
    ax_smad_por = []
    ax_int_lefty = []
    ax_np_wt = []
    ax_np_wnt11 = []
    ax_n_prod_wt = []
    ax_n_prod_wnt11 = []
    ax_lfty_deg_wt = []
    ax_lfty_deg_wnt11 = []
    ax_lamba_t = []

    for (n,pv) in enumerate(pv_list)

        p_tuple,p_cp_tuple,p_lm_tuple,p_ro_tuple = get_params_ro(pv)

        p_orig,_,_ = get_params(pv_orig)

        ax = Axis(fig[1,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Position of 20% max } c_N \text{ } (μm) ",ygridvisible = false,xgridvisible = false,title = "w/ cmax 0.2")
        ax_por = Axis(fig[1,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Average porosity, } \phi", yaxisposition = :right,ylabelcolor = :red,yticklabelcolor = :red,ygridvisible = false,xgridvisible = false,xticksvisible = false)

        (t_grid_alpha,dyn_alpha),(t_plot,(level_x_wt_rescaled,level_x_cp_rescaled,level_x_lm_rescaled )),(porosity_dyn,porosity_dyn_cp),c_level,(sol,sol_cp,sol_lm),(λhalf_t,λhalf_ro_t),(λhalf,λhalf_ro),λ_trange = get_alpha_xmax_lambda(pv,prob,0.2);

        hidexdecorations!(ax_por)

        lines!(ax,t_plot,level_x_wt_rescaled,color = :black,label = L"\text{WT}")
        lines!(ax,t_plot,level_x_cp_rescaled,color = :orange,label = L"\text{Wnt11}")
        lines!(ax,t_plot,level_x_lm_rescaled,color = :pink,label = L"\text{Lefty mutant}")

        lines!(ax_por,t_plot,porosity_dyn ,linestyle = :dash,color = :black,label = L"\text{ϕ}")
        lines!(ax_por,t_plot,porosity_dyn_cp ,linestyle = :dash,color = :orange,label = L"\text{ϕ}")

        axislegend(ax,position = :lt)

        ylims!(ax_por,0.,0.14)

        ax.xticks = (0:0.5:3.5,string.(0:0.5:3.5))

        push!(ax_smad,ax)
        push!(ax_smad_por,ax_por)

        orig_metrics = get_summary_metrics(pv,prob,data,alpha_data,0.2)

        t_plot_int = LinRange(0,3*orig_metrics[:wt_t0],t_plot_N)

        νN_int_cp,νN_int = get_integrated_lefty_prod_values(sol,sol_cp,t_plot_int,p_tuple,p_cp_tuple)

        # nodal_prod_cp,nodal_prod = get_integrated_nodal_prod_values(sol,sol_cp,t_plot_int)

        ax = Axis(fig[2,n], xlabel = L"\text{Rescaled time} t^{'}", ylabel= L"\int_0^L ν_L(N(t,x)) dx" )

        lines!(ax,LinRange(0,3,1000),νN_int_cp,linestyle = :dash,color = :grey, label = L"\text{Wnt11}")
        lines!(ax,LinRange(0,3,1000),νN_int,color = :grey, label = L"\text{WT}")
        
        axislegend(ax,position = :rt)

        ax = Axis(fig[3,n], xlabel = L"\text{Position}", ylabel= L"α(t)", title = "Beta = " * string(β) )

        for (n,d) in enumerate(dyn_alpha[1:end-1])
            lines!(ax,alpha_x,d, label = string(alpha_data_times_norm[n])* "* t_wt")
            scatter!(ax,alpha_x,alpha_data[:,n+1])
        end

        axislegend(ax,position = :rb)

        push!(ax_int_lefty,ax)

        alpha_data_times_norm_v = [0.,0.3,0.6,0.9,1.2,1.5,1.8]

        t_check =  alpha_data_times_norm_v[2:end] .* orig_metrics[:wt_t0]

        prob_finite = remake(prob,tspan = (0, alpha_data_times_norm_v[end] * orig_metrics[:wt_t0]))

        sol_profiles = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);
        sol_profiles_cp = solve(prob_finite, p = p_cp_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);

        dyn_N = [sol[:,1] for sol in sol_profiles.u[1:length(t_check)]]
        dyn_L = [sol[:,2] for sol in sol_profiles.u[1:length(t_check)]];
        dyn_α = [sol[:,4] for sol in sol_profiles.u[1:length(t_check)]];

        dyn_N_cp = [sol[:,1] for sol in sol_profiles_cp.u[1:length(t_check)]];
        dyn_L_cp = [sol[:,2] for sol in sol_profiles_cp.u[1:length(t_check)]];
        dyn_α_cp = [sol[:,4] for sol in sol_profiles_cp.u[1:length(t_check)]];

        lefty_prod_profiles_cp = [ν.(cN,σ.(p_cp_tuple[:σL0],ϕ0,ϕ.(α)),p_cp_tuple[:NL],p_cp_tuple[:mL]) for (cN,α) in zip(dyn_N_cp,dyn_α_cp)];
        lefty_prod_profiles = [ν.(cN,σ.(p_tuple[:σL0],ϕ0,ϕ.(α)),p_tuple[:NL],p_tuple[:mL]) for (cN,α) in zip(dyn_N,dyn_α)];

        nodal_prod_profiles_cp = [ν.(cN,σ.(p_cp_tuple[:σN0],ϕ0,ϕ.(α)),p_cp_tuple[:Na],p_cp_tuple[:mN]) for (cN,α) in zip(dyn_N_cp,dyn_α_cp)];
        nodal_prod_profiles = [ν.(cN,σ.(p_tuple[:σN0],ϕ0,ϕ.(α)),p_tuple[:Na],p_tuple[:mN]) for (cN,α) in zip(dyn_N,dyn_α)];

        # dyn_N = [sol[:,1] for sol in sol_profiles.u]
        # dyn_L = [sol[:,2] for sol in sol_profiles.u];

        # dyn_N_cp = [sol[:,1] for sol in sol_profiles_cp.u];
        # dyn_L_cp = [sol[:,2] for sol in sol_profiles_cp.u];

        ax1 = Axis(fig[4,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = "WT")

        for (i,N) in enumerate(dyn_N)
            lines!(ax1,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        push!(ax_np_wt,ax1)

        ax2 = Axis(fig[5,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = "Wnt11")

        for (i,N) in enumerate(dyn_N_cp)
            lines!(ax2,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax2,position = :rt)

        push!(ax_np_wnt11,ax2)

        # ax3 = Axis(fig[2,3], xlabel = L"\text{Position}", ylabel= L"\text{Lefty Concentration}", title = "WT")

        # for (i,N) in enumerate(dyn_L)
        #     lines!(ax3,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        # end

        # axislegend(ax3,position = :rt)

        # ax4 = Axis(fig[3,1], xlabel = L"\text{Position}", ylabel= L"\text{Lefty Concentration}", title = "Wnt11")

        # for (i,N) in enumerate(dyn_L_cp)
        #     lines!(ax4,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        # end

        # axislegend(ax4,position = :rt)

        axt1 = Axis(fig[11,n])
        axt2 = Axis(fig[12,n])

        hidedecorations!(axt1)
        hidedecorations!(axt2)

        text_pos = [Point(-0.5,i) for i in LinRange(-1,1,6)]

        for (n,p) in enumerate(p_names[1:6])

            if p_tuple[p] != p_orig[p]
                col = :red
            else
                col = :black
            end

            if p == :s0
                text!(axt1,text_pos[n], text = p_names_string[p] * " = " * string(round(2*p_tuple[p],digits = 8)),color = col)
            else
                text!(axt1,text_pos[n], text = p_names_string[p] * " = " * string(round(p_tuple[p],digits = 8)),color = col)
            end
        end

        # text!(axt1,-0.5,-1.5,text = "normal : " * string(round(λhalf_t,digits = 2)),color = :blue)

        ylims!(axt1,-1.8,1.8)
        ylims!(axt2,-1.8,1.8)
        xlims!(axt1,-1.5,1.5)
        xlims!(axt2,-1.5,1.5)

        text_pos = [Point(-0.5,i) for i in LinRange(-1,1,7)]

        for (n,p) in enumerate(p_names[7:end])

            if p_tuple[p] != p_orig[p]
                col = :red
            else
                col = :black
            end

            if p == :s0
                text!(axt2,text_pos[n], text = p_names_string[p] * " = " * string(round(2*p_tuple[p],digits = 8)),color = col)
            else
                text!(axt2,text_pos[n], text = p_names_string[p] * " = " * string(round(p_tuple[p],digits = 8)),color = col)
            end
        end

        # text!(axt2,-0.5,-1.5,text = "ro : " * string(round(λhalf_ro_t,digits = 2)),color = :green)

        # text!(axt2,-0.5,-1.7,text = "ratio ro/normal : " * string(round(λhalf_ro_t / λhalf_t,digits = 2)),color = :orange)

        #######

        ax1 = Axis(fig[6,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal production}", title = "WT")

        for (i,N) in enumerate(nodal_prod_profiles)
            lines!(ax1,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        push!(ax_n_prod_wt,ax1)

        ax2 = Axis(fig[7,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal production}", title = "Wnt11")

        for (i,N) in enumerate(nodal_prod_profiles_cp)
            lines!(ax2,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax2,position = :rt)

        push!(ax_n_prod_wnt11,ax2)

        # ax3 = Axis(fig[8,n], xlabel = L"\text{Position}", ylabel= L"\text{Lefty production}", title = "WT")

        # for (i,N) in enumerate(lefty_prod_profiles)
        #     lines!(ax3,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        # end

        # axislegend(ax3,position = :rt)

        # ax4 = Axis(fig[5,1], xlabel = L"\text{Position}", ylabel= L"\text{Lefty production}", title = "Wnt11")

        # for (i,N) in enumerate(lefty_prod_profiles_cp)
        #     lines!(ax4,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        # end

        nodal_degr = [p_tuple[:kNL] .*(1 ./ (1 .+ (p_tuple[:LN] ./ L) .^ p_tuple[:mNL])) for L in dyn_L]
        nodal_degr_cp = [p_cp_tuple[:kNL].*(1 ./ (1 .+ (p_cp_tuple[:LN] ./ L) .^ p_cp_tuple[:mNL])) for L in dyn_L_cp]

        ax1 = Axis(fig[8,n], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = "WT")

        for (i,N) in enumerate(nodal_degr)
            lines!(ax1,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        push!(ax_lfty_deg_wt,ax1)

        ax2 = Axis(fig[9,n], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = "Wnt11")

        for (i,N) in enumerate(nodal_degr_cp)
            lines!(ax2,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        push!(ax_lfty_deg_wnt11,ax2)

        ts = max(λhalf_t,λhalf_ro_t)
    
        prob_finite = remake(prob,tspan = (0, 3*ts))

        λ_trange = LinRange(0.,3*ts,10000)

        sol1  = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);
        sol_ro1  = solve(prob_finite, p = p_ro_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);

        N0t = [sol1(t)[1,1] for t in λ_trange]    
        
        λhalf,λhalf_max_t = get_lambda_half(sol1,λ_trange)
        λhalf_ro,λhalf_max_t_ro = get_lambda_half(sol_ro1,λ_trange)

        ax = Axis(fig[10,n], xlabel = L"\text{time (mins)}", ylabel= L"\lambda_{1/2} \quad (\mu m)",ygridvisible = false,xgridvisible = false)

        ax_N = Axis(fig[10,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"N(x=0,t)", yaxisposition = :right,ylabelcolor = :orange,yticklabelcolor = :orange,ygridvisible = false,xgridvisible = false,xticksvisible = false)

        lines!(ax_N,λ_trange,N0t,color = :orange, label = "WT")

        text!(axt1,-0.5,-1.5,text = "normal : " * string(round(λhalf_max_t,digits = 2)),color = :blue)
        text!(axt2,-0.5,-1.5,text = "ro : " * string(round(λhalf_max_t_ro,digits = 2)),color = :green)

        text!(axt2,-0.5,-1.7,text = "ratio ro/normal : " * string(round(λhalf_max_t_ro / λhalf_max_t,digits = 2)),color = :orange)

        vlines!(ax,λhalf_max_t,color = :blue,linestyle = :dash)
        vlines!(ax,λhalf_max_t_ro,color = :green,linestyle = :dash)

        lines!(ax,λ_trange,λhalf,color = :blue, label = "WT")

        lines!(ax,λ_trange,λhalf_ro,color = :green, label = "relay off (WT)")

        push!(ax_lamba_t,ax)

    end

    for ax_set in  [ax_smad,ax_smad_por,ax_int_lefty,ax_np_wt,ax_np_wnt11,ax_n_prod_wt,ax_n_prod_wnt11,ax_lfty_deg_wt,ax_lfty_deg_wnt11,ax_lamba_t]
        linkyaxes!(ax_set...)
    end

    return fig

    # catch DomainError
    #     return "Numerical instability. Increase abstol/reltol"
    # end

end

# binned nodal profile version 

function plot_summary_newtimes_compare_binned!(fig,pv_list,prob,bin_size)

    legend_labels = [L"0.3t_{\text{WT}}",L"0.6t_{\text{WT}}",L"0.9t_{\text{WT}}",L"1.2t_{\text{WT}}",L"1.5t_{\text{WT}}",L"1.8t_{\text{WT}}"]

    ax_smad = []
    ax_smad_por = []
    ax_int_lefty = []
    ax_np_wt = []
    ax_np_wnt11 = []
    ax_n_prod_wt = []
    ax_n_prod_wnt11 = []
    ax_lfty_deg_wt = []
    ax_lfty_deg_wnt11 = []
    ax_lamba_t = []

    bin_x = [i for i in 0:bin_size:150-bin_size];

    bin_x = bin_x .+ bin_size ./ 2;

    data_dict = Dict(i => Dict() for i in 1:length(pv_list))

    for (n,pv) in enumerate(pv_list)

        p_tuple,p_cp_tuple,p_lm_tuple,p_ro_tuple = get_params_ro(pv)

        p_orig,_,_ = get_params(pv_orig)

        ax = Axis(fig[1,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Position of 20% max } c_N \text{ } (μm) ",ygridvisible = false,xgridvisible = false,title = "w/ cmax 0.2")
        ax_por = Axis(fig[1,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Average porosity, } \phi", yaxisposition = :right,ylabelcolor = :red,yticklabelcolor = :red,ygridvisible = false,xgridvisible = false,xticksvisible = false)

        (t_grid_alpha,dyn_alpha),(t_plot,(level_x_wt_rescaled,level_x_cp_rescaled,level_x_lm_rescaled )),(porosity_dyn,porosity_dyn_cp),c_level,(sol,sol_cp,sol_lm),(λhalf_t,λhalf_ro_t),(λhalf,λhalf_ro),λ_trange = get_alpha_xmax_lambda(pv,prob,0.2);

        hidexdecorations!(ax_por)

        rs_wt_max = maximum(level_x_wt_rescaled)

        lines!(ax,t_plot,level_x_wt_rescaled ./ rs_wt_max,color = :black,label = L"\text{WT}")
        lines!(ax,t_plot,level_x_cp_rescaled ./ rs_wt_max,color = :orange,label = L"\text{Wnt11}")
        lines!(ax,t_plot,level_x_lm_rescaled ./ rs_wt_max,color = :pink,label = L"\text{Lefty mutant}")

        lines!(ax_por,t_plot,porosity_dyn ,linestyle = :dash,color = :black,label = L"\text{ϕ}")
        lines!(ax_por,t_plot,porosity_dyn_cp ,linestyle = :dash,color = :orange,label = L"\text{ϕ}")

        axislegend(ax,position = :lt)

        ylims!(ax_por,0.,0.14)

        ax.xticks = (0:0.5:3.5,string.(0:0.5:3.5))

        push!(ax_smad,ax)
        push!(ax_smad_por,ax_por)

        orig_metrics = get_summary_metrics(pv,prob,data,alpha_data,0.2)

        t_plot_int = LinRange(0,3*orig_metrics[:wt_t0],t_plot_N)

        νN_int_cp,νN_int = get_integrated_lefty_prod_values(sol,sol_cp,t_plot_int,p_tuple,p_cp_tuple)

        # nodal_prod_cp,nodal_prod = get_integrated_nodal_prod_values(sol,sol_cp,t_plot_int)

        ax = Axis(fig[2,n], xlabel = L"\text{Rescaled time} t^{'}", ylabel= L"\int_0^L ν_L(N(t,x)) dx" )

        lines!(ax,LinRange(0,3,1000),νN_int_cp,linestyle = :dash,color = :grey, label = L"\text{Wnt11}")
        lines!(ax,LinRange(0,3,1000),νN_int,color = :grey, label = L"\text{WT}")
        
        axislegend(ax,position = :rt)

        ax = Axis(fig[3,n], xlabel = L"\text{Position}", ylabel= L"α(t)", title = "Beta = " * string(β) )

        for (n,d) in enumerate(dyn_alpha[1:end-1])
            lines!(ax,alpha_x,d, label = string(alpha_data_times_norm[n])* "* t_wt")
            scatter!(ax,alpha_x,alpha_data[:,n+1])
        end

        axislegend(ax,position = :rb)

        push!(ax_int_lefty,ax)

        alpha_data_times_norm_v = [0.,0.3,0.6,0.9,1.2,1.5,1.8]

        t_check =  alpha_data_times_norm_v[2:end] .* orig_metrics[:wt_t0]

        prob_finite = remake(prob,tspan = (0, alpha_data_times_norm_v[end] * orig_metrics[:wt_t0]))

        sol_profiles = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);
        sol_profiles_cp = solve(prob_finite, p = p_cp_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);

        dyn_N = [sol[:,1] for sol in sol_profiles.u[1:length(t_check)]]
        dyn_L = [sol[:,2] for sol in sol_profiles.u[1:length(t_check)]];
        dyn_α = [sol[:,4] for sol in sol_profiles.u[1:length(t_check)]];

        dyn_N_binned = [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in dyn_N]

        dyn_N_cp = [sol[:,1] for sol in sol_profiles_cp.u[1:length(t_check)]];
        dyn_L_cp = [sol[:,2] for sol in sol_profiles_cp.u[1:length(t_check)]];
        dyn_α_cp = [sol[:,4] for sol in sol_profiles_cp.u[1:length(t_check)]];

        dyn_N_cp_binned = [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in dyn_N_cp]

        lefty_prod_profiles_cp = [ν.(cN,σ.(p_cp_tuple[:σL0],ϕ0,ϕ.(α)),p_cp_tuple[:NL],p_cp_tuple[:mL]) for (cN,α) in zip(dyn_N_cp,dyn_α_cp)];
        lefty_prod_profiles = [ν.(cN,σ.(p_tuple[:σL0],ϕ0,ϕ.(α)),p_tuple[:NL],p_tuple[:mL]) for (cN,α) in zip(dyn_N,dyn_α)];

        nodal_prod_profiles_cp = [ν.(cN,σ.(p_cp_tuple[:σN0],ϕ0,ϕ.(α)),p_cp_tuple[:Na],p_cp_tuple[:mN]) for (cN,α) in zip(dyn_N_cp,dyn_α_cp)];
        nodal_prod_profiles = [ν.(cN,σ.(p_tuple[:σN0],ϕ0,ϕ.(α)),p_tuple[:Na],p_tuple[:mN]) for (cN,α) in zip(dyn_N,dyn_α)];

        nodal_prod_profiles_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in nodal_prod_profiles]
        nodal_prod_profiles_cp_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in nodal_prod_profiles_cp]

        lefty_prod_profiles_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in lefty_prod_profiles]
        lefty_prod_profiles_cp_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in lefty_prod_profiles_cp]

        # dyn_N = [sol[:,1] for sol in sol_profiles.u]
        # dyn_L = [sol[:,2] for sol in sol_profiles.u];

        # dyn_N_cp = [sol[:,1] for sol in sol_profiles_cp.u];
        # dyn_L_cp = [sol[:,2] for sol in sol_profiles_cp.u];

        ax1 = Axis(fig[4,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = "WT")

        for (i,N) in enumerate(dyn_N_binned)
            lines!(ax1,bin_x,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        push!(ax_np_wt,ax1)

        ax2 = Axis(fig[5,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = "Wnt11")

        for (i,N) in enumerate(dyn_N_cp_binned)
            lines!(ax2,bin_x,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax2,position = :rt)

        push!(ax_np_wnt11,ax2)

        # ax3 = Axis(fig[2,3], xlabel = L"\text{Position}", ylabel= L"\text{Lefty Concentration}", title = "WT")

        # for (i,N) in enumerate(dyn_L)
        #     lines!(ax3,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        # end

        # axislegend(ax3,position = :rt)

        # ax4 = Axis(fig[3,1], xlabel = L"\text{Position}", ylabel= L"\text{Lefty Concentration}", title = "Wnt11")

        # for (i,N) in enumerate(dyn_L_cp)
        #     lines!(ax4,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        # end

        # axislegend(ax4,position = :rt)

        axt1 = Axis(fig[11,n])
        axt2 = Axis(fig[12,n])

        hidedecorations!(axt1)
        hidedecorations!(axt2)

        text_pos = [Point(-0.5,i) for i in LinRange(-1,1,6)]

        for (n,p) in enumerate(p_names[1:6])

            if p_tuple[p] != p_orig[p]
                col = :red
            else
                col = :black
            end

            if p == :s0
                text!(axt1,text_pos[n], text = p_names_string[p] * " = " * string(round(2*p_tuple[p],digits = 8)),color = col)
            else
                text!(axt1,text_pos[n], text = p_names_string[p] * " = " * string(round(p_tuple[p],digits = 8)),color = col)
            end
        end

        # text!(axt1,-0.5,-1.5,text = "normal : " * string(round(λhalf_t,digits = 2)),color = :blue)

        ylims!(axt1,-1.8,1.8)
        ylims!(axt2,-1.8,1.8)
        xlims!(axt1,-1.5,1.5)
        xlims!(axt2,-1.5,1.5)

        text_pos = [Point(-0.5,i) for i in LinRange(-1,1,7)]

        for (n,p) in enumerate(p_names[7:end])

            if p_tuple[p] != p_orig[p]
                col = :red
            else
                col = :black
            end

            if p == :s0
                text!(axt2,text_pos[n], text = p_names_string[p] * " = " * string(round(2*p_tuple[p],digits = 8)),color = col)
            else
                text!(axt2,text_pos[n], text = p_names_string[p] * " = " * string(round(p_tuple[p],digits = 8)),color = col)
            end
        end

        # text!(axt2,-0.5,-1.5,text = "ro : " * string(round(λhalf_ro_t,digits = 2)),color = :green)

        # text!(axt2,-0.5,-1.7,text = "ratio ro/normal : " * string(round(λhalf_ro_t / λhalf_t,digits = 2)),color = :orange)

        #######

        ax1 = Axis(fig[6,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal production}", title = "WT")

        for (i,N) in enumerate(nodal_prod_profiles_binned)
            lines!(ax1,bin_x,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        push!(ax_n_prod_wt,ax1)

        ax2 = Axis(fig[7,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal production}", title = "Wnt11")

        for (i,N) in enumerate(nodal_prod_profiles_cp_binned)
            lines!(ax2,bin_x,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax2,position = :rt)

        push!(ax_n_prod_wnt11,ax2)

        # ax3 = Axis(fig[8,n], xlabel = L"\text{Position}", ylabel= L"\text{Lefty production}", title = "WT")

        # for (i,N) in enumerate(lefty_prod_profiles)
        #     lines!(ax3,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        # end

        # axislegend(ax3,position = :rt)

        # ax4 = Axis(fig[5,1], xlabel = L"\text{Position}", ylabel= L"\text{Lefty production}", title = "Wnt11")

        # for (i,N) in enumerate(lefty_prod_profiles_cp)
        #     lines!(ax4,N,label = string( alpha_data_times_norm_v[2:end][i]) * "* t_wt",colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        # end

        nodal_degr = [p_tuple[:kNL] .*(1 ./ (1 .+ (p_tuple[:LN] ./ L) .^ p_tuple[:mNL])) for L in dyn_L]
        nodal_degr_cp = [p_cp_tuple[:kNL].*(1 ./ (1 .+ (p_cp_tuple[:LN] ./ L) .^ p_cp_tuple[:mNL])) for L in dyn_L_cp]

        nodal_degr_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in nodal_degr]
        nodal_degr_cp_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in nodal_degr_cp]

        ax1 = Axis(fig[8,n], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = "WT")

        for (i,N) in enumerate(nodal_degr)
            lines!(ax1,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        axislegend(ax1,position = :rt)

        push!(ax_lfty_deg_wt,ax1)

        ax2 = Axis(fig[9,n], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = "Wnt11")

        for (i,N) in enumerate(nodal_degr_cp)
            lines!(ax2,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
        end

        push!(ax_lfty_deg_wnt11,ax2)

        ts = max(λhalf_t,λhalf_ro_t)
    
        prob_finite = remake(prob,tspan = (0, 3*ts))

        λ_trange = LinRange(0.,3*ts,10000)

        sol1  = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);
        sol_ro1  = solve(prob_finite, p = p_ro_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);

        N0t = [sol1(t)[1,1] for t in λ_trange]    
        
        λhalf,λhalf_max_t = get_lambda_half(sol1,λ_trange)
        λhalf_ro,λhalf_max_t_ro = get_lambda_half(sol_ro1,λ_trange)

        ax = Axis(fig[10,n], xlabel = L"\text{time (mins)}", ylabel= L"\lambda_{1/2} \quad (\mu m)",ygridvisible = false,xgridvisible = false)

        # ax_N = Axis(fig[10,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"N(x=0,t)", yaxisposition = :right,ylabelcolor = :orange,yticklabelcolor = :orange,ygridvisible = false,xgridvisible = false,xticksvisible = false)

        # hidexdecorations!(ax_N)

        # lines!(ax_N,λ_trange,N0t,color = :orange, label = "WT")

        text!(axt1,-0.5,-1.5,text = "normal lambda_1/2 : " * string(round(maximum(λhalf),digits = 2)),color = :blue)
        text!(axt2,-0.5,-1.5,text = "ro lambda_1/2: " * string(round(maximum(λhalf_ro),digits = 2)),color = :green)

        text!(axt2,-0.5,-1.7,text = "ratio: " * string(round(maximum(λhalf_ro) / maximum(λhalf),digits = 2)),color = :orange)

        vlines!(ax,λhalf_max_t,color = :blue,linestyle = :dash)
        vlines!(ax,λhalf_max_t_ro,color = :green,linestyle = :dash)

        lines!(ax,λ_trange,λhalf,color = :blue, label = "WT")

        lines!(ax,λ_trange,λhalf_ro,color = :green, label = "relay off (WT)")

        push!(ax_lamba_t,ax)

    end

    for ax_set in  [ax_smad,ax_smad_por,ax_int_lefty,ax_np_wt,ax_np_wnt11,ax_n_prod_wt,ax_n_prod_wnt11,ax_lfty_deg_wt,ax_lfty_deg_wnt11,ax_lamba_t]
        linkyaxes!(ax_set...)
    end

    return fig

    # catch DomainError
    #     return "Numerical instability. Increase abstol/reltol"
    # end

end

# binned nodal profile version with reduced set of plots for SM 

function plot_summary_newtimes_compare_binned_reduced!(fig,pv_list,prob,bin_size)

    # legend_labels = [L"0.3t_{\text{WT}}",L"0.6t_{\text{WT}}",L"0.9t_{\text{WT}}",L"1.2t_{\text{WT}}",L"1.5t_{\text{WT}}",L"1.8t_{\text{WT}}"]

    legend_labels = [L"0.5t_{\text{WT}}",L"t_{\text{WT}}",L"1.5t_{\text{WT}}",L"2t_{\text{WT}}",L"2.5t_{\text{WT}}",L"3t_{\text{WT}}"]

    ax_smad = []
    ax_smad_por = []
    ax_int_lefty = []
    ax_np_wt = []
    ax_np_wnt11 = []
    ax_n_prod_wt = []
    ax_n_prod_wnt11 = []
    ax_lfty_deg_wt = []
    ax_lfty_deg_wnt11 = []
    ax_lamba_t = []

    bin_x = [i for i in 0:bin_size:150-bin_size];

    bin_x = bin_x .+ bin_size ./ 2;

    data_dict = Dict(i => Dict() for i in 1:length(pv_list))

    for (n,pv) in enumerate(pv_list)

        p_tuple,p_cp_tuple,p_lm_tuple,p_ro_tuple = get_params_ro(pv)

        p_orig,_,_ = get_params(pv_orig)

        ax = Axis(fig[1,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Position of 20% max } c_N \text{ } (μm) ",ygridvisible = false,xgridvisible = false)
        ax_por = Axis(fig[1,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Average porosity, } \phi", yaxisposition = :right,ylabelcolor = :red,yticklabelcolor = :red,ygridvisible = false,xgridvisible = false,xticksvisible = false)

        (t_grid_alpha,dyn_alpha),(t_plot,(level_x_wt_rescaled,level_x_cp_rescaled,level_x_lm_rescaled )),(porosity_dyn,porosity_dyn_cp),c_level,(sol,sol_cp,sol_lm),(λhalf_t,λhalf_ro_t),(λhalf,λhalf_ro),λ_trange = get_alpha_xmax_lambda(pv,prob,0.2);

        hidexdecorations!(ax_por)

        rs_wt_max = maximum(level_x_wt_rescaled)

        lines!(ax,t_plot,level_x_wt_rescaled ./ rs_wt_max,color = :black,label = L"\text{with feedback}")
        lines!(ax,t_plot,level_x_cp_rescaled ./ rs_wt_max,color = :orange,label = L"\text{without feedback}")
        # lines!(ax,t_plot,level_x_lm_rescaled ./ rs_wt_max,color = :pink,label = L"\text{Lefty mutant}")

        lines!(ax_por,t_plot,porosity_dyn ,linestyle = :dash,color = :black,label = L"\text{ϕ}")
        lines!(ax_por,t_plot,porosity_dyn_cp ,linestyle = :dash,color = :orange,label = L"\text{ϕ}")

        data_dict[n]["SMAD_data"] = Dict("WT"=>level_x_wt_rescaled ./ rs_wt_max ,"Wnt11" => level_x_cp_rescaled ./ rs_wt_max,"Lefty"=> level_x_lm_rescaled ./rs_wt_max, "porosity_wt" => porosity_dyn,"porosity_wnt11" => porosity_dyn_cp, "t_plot" => t_plot)

        axislegend(ax,position = :lt)

        ylims!(ax_por,0.,0.14)

        ax.xticks = (0:0.5:3.5,string.(0:0.5:3.5))

        push!(ax_smad,ax)
        push!(ax_smad_por,ax_por)

        orig_metrics = get_summary_metrics(pv,prob,data,alpha_data,0.2)

        t_plot_int = LinRange(0,3*orig_metrics[:wt_t0],t_plot_N)

        νN_int_cp,νN_int = get_integrated_lefty_prod_values(sol,sol_cp,t_plot_int,p_tuple,p_cp_tuple)

        ax = Axis(fig[2,n], xlabel = L"\text{Rescaled time, } t", ylabel= L"\int_0^L ν_L(N(t,x)) dx" ,ygridvisible = false,xgridvisible = false)

        lines!(ax,LinRange(0,3,1000),νN_int_cp,linestyle = :dash,color = :grey, label = L"\text{without feedback}")
        lines!(ax,LinRange(0,3,1000),νN_int,color = :grey, label = L"\text{with feedback}")

        data_dict[n]["IntLefty"] = Dict("WT"=>νN_int,"Wnt11" => νN_int_cp,"t_plot" => LinRange(0,3,1000))
        
        axislegend(ax,position = :rt)

        push!(ax_int_lefty,ax)

        alpha_data_times_norm_v = [0.,0.5,1,1.5,2,2.5,3]

        t_check =  alpha_data_times_norm_v[2:end] .* orig_metrics[:wt_t0]

        prob_finite = remake(prob,tspan = (0, alpha_data_times_norm_v[end] * orig_metrics[:wt_t0]))

        sol_profiles = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);
        sol_profiles_cp = solve(prob_finite, p = p_cp_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = t_check);

        dyn_N = [sol[:,1] for sol in sol_profiles.u[1:length(t_check)]]
        dyn_L = [sol[:,2] for sol in sol_profiles.u[1:length(t_check)]];

        dyn_N_binned = [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in dyn_N]

        max_N = maximum([maximum(v) for v in dyn_N_binned])

        dyn_N_cp = [sol[:,1] for sol in sol_profiles_cp.u[1:length(t_check)]];
        dyn_L_cp = [sol[:,2] for sol in sol_profiles_cp.u[1:length(t_check)]];

        dyn_N_cp_binned = [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in dyn_N_cp]

        ax1 = Axis(fig[3,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = L"\text{with feedback}",ygridvisible = false,xgridvisible = false)

        nodal_profile_dict = Dict()

        for (i,N) in enumerate(dyn_N_binned)
            lines!(ax1,bin_x,N ./ max_N ,label =legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
            nodal_profile_dict[legend_labels[i]] = N  ./ max_N
        end

        Colorbar(fig,colormap = :viridis,colorrange = (1,length(t_check)),bbox=ax1.scene.px_area,halign = :right, width = 10, height = 40,valign = :top,ticksvisible = false,ticklabelsvisible  = false)

        data_dict[n]["NodalProfiles"] = nodal_profile_dict

        # axislegend(ax1,position = :rt)

        push!(ax_np_wt,ax1)

        ax2 = Axis(fig[4,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = L"\text{without feedback}",ygridvisible = false,xgridvisible = false)

        nodal_profile_dict_cp = Dict()

        for (i,N) in enumerate(dyn_N_cp_binned)
            lines!(ax2,bin_x,N ./ max_N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
            nodal_profile_dict_cp[legend_labels[i]] = N  ./ max_N
        end

        Colorbar(fig,colormap = :viridis,colorrange = (1,length(t_check)),bbox=ax2.scene.px_area,halign = :right, width = 10, height = 40,valign = :top,ticksvisible = false,ticklabelsvisible  = false)

        data_dict[n]["NodalProfiles_wnt11"] = nodal_profile_dict_cp

        # axislegend(ax2,position = :rt)

        push!(ax_np_wnt11,ax2)

        nodal_degr = [p_tuple[:kNL] .*(1 ./ (1 .+ (p_tuple[:LN] ./ L) .^ p_tuple[:mNL])) for L in dyn_L]
        nodal_degr_cp = [p_cp_tuple[:kNL].*(1 ./ (1 .+ (p_cp_tuple[:LN] ./ L) .^ p_cp_tuple[:mNL])) for L in dyn_L_cp]

        # nodal_degr_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in nodal_degr]
        # nodal_degr_cp_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in nodal_degr_cp]

        ax1 = Axis(fig[5,n], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = L"\text{with feedback}",ygridvisible = false,xgridvisible = false)

        nodal_degr_dict = Dict()

        for (i,N) in enumerate(nodal_degr)
            lines!(ax1,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
            nodal_degr_dict[legend_labels[i]] = N
        end

        Colorbar(fig,colormap = :viridis,colorrange = (1,length(t_check)),bbox=ax1.scene.px_area,halign = :right, width = 10, height = 40,valign = :top,ticksvisible = false,ticklabelsvisible  = false)

        data_dict[n]["LeftyDegradationN"] = nodal_degr_dict

        # axislegend(ax1,position = :rt)

        push!(ax_lfty_deg_wt,ax1)

        ax2 = Axis(fig[6,n], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = L"\text{without feedback}",ygridvisible = false,xgridvisible = false)

        nodal_degr_dict_cp = Dict()

        for (i,N) in enumerate(nodal_degr_cp)
            lines!(ax2,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
            nodal_degr_dict_cp[legend_labels[i]] = N
        end

        Colorbar(fig,colormap = :viridis,colorrange = (1,length(t_check)),bbox=ax2.scene.px_area,halign = :right, width = 10, height = 40,valign = :top,ticksvisible = false,ticklabelsvisible  = false)

        data_dict[n]["LeftyDegradationN_wnt11"] = nodal_degr_dict_cp

        push!(ax_lfty_deg_wnt11,ax2)

        ts = max(λhalf_t,λhalf_ro_t)
    
        prob_finite = remake(prob,tspan = (0, 3*ts))

        λ_trange = LinRange(0.,3*ts,10000)

        sol1  = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);
        sol_ro1  = solve(prob_finite, p = p_ro_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);

        # N0t = [sol1(t)[1,1] for t in λ_trange]    
        
        λhalf,λhalf_max_t = get_lambda_half(sol1,λ_trange)
        λhalf_ro,λhalf_max_t_ro = get_lambda_half(sol_ro1,λ_trange)

        ax = Axis(fig[7,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\lambda_{1/2} \quad (\mu m)",ygridvisible = false,xgridvisible = false)

        vlines!(ax,λhalf_max_t / orig_metrics[:wt_t0],color = :blue,linestyle = :dash)
        vlines!(ax,λhalf_max_t_ro / orig_metrics[:wt_t0],color = :green,linestyle = :dash)
    
        lines!(ax,λ_trange ./ orig_metrics[:wt_t0] ,λhalf,color = :blue, label = L"\text{relay on}")

        lines!(ax,λ_trange  ./ orig_metrics[:wt_t0] ,λhalf_ro,color = :green, label = L"\text{relay off}")

        data_dict[n]["Lambda1/2"] = Dict("RelayOff"=>λhalf_ro, "RelayOn"=>λhalf)

        axislegend(ax,position = :rt)

        push!(ax_lamba_t,ax)

    end

    for ax_set in  [ax_smad,ax_smad_por,ax_int_lefty,ax_np_wt,ax_np_wnt11,ax_lfty_deg_wt,ax_lfty_deg_wnt11,ax_lamba_t]
        linkyaxes!(ax_set...)
    end

    return fig, data_dict

    # catch DomainError
    #     return "Numerical instability. Increase abstol/reltol"
    # end

end

# binned nodal profile version with reduced set of plots for SM, but uses existing data to speed up (pass data_dict)

function plot_summary_newtimes_compare_binned_reduced!(fig,pv_list,prob,bin_size, data_dict)

    # legend_labels = [L"0.3t_{\text{WT}}",L"0.6t_{\text{WT}}",L"0.9t_{\text{WT}}",L"1.2t_{\text{WT}}",L"1.5t_{\text{WT}}",L"1.8t_{\text{WT}}"]

    legend_labels = [L"0.5t_{\text{WT}}",L"t_{\text{WT}}",L"1.5t_{\text{WT}}",L"2t_{\text{WT}}",L"2.5t_{\text{WT}}",L"3t_{\text{WT}}"]

    ax_smad = []
    ax_smad_por = []
    ax_int_lefty = []
    ax_np_wt = []
    ax_np_wnt11 = []
    ax_n_prod_wt = []
    ax_n_prod_wnt11 = []
    ax_lfty_deg_wt = []
    ax_lfty_deg_wnt11 = []
    ax_lamba_t = []

    bin_x = [i for i in 0:bin_size:150-bin_size];

    bin_x = bin_x .+ bin_size ./ 2;

    data_dict = Dict(i => Dict() for i in 1:length(pv_list))

    for (n,pv) in enumerate(pv_list)

        p_tuple,p_cp_tuple,p_lm_tuple,p_ro_tuple = get_params_ro(pv)

        p_orig,_,_ = get_params(pv_orig)

        ax = Axis(fig[1,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Position of 20% max } c_N \text{ } (μm) ",ygridvisible = false,xgridvisible = false,title = "w/ cmax 0.2")
        ax_por = Axis(fig[1,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\text{Average porosity, } \phi", yaxisposition = :right,ylabelcolor = :red,yticklabelcolor = :red,ygridvisible = false,xgridvisible = false,xticksvisible = false)

        hidexdecorations!(ax_por)

        lines!(ax,data_dict[n]["SMAD_data"]["t_plot"],data_dict[n]["SMAD_data"]["WT"],color = :black,label = L"\text{WT}")
        lines!(ax,data_dict[n]["SMAD_data"]["t_plot"],data_dict[n]["SMAD_data"]["Wnt11"],color = :orange,label = L"\text{Wnt11}")
        lines!(ax,data_dict[n]["SMAD_data"]["t_plot"],data_dict[n]["SMAD_data"]["Lefty"],color = :pink,label = L"\text{Lefty mutant}")

        lines!(ax_por,data_dict[n]["SMAD_data"]["t_plot"],data_dict[n]["SMAD_data"]["porosity_wt"],linestyle = :dash,color = :black,label = L"\text{ϕ}")
        # lines!(ax_por,t_plot,data_dict[n]["SMAD_data"]["porosity_wnt11"] ,linestyle = :dash,color = :orange,label = L"\text{ϕ}")

        data_dict[n]["SMAD_data"] = Dict("WT"=>level_x_wt_rescaled ./ rs_wt_max ,"Wnt11" => level_x_cp_rescaled ./ rs_wt_max,"Lefty"=> level_x_lm_rescaled ./rs_wt_max, "porosity_wt" => porosity_dyn, "t_plot" => t_plot)

        axislegend(ax,position = :lt)

        ylims!(ax_por,0.,0.14)

        ax.xticks = (0:0.5:3.5,string.(0:0.5:3.5))

        push!(ax_smad,ax)
        push!(ax_smad_por,ax_por)

         # plot nodal profiles

        t_plot_int = LinRange(0,3*orig_metrics[:wt_t0],t_plot_N)

        νN_int_cp,νN_int = get_integrated_lefty_prod_values(sol,sol_cp,t_plot_int,p_tuple,p_cp_tuple)

        ax = Axis(fig[2,n], xlabel = L"\text{Rescaled time} t^{'}", ylabel= L"\int_0^L ν_L(N(t,x)) dx" )

        lines!(ax,LinRange(0,3,1000),νN_int_cp,linestyle = :dash,color = :grey, label = L"\text{Without feedback}")
        lines!(ax,LinRange(0,3,1000),νN_int,color = :grey, label = L"\text{With feedback}")

        data_dict[n]["IntLefty"] = Dict("WT"=>νN_int,"Wnt11" => νN_int_cp,"t_plot" => LinRange(0,3,1000))
        
        axislegend(ax,position = :rt)

        push!(ax_int_lefty,ax)

        alpha_data_times_norm_v = [0.,0.5,1,1.5,2,2.5,3]

        t_check =  alpha_data_times_norm_v[2:end] .* orig_metrics[:wt_t0]

        # plot nodal profiles

        ax1 = Axis(fig[3,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = L"\text{With feedback}")

        nodal_profile_dict = Dict()

        for (i,N) in enumerate(dyn_N_binned)
            lines!(ax1,bin_x,N ./ max_N ,label =legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
            nodal_profile_dict[legend_labels[i]] = N  ./ max_N
        end

        Colorbar(fig,colormap = :viridis,colorrange = (1,length(t_check)),bbox=ax1.scene.px_area,halign = :right, width = 10, height = 40,valign = :top,ticksvisible = false,ticklabelsvisible  = false)

        data_dict[n]["NodalProfiles"] = nodal_profile_dict

        axislegend(ax1,position = :rt)

        push!(ax_np_wt,ax1)

        ax2 = Axis(fig[4,n], xlabel = L"\text{Position}", ylabel= L"\text{Nodal Concentration}", title = L"\text{Without feedback}")

        nodal_profile_dict_cp = Dict()

        for (i,N) in enumerate(dyn_N_cp_binned)
            lines!(ax2,bin_x,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
            nodal_profile_dict_cp[legend_labels[i]] = N  ./ max_N
        end

        Colorbar(fig,colormap = :viridis,colorrange = (1,length(t_check)),bbox=ax2.scene.px_area,halign = :right, width = 10, height = 40,valign = :top,ticksvisible = false,ticklabelsvisible  = false)

        data_dict[n]["NodalProfiles_wnt11"] = nodal_profile_dict_cp

        axislegend(ax2,position = :rt)

        push!(ax_np_wnt11,ax2)

        nodal_degr = [p_tuple[:kNL] .*(1 ./ (1 .+ (p_tuple[:LN] ./ L) .^ p_tuple[:mNL])) for L in dyn_L]
        nodal_degr_cp = [p_cp_tuple[:kNL].*(1 ./ (1 .+ (p_cp_tuple[:LN] ./ L) .^ p_cp_tuple[:mNL])) for L in dyn_L_cp]

        # nodal_degr_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in nodal_degr]
        # nodal_degr_cp_binned =  [[mean(N[n+1:n+bin_size]) for n in 0:bin_size:150-bin_size] for N in nodal_degr_cp]

        ax1 = Axis(fig[5,n], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = L"\text{With feedback}")

        nodal_degr_dict = Dict()

        for (i,N) in enumerate(nodal_degr)
            lines!(ax1,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
            nodal_degr_dict[legend_labels[i]] = N
        end

        Colorbar(fig,colormap = :viridis,colorrange = (1,length(t_check)),bbox=ax1.scene.px_area,halign = :right, width = 10, height = 40,valign = :top,ticksvisible = false,ticklabelsvisible  = false)

        data_dict[n]["LeftyDegradationN"] = nodal_degr_dict

        axislegend(ax1,position = :rt)

        push!(ax_lfty_deg_wt,ax1)

        ax2 = Axis(fig[6,n], xlabel = L"\text{Position}", ylabel= L"k_{NL} (1 / (1 + (L_N / L)^{m_{NL}})", title = L"\text{Without feedback}")

        nodal_degr_dict_cp = Dict()

        for (i,N) in enumerate(nodal_degr_cp)
            lines!(ax2,N,label = legend_labels[i],colormap = :viridis, color = i, colorrange = (1,length(t_check)))
            nodal_degr_dict_cp[legend_labels[i]] = N
        end

        Colorbar(fig,colormap = :viridis,colorrange = (1,length(t_check)),bbox=ax2.scene.px_area,halign = :right, width = 10, height = 40,valign = :top,ticksvisible = false,ticklabelsvisible  = false)

        data_dict[n]["LeftyDegradationN_wnt11"] = nodal_degr_dict_cp

        push!(ax_lfty_deg_wnt11,ax2)

        ts = max(λhalf_t,λhalf_ro_t)
    
        prob_finite = remake(prob,tspan = (0, 3*ts))

        λ_trange = LinRange(0.,3*ts,10000)

        sol1  = solve(prob_finite, p = p_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);
        sol_ro1  = solve(prob_finite, p = p_ro_tuple, default_solver,abstol = de_abstol,reltol = de_reltol, maxiters = 1e6,isoutofdomain = (u,p,t) -> any(x->x<0, u), saveat = λ_trange);

        # N0t = [sol1(t)[1,1] for t in λ_trange]    
        
        λhalf,λhalf_max_t = get_lambda_half(sol1,λ_trange)
        λhalf_ro,λhalf_max_t_ro = get_lambda_half(sol_ro1,λ_trange)

        ax = Axis(fig[7,n], xlabel = L"\text{Rescaled time, t}", ylabel= L"\lambda_{1/2} \quad (\mu m)",ygridvisible = false,xgridvisible = false)

        vlines!(ax,λhalf_max_t / orig_metrics[:wt_t0],color = :blue,linestyle = :dash)
        vlines!(ax,λhalf_max_t_ro / orig_metrics[:wt_t0],color = :green,linestyle = :dash)
    
        lines!(ax,λ_trange ./ orig_metrics[:wt_t0] ,λhalf,color = :blue, label = L"\text{relay on}")

        lines!(ax,λ_trange  ./ orig_metrics[:wt_t0] ,λhalf_ro,color = :green, label = L"\text{relay off}")

        data_dict[n]["Lambda1/2"] = Dict("RelayOff"=>λhalf_ro, "RelayOn"=>λhalf)

        axislegend(ax,position = :rt)

        push!(ax_lamba_t,ax)

    end

    for ax_set in  [ax_smad,ax_smad_por,ax_int_lefty,ax_np_wt,ax_np_wnt11,ax_lfty_deg_wt,ax_lfty_deg_wnt11,ax_lamba_t]
        linkyaxes!(ax_set...)
    end

    return fig, data_dict

    # catch DomainError
    #     return "Numerical instability. Increase abstol/reltol"
    # end

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