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
