using LaTeXStrings
using Juqbox
using Plots
using DelimitedFiles
using StatsPlots
using Plots.PlotMeasures

# Function to loop over ϵ values and create a plot of the objective function for pertrubed Hamiltonians
function ep_plot(pcof::Vector{Float64}, params:: Juqbox.objparams, wa::Working_Arrays, ep_vals::AbstractArray)
    results = zeros(length(ep_vals),4)

    for i =1:length(ep_vals)
        epsh = ep_vals[i]
        
        # Additive noise
        for j = 2:size(params.Hconst,2)
            params.Hconst[j,j] += 0.01*epsh*(10.0^(j-2))
        end

        obj, _, _, secondaryobjf, traceinfid = Juqbox.traceobjgrad(pcof,params,wa,false, true)
        results[i,1] = obj
        results[i,2] = secondaryobjf
        results[i,3] = traceinfid
        results[i,4] = traceinfid+secondaryobjf
        
        # Reset
        for j = 2:size(params.Hconst,2)
            params.Hconst[j,j] -= 0.01*epsh*(10.0^(j-2))
        end

    end

    pl1 = Plots.plot(ep_vals,results[:,1],yaxis=:log,xlabel = L"\epsilon",ylabel = "Objective Function")
    return results,pl1
end


function evalctrl_no_carrier(params::objparams, pcof0:: Array{Float64, 1}, jFunc:: Int64, D1::Int64) 
    
    # Evaluate the ctrl functions on this grid in time
    nplot = round(Int64, params.T*32)
    # is this resolution sufficient for the lab frame ctrl functions so we can get meaningful FFTs?
    td = collect(range(0, stop = params.T, length = nplot+1))

    nfreq = length(params.Cfreq)
    offset = (jFunc-1)*2*D1
    nCoeff = 2*D1
    pcof = copy(pcof0[offset+1:offset+nCoeff])

    if (params.use_bcarrier)
        # B-splines with carrier waves (zero out the frequency to plot just the splines)
        splinepar = Juqbox.bcparams(params.T, D1, params.Ncoupled, params.Nunc, zeros(1,1), pcof)
    else
        # regular B-splines
        splinepar = Juqbox.splineparams(params.T, D1, 2*(params.Ncoupled + params.Nunc), pcof)
    end

    # define inline function to enable vectorization over t
    controlplot(t, splinefunc) = Juqbox.controlfunc(t, splinepar, splinefunc)

    fact = 1.0/(2*pi) # conversion factor to GHz
    fact = 1.0 # conversion factor to rad/ns
    pj = fact.*controlplot.(td, 0)
    qj = fact.*controlplot.(td, 1)
    return pj, qj
    
end


# Maximum shift in Hamiltonian (in rad*GHz)
ep_max = 2*pi*2e-2

# For plotting purposes 
max_ep_sweep = 2*pi*3e-2
len = 1001
ep_vals = range(-max_ep_sweep,stop=max_ep_sweep,length=len)
scalefactor = 1000/(2*pi)
unitStr = "MHz"

fnt = Plots.font("Computer Modern", 12)
lfnt = Plots.font("Computer Modern", 12)
Plots.default(fontfamily="Computer Modern", titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=lfnt, linewidth=1.5, size=(750, 400), dpi=300)

freshOptim = false
#### Usual optimization, control & data stored in pcof_old & data 
nquad = 1
switch_loss = 1
if(freshOptim)
    include("swap-02-risk-neutral.jl")
    pcof_NF = Juqbox.run_optimizer(prob, pcof0)
    results,pl1 = ep_plot(pcof_NF, params, wa, ep_vals)
    data_NF = zeros(size(results,1),1)
    data_NF[:,1] = results[:,3]
    writedlm("NF_control.dat", pcof_NF)
    writedlm("NF_optim.dat", data_NF)
else
    pcof_NF = vec(readdlm("NF_control.dat"))
    data_NF = readdlm("NF_optim.dat")
end

freshOptim = false
#### CVaR with beta = 0.9
nquad = 100
switch_loss = 4
hyper = 0.9
if(freshOptim)
    include("swap-02-risk-neutral.jl")
    pcof_cvar4 = Juqbox.run_optimizer(prob, pcof0)
    results,pl1 = ep_plot(pcof_cvar4, params, wa, ep_vals)
    data_cvar4 = zeros(size(results,1),1)
    data_cvar4[:,1] = results[:,3]
    writedlm("cvar_control_quad100_beta90_itr300.dat", pcof_cvar4)
    writedlm("cvar_optim_quad100_beta90_itr300.dat", data_cvar4)
else
    pcof_cvar4 = vec(readdlm("cvar_control_quad100_beta90_itr300.dat"))
    data_cvar4 = readdlm("cvar_optim_quad100_beta90_itr300.dat")
end

freshOptim = false
#### CVaR with beta = 0.9
nquad = 500
switch_loss = 4
hyper = 0.9
if(freshOptim)
    include("swap-02-risk-neutral.jl")
    pcof_cvar = Juqbox.run_optimizer(prob, pcof0)
    results,pl1 = ep_plot(pcof_cvar, params, wa, ep_vals)
    data_cvar = zeros(size(results,1),1)
    data_cvar[:,1] = results[:,3]
    writedlm("cvar_control_quad500_beta90_itr300.dat", pcof_cvar)
    writedlm("cvar_optim_quad500_beta90_itr300.dat", data_cvar)
else
    pcof_cvar = vec(readdlm("cvar_control_quad500_beta90_itr300.dat"))
    data_cvar = readdlm("cvar_optim_quad500_beta90_itr300.dat")
end
#=
freshOptim = true
#### CVaR with beta = 0.9
nquad = 35
switch_loss = 4
hyper = 0.9
if(freshOptim)
    include("swap-02-risk-neutral.jl")
    pcof_cvar2 = Juqbox.run_optimizer(prob, pcof0)
    results,pl1 = ep_plot(pcof_cvar2, params, wa, ep_vals)
    data_cvar2 = zeros(size(results,1),4)
    data_cvar2[:,1] = ep_vals
    for j = 1:3
        data_cvar2[:,j+1] = results[:,j]
    end
    writedlm("cvar_control_35rand_beta90_itr300.dat", pcof_cvar2)
    writedlm("cvar_optim_35rand_beta90_itr300.dat", data_cvar2)
else
    pcof_cvar2 = vec(readdlm("cvar_control_51quad_beta90.dat"))
    data_cvar2 = readdlm("cvar_optim_51quad_beta90.dat")
end
=#
freshOptim = true
#### CVaR with beta = 0.9
nquad = 1000
switch_loss = 4
hyper = 0.9
if(freshOptim)
    include("swap-02-risk-neutral.jl")
    pcof_cvar3 = Juqbox.run_optimizer(prob, pcof0)
    results,pl1 = ep_plot(pcof_cvar3, params, wa, ep_vals)
    data_cvar3 = zeros(size(results,1),1)
    data_cvar3[:,1] = results[:,3]
    writedlm("cvar_control_quad1000_beta90_itr300.dat", pcof_cvar3)
    writedlm("cvar_optim_quad1000_beta90_itr300.dat", data_cvar3)
else
    pcof_cvar3 = vec(readdlm("cvar_control_quad1000_beta90_itr300.dat"))
    data_cvar3 = readdlm("cvar_optim_quad1000_beta90_itr300.dat")
end


freshOptim = false
#### RN
switch_loss = 2
nquad = 9
if(freshOptim)
    include("swap-02-risk-neutral.jl")
    pcof_RN = Juqbox.run_optimizer(prob, pcof0)
    results,pl1 = ep_plot(pcof_RN, params, wa, ep_vals)
    data_RN = zeros(size(results,1),4)
    data_RN[:,1] = results[:,3]
    writedlm("RN_control_quad9.dat", pcof_RN)
    writedlm("RN_optim_quad9.dat", data_RN)
else
    pcof_RN = vec(readdlm("RN_control_51quad.dat"))
    data_RN = readdlm("RN_optim_51quad.dat")
end


#### comparison between loss functions (NF, cvar)
# Plot objective functions for various Hamiltonian perturbations
plc = Plots.plot(scalefactor.*ep_vals,data_NF[:,1],yaxis=:log,xlabel = "Hamiltonian Perturbation 𝜖 [MHz]", ylabel = "Loss Function", lab="NF ctrl",
                  title = "Robustness against noise, beta=0.9, Monte Carlo sampling", legend= :outerright, color = "green", legendfontsize=10, margin = 20px)
Plots.plot!(plc,scalefactor.*ep_vals,data_cvar4[:,1],yaxis=:log,xlabel = "𝜖: Hamiltonian Perturbation [MHz]", ylabel = "Loss Function", lab="CVaR ctrl, Nsample=100", color = "orange")
Plots.plot!(plc,scalefactor.*ep_vals,data_cvar[:,1],yaxis=:log,xlabel = "𝜖: Hamiltonian Perturbation [MHz]", ylabel = "Loss Function", lab="CVaR ctrl, Nsample=500", color = "red")
# Plots.plot!(plc,scalefactor.*ep_vals,data_cvar2[:,4],yaxis=:log,xlabel = "𝜖: Hamiltonian Perturbation [MHz]", ylabel = "Loss Function", lab="CVaR ctrl, nquad=35", color = "blue")
Plots.plot!(plc,scalefactor.*ep_vals,data_cvar3[:,1],yaxis=:log,xlabel = "𝜖: Hamiltonian Perturbation [MHz]", ylabel = "Loss Function", lab="CVaR ctrl, Nsample=1000", color = "purple")
# Plots.plot!(plc,scalefactor.*ep_vals,data_RN[:,4],yaxis=:log,xlabel = "𝜖: Hamiltonian Perturbation [MHz]", ylabel = "Loss Function", lab="RN ctrl", color = "blue")
# vspan!([-10, 10], linecolor = :yellow, fillcolor = :yellow, fillalpha = 0.2, lab="noise", legendfontsize=16)

# Save plots to file
Plots.savefig(plc,  "loss_comparison_MC_eps4_itr300.png")
