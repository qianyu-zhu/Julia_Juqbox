using LaTeXStrings
using Juqbox
using Plots
using DelimitedFiles
using StatsPlots
using Plots.PlotMeasures
using Statistics
using KernelDensity

# Function to loop over ϵ values and create a plot of the objective function for perturbed Hamiltonians
function ep_plot(pcof::Vector{Float64}, params:: Juqbox.objparams, wa::Working_Arrays, ep_vals::AbstractArray)
    results = zeros(size(ep_vals, 2),4)

    for i =1:size(ep_vals, 2)
        ep = ep_vals[:,i]
        
        # Additive noise
        for j = 1:size(params.Hconst,2)
            # params.Hconst[j,j] += H0_old[j,j] + 0.01*ep*(10.0^(j-2))
            params.Hconst[j,j] += 0.0001*ep[j]*(10.0^(j))
        end
        # params.Hconst[1,2] += 0.001*ep[5]
        # params.Hconst[2,1] += 0.001*ep[5]
        # params.Hconst[2,3] += 0.001*ep[6]
        # params.Hconst[3,2] += 0.001*ep[6]
        # params.Hconst[3,4] += 0.001*ep[7]
        # params.Hconst[4,3] += 0.001*ep[7]

        obj, _, _, secondaryobjf, traceinfid = Juqbox.traceobjgrad(pcof,params,wa,false, true)
        results[i,1] = obj
        results[i,2] = secondaryobjf
        results[i,3] = traceinfid
        results[i,4] = traceinfid+secondaryobjf
        
        # Reset
        for j = 1:size(params.Hconst,2)
            # params.Hconst[j,j] += H0_old[j,j] + 0.01*ep*(10.0^(j-2))
            params.Hconst[j,j] -= 0.0001*ep[j]*(10.0^(j))
        end
        # params.Hconst[1,2] -= 0.001*ep[5]
        # params.Hconst[2,1] -= 0.001*ep[5]
        # params.Hconst[2,3] -= 0.001*ep[6]
        # params.Hconst[3,2] -= 0.001*ep[6]
        # params.Hconst[3,4] -= 0.001*ep[7]
        # params.Hconst[4,3] -= 0.001*ep[7]

    end

    # pl1 = Plots.plot(ep_vals,results[:,1],yaxis=:log,xlabel = L"\epsilon",ylabel = "Objective Function")
    return results  # ,pl1
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
len = 10000 # number of samples
dim_noise = 4
ep_vals = (rand(Float64, (dim_noise, len)).- 0.5) .* (2*max_ep_sweep)
writedlm("ep_vals_10000.dat", ep_vals)
# ep_vals = readdlm("ep_vals.dat")


scalefactor = 1000/(2*pi)
unitStr = "MHz"

fnt = Plots.font("Computer Modern", 12)
lfnt = Plots.font("Computer Modern", 12)
Plots.default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=lfnt, linewidth=1.5, size=(650, 400), margin = 10px, dpi=300)
trajectory = 3
ctrl_var = 48

freshOptim = false
#### Noise-Free optimization, control & data stored in pcof_NF & data_NF
nquad = 1
switch_loss = 1

if(freshOptim)
    data_NF_his = zeros(size(ep_vals, 2), trajectory)
    ctrl_NF_his = zeros(ctrl_var, trajectory)
    for i = 1:trajectory
        include("swap-02-risk-neutral.jl")
        global pcof_NF = Juqbox.run_optimizer(prob, pcof0)
        global results = ep_plot(pcof_NF, params, wa, ep_vals)
        # data_NF = zeros(size(results,1),3)
        data_NF_his[:,i] = results[:,3]
        ctrl_NF_his[:,i] = pcof_NF
        # writedlm("NF_control.dat", pcof_NF)
        # writedlm("NF_optim.dat", data_NF)
    end
    writedlm("NF_optim_his_fixMC.dat", data_NF_his)
    writedlm("NF_ctrl_his_fixMC.dat", ctrl_NF_his)
else
    data_NF_his = readdlm("NF_optim_his_fixMC.dat")
    ctrl_NF_his = readdlm("NF_ctrl_his_fixMC.dat")
    # pcof_NF = vec(readdlm("NF_control.dat"))
    # data_NF = readdlm("NF_optim.dat")
end


freshOptim = false
#### Risk-neutral optimization, data stored in data_RN
nquad = 100
switch_loss = 1

if(freshOptim)
    data_RN_his = zeros(size(ep_vals, 2), trajectory)
    ctrl_RN_his = zeros(ctrl_var, trajectory)
    for i = 1:trajectory
        include("swap-02-risk-neutral.jl")
        global pcof_RN = Juqbox.run_optimizer(prob, pcof0)
        global results2 = ep_plot(pcof_RN, params, wa, ep_vals)
        # data_RN = zeros(size(results2,1),4)
        data_RN_his[:,i] = results2[:,3]
        ctrl_RN_his[:,i] = pcof_RN
        # writedlm("RN_control.dat", pcof_RN)
        # writedlm("RN_optim.dat", data_RN)
    end
    writedlm("RN_optim_his_fixMC.dat", data_RN_his)
    writedlm("RN_ctrl_his_fixMC.dat", ctrl_RN_his)
else
    data_RN_his = readdlm("RN_optim_his_fixMC.dat")
    ctrl_RN_his = readdlm("RN_ctrl_his_fixMC.dat")
    # pcof_RN = vec(readdlm("RN_control.dat"))
    # data_RN = readdlm("RN_optim.dat")
end




freshOptim = false
#### tuning hyper-parameter in Risk-sensitive optimization, data stored in data1 & data2
switch_loss = 2
hyper = 12.0
# hyper_list = [1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 12.0, 15.0, 17.0, 20.0]
# hyper_list = [1.0, 3.0, 5.0, 7.0, 10.0, 12.0, 15.0, 17.0, 20.0]

if(freshOptim)
    data_RS_his = zeros(size(ep_vals, 2), trajectory)
    ctrl_RS_his = zeros(ctrl_var, trajectory)
    for i in 1:trajectory
        # global hyper = hyper_list[k] # change hyper parameter mu in RS
        include("swap-02-risk-neutral.jl")
        global pcof_RS = Juqbox.run_optimizer(prob, pcof0) # solve optimal control problem
        global results3 = ep_plot(pcof_RS, params, wa, ep_vals) # err-perturbation
        # data1[:,k] = results3[:,3] # infidelity err
        # data2[:,k] = results3[:,2] # leak err
        data_RS_his[:,i] = results3[:,3]
        ctrl_RS_his[:,i] = pcof_RS

        # i = Int(hyper)
        # writedlm("RS_control_$i.dat", pcof_RS)
    end
    # writedlm("RS_optim_infidelity.dat", data1)
    # writedlm("RS_optim_leak.dat", data2)
    writedlm("RS_optim_his_fixMC.dat", data_RS_his)
    writedlm("RS_ctrl_his_fixMC.dat", ctrl_RS_his)
else
    data_RS_his = readdlm("RS_optim_his_fixMC.dat")
    ctrl_RS_his = readdlm("RS_ctrl_his_fixMC.dat")
    # data1 = readdlm("RS_optim_infidelity.dat") # use the 7_th column for mu=10
    # data2 = readdlm("RS_optim_leak.dat")
end



# Plot objective functions for various Hamiltonian perturbations
# plc = Plots.plot(scalefactor.*ep_vals,data1[:,1],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=1",
#                  title = "RS sensitivity w.r.t. 𝜇", legend= :outerright)
# Plots.plot!(plc,scalefactor.*ep_vals,data1[:,2],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=3")
# Plots.plot!(plc,scalefactor.*ep_vals,data1[:,3],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=5")
# Plots.plot!(plc,scalefactor.*ep_vals,data1[:,4],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=7")
# Plots.plot!(plc,scalefactor.*ep_vals,data1[:,5],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=10")
# Plots.plot!(plc,scalefactor.*ep_vals,data1[:,6],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=12")
# Plots.plot!(plc,scalefactor.*ep_vals,data1[:,7],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=15")
# Plots.plot!(plc,scalefactor.*ep_vals,data1[:,8],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=17")
# Plots.plot!(plc,scalefactor.*ep_vals,data1[:,9],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RS Objective Function", lab="Infidelity, 𝜇=20")

# Save plots to file
# Plots.savefig(plc,  "RS_sensitivity.png")



#=
for k in hyper_list
    i = Int(k)
    global pcof = vec(readdlm("RS_control_$i.dat"))
    obj, _, _, secondaryobjf, traceinfid = Juqbox.traceobjgrad(pcof,params,wa,false, true)
    println(k, traceinfid)
end
=#


freshOptim = false
#### tuning hyper-parameter in Risk-averse optimization
switch_loss = 3
hyper = 100.0
# hyper_list = [1.0, 5.0, 10.0, 20.0, 40.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]
# hyper_list = [1.0, 5.0, 10.0, 20.0, 40.0, 100.0, 200.0, 400.0, 800.0]

if(freshOptim)
    data_RA_his = zeros(size(ep_vals, 2), trajectory)
    ctrl_RA_his = zeros(ctrl_var, trajectory)
    for i in 1:trajectory
        # global hyper = hyper_list[k] # change hyper parameter mu in RA
        include("swap-02-risk-neutral.jl")
        global pcof_RA = Juqbox.run_optimizer(prob, pcof0) # solve optimal control
        global results4 = ep_plot(pcof_RA, params, wa, ep_vals) # err-perturbation
        # data3[:,k] = results4[:,3] # infidelity err
        # data4[:,k] = results4[:,2] # leak err
        data_RA_his[:,i] = results4[:,3]
        ctrl_RA_his[:,i] = pcof_RA
        
        # i = Int(hyper)
        # writedlm("RA_control_$i.dat", pcof_RA)
    end
    # writedlm("RA_optim_infidelity.dat", data3)
    # writedlm("RA_optim_leak.dat", data4)
    writedlm("RA_optim_his_fixMC.dat", data_RA_his)
    writedlm("RA_ctrl_his_fixMC.dat", ctrl_RA_his)
else
    data_RA_his = readdlm("RA_optim_his_fixMC.dat")
    ctrl_RA_his = readdlm("RA_ctrl_his_fixMC.dat")
    # data3 = readdlm("RA_optim_infidelity.dat") # use the 6_th column for theta=100
    # data4 = readdlm("RA_optim_leak.dat")
end


#### Plot objective functions for various Hamiltonian perturbations
# plc = Plots.plot(scalefactor.*ep_vals,data3[:,1],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=1",
#                  title = "RA sensitivity w.r.t. 𝜃", legend= :outerright)
# Plots.plot!(plc,scalefactor.*ep_vals,data3[:,2],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=5")
# Plots.plot!(plc,scalefactor.*ep_vals,data3[:,3],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=10")
# Plots.plot!(plc,scalefactor.*ep_vals,data3[:,4],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=20")
# Plots.plot!(plc,scalefactor.*ep_vals,data3[:,5],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=40")
# Plots.plot!(plc,scalefactor.*ep_vals,data3[:,6],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=100")
# Plots.plot!(plc,scalefactor.*ep_vals,data3[:,7],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=200")
# Plots.plot!(plc,scalefactor.*ep_vals,data3[:,8],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=400")
# Plots.plot!(plc,scalefactor.*ep_vals,data3[:,9],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "RA Objective Function", lab="Infidelity, 𝜃=800")

# Save plots to file
# Plots.savefig(plc,  "RA_sensitivity.png")



#=
#### comparison between loss functions (NF, RN, RS, RA)
# Plot objective functions for various Hamiltonian perturbations
plc = Plots.plot(scalefactor.*ep_vals,data_NF[:,3],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "Objective Function", lab="NF Infidelity",
                 title = "Robustness against noise", legend= :outerright, color = "aquamarine3")
Plots.plot!(plc,scalefactor.*ep_vals,data_RN[:,3],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "Objective Function", lab="RN Infidelity", color = "red")
Plots.plot!(plc,scalefactor.*ep_vals,data1[:,1],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "Objective Function", lab="RS Infidelity", color = "green")
Plots.plot!(plc,scalefactor.*ep_vals,data3[:,1],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "Objective Function", lab="RA Infidelity", color = "orange")

Plots.plot!(plc,scalefactor.*ep_vals,data_NF[:,2],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "Objective Function", lab="NF Guard Level Leak",linestyle=:dash, color = "aquamarine3")
Plots.plot!(plc,scalefactor.*ep_vals,data_RN[:,2],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "Objective Function", lab="RN Guard Level Leak",linestyle=:dash, color = "red")
Plots.plot!(plc,scalefactor.*ep_vals,data2[:,1],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "Objective Function", lab="RS Guard Level Leak",linestyle=:dash, color = "green")
Plots.plot!(plc,scalefactor.*ep_vals,data4[:,1],yaxis=:log,xlabel = "Hamiltonian Perturbation [MHz]",ylabel = "Objective Function", lab="RA Guard Level Leak",linestyle=:dash, color = "orange")

# Save plots to file
Plots.savefig(plc,  "loss_comparison.png")
=#


#=
#### plot control function, NF
nplot = round(Int64, params.T*32)
td = collect(range(0, stop = params.T, length = nplot+1))
p_NF_1,q_NF_1 = evalctrl_no_carrier(params, pcof_NF, 1, D1) 
p_NF_2,q_NF_2 = evalctrl_no_carrier(params, pcof_NF, 2, D1) 

pfunc1 = scalefactor .* p_NF_1
qfunc1 = scalefactor .* q_NF_1
pmax1 = maximum(abs.(pfunc1))
qmax1 = maximum(abs.(qfunc1))
pfunc2 = scalefactor .* p_NF_2
qfunc2 = scalefactor .* q_NF_2
pmax2 = maximum(abs.(pfunc2))
qmax2 = maximum(abs.(qfunc2))
pmax = maximum([pmax1,pmax2])
qmax = maximum([qmax1,qmax2])

titlestr = "Rotating frame NF ctrl " * " Max-p=" *@sprintf("%.3e", pmax) * " Max-q=" *@sprintf("%.3e", qmax) * " " * unitStr

pl_ctrl_NF = Plots.plot(td, pfunc1, lab=L"p_{1,1}(t)", title = titlestr, xlabel="Time [ns]",
                                  ylabel=unitStr, legend= :outerright, linewidth=1.5, legendfontsize=12)
# add in the control function for the anti-symmetric Hamiltonian
Plots.plot!(pl_ctrl_NF,td, qfunc1, lab=L"q_{1,1}(t)", linewidth=1.5, legendfontsize=12)
Plots.plot!(pl_ctrl_NF,td, pfunc2, lab=L"p_{1,2}(t)", linewidth=1.5, legendfontsize=12)
Plots.plot!(pl_ctrl_NF,td, qfunc2, lab=L"q_{1,2}(t)", linewidth=1.5, legendfontsize=12)


# save plots of control functions in rotating frame without carrier waves
Plots.savefig(pl_ctrl_NF,  "NF_ctrl.png")


#### plot control function, RN
p_RN_1,q_RN_1 = evalctrl_no_carrier(params, pcof_RN, 1, D1) 
p_RN_2,q_RN_2 = evalctrl_no_carrier(params, pcof_RN, 2, D1) 
pfunc1 = scalefactor .* p_RN_1
qfunc1 = scalefactor .* q_RN_1
pmax1 = maximum(abs.(pfunc2))
qmax1 = maximum(abs.(qfunc2))
pfunc2 = scalefactor .* p_RN_2
qfunc2 = scalefactor .* q_RN_2
pmax1 = maximum(abs.(pfunc2))
qmax1 = maximum(abs.(qfunc2))
pmax = maximum([pmax1,pmax2])
qmax = maximum([qmax1,qmax2])

titlestr = "Rotating frame RN ctrl " * " Max-p=" *@sprintf("%.3e", pmax) * " Max-q=" *@sprintf("%.3e", qmax) * " " * unitStr
pl_ctrl_RN = Plots.plot(td, pfunc1, lab=L"p_{1,1}(t)", title = titlestr, xlabel="Time [ns]",
                                  ylabel=unitStr, legend= :topleft, linewidth=1.5, legendfontsize=12)
# add in the control function for the anti-symmetric Hamiltonian
Plots.plot!(pl_ctrl_RN,td, qfunc1, lab=L"q_{1,1}(t)", linewidth=1.5, legendfontsize=12, legend= :topleft)
Plots.plot!(pl_ctrl_RN,td, pfunc2, lab=L"p_{1,2}(t)", linewidth=1.5, legendfontsize=12, legend= :topleft)
Plots.plot!(pl_ctrl_RN,td, qfunc2, lab=L"q_{1,2}(t)", linewidth=1.5, legendfontsize=12, legend= :topleft)
Plots.savefig(pl_ctrl_RN,  "RN_ctrl.png")


#### plot control function, RS
p_RS_1,q_RS_1 = evalctrl_no_carrier(params, pcof_RS, 1, D1) 
p_RS_2,q_RS_2 = evalctrl_no_carrier(params, pcof_RS, 2, D1) 
pfunc1 = scalefactor .* p_RS_1
qfunc1 = scalefactor .* q_RS_1
pmax1 = maximum(abs.(pfunc2))
qmax1 = maximum(abs.(qfunc2))
pfunc2 = scalefactor .* p_RS_2
qfunc2 = scalefactor .* q_RS_2
pmax1 = maximum(abs.(pfunc2))
qmax1 = maximum(abs.(qfunc2))
pmax = maximum([pmax1,pmax2])
qmax = maximum([qmax1,qmax2])

titlestr = "Rotating frame RS ctrl " * " Max-p=" *@sprintf("%.3e", pmax) * " Max-q=" *@sprintf("%.3e", qmax) * " " * unitStr
pl_ctrl_RS = Plots.plot(td, pfunc1, lab=L"p_{1,1}(t)", title = titlestr, xlabel="Time [ns]",
                                  ylabel=unitStr, legend= :outerright, linewidth=1.5, legendfontsize=12)
# add in the control function for the anti-symmetric Hamiltonian
Plots.plot!(pl_ctrl_RS,td, qfunc1, lab=L"q_{1,1}(t)", linewidth=1.5, legendfontsize=12)
Plots.plot!(pl_ctrl_RS,td, pfunc2, lab=L"p_{1,2}(t)", linewidth=1.5, legendfontsize=12)
Plots.plot!(pl_ctrl_RS,td, qfunc2, lab=L"q_{1,2}(t)", linewidth=1.5, legendfontsize=12)
Plots.savefig(pl_ctrl_RS,  "RS_ctrl.png")


#### plot control function, RA
p_RA_1,q_RA_1 = evalctrl_no_carrier(params, pcof_RA, 1, D1) 
p_RA_2,q_RA_2 = evalctrl_no_carrier(params, pcof_RA, 2, D1) 
pfunc1 = scalefactor .* p_RA_1
qfunc1 = scalefactor .* q_RA_1
pmax1 = maximum(abs.(pfunc2))
qmax1 = maximum(abs.(qfunc2))
pfunc2 = scalefactor .* p_RA_2
qfunc2 = scalefactor .* q_RA_2
pmax1 = maximum(abs.(pfunc2))
qmax1 = maximum(abs.(qfunc2))
pmax = maximum([pmax1,pmax2])
qmax = maximum([qmax1,qmax2])

titlestr = "Rotating frame RA ctrl " * " Max-p=" *@sprintf("%.3e", pmax) * " Max-q=" *@sprintf("%.3e", qmax) * " " * unitStr
pl_ctrl_RA = Plots.plot(td, pfunc1, lab=L"p_{1,1}(t)", title = titlestr, xlabel="Time [ns]",
                                  ylabel=unitStr, legend= :outerright, linewidth=1.5, legendfontsize=12)
# add in the control function for the anti-symmetric Hamiltonian
Plots.plot!(pl_ctrl_RA,td, qfunc1, lab=L"q_{1,1}(t)", linewidth=1.5, legendfontsize=12)
Plots.plot!(pl_ctrl_RA,td, pfunc2, lab=L"p_{1,2}(t)", linewidth=1.5, legendfontsize=12)
Plots.plot!(pl_ctrl_RA,td, qfunc2, lab=L"q_{1,2}(t)", linewidth=1.5, legendfontsize=12)
Plots.savefig(pl_ctrl_RA,  "RA_ctrl.png")
=#

#### plot density function for error
# plt_density = Plots.density(log10.(data_NF[:,3]), label="NF-error", 
#                             title = "Error distribution for different methods", xlabel="log10(loss)",
#                             ylabel="density", legend= :topleft, linewidth=1.5, legendfontsize=10)
# Plots.density!(log10.(data_RN[:,3]), label="RN-error")
# Plots.density!(log10.(data1[:,1]), label="RS-error")
# Plots.density!(log10.(data3[:,1]), label="RA-error")

# Plots.savefig(plt_density,  "err_density.png")

#### plot density function(3 trajectories) for error

for i in 1:trajectory
    if i == 1
        println(i)
        global plt_density = Plots.density(data_NF_his[:,i], label="noise-free", boundary=(0,5),
                                    title = "Error distribution for different methods", xlabel="loss",
                                    ylabel="density of loss", legend= :topleft, linewidth=1.5, legendfontsize=10, color = "skyblue", minorgrid=true)
        # Plots.density!(data_RN_his[:,i], xaxis=:log, label = false, color = "tomato")
        # Plots.density!(data_RS_his[:,i], label="RS-error", color = "yellowgreen")
        # Plots.density!(data_RA_his[:,i], label="RA-error", color = "goldenrod1")
        # plot!([mean(log10.(data_NF_his[:,i]))], label = false, seriestype="vline", color = "skyblue")
        # plot!([mean(log10.(data_RN_his[:,i]))], label = false, seriestype="vline", color = "orchid")
        # plot!([mean(log10.(data_RS_his[:,i]))], label = false, seriestype="vline", color = "yellowgreen")
        # plot!([mean(log10.(data_RA_his[:,i]))], label = false, seriestype="vline", color = "goldenrod1")
    elseif i == 2
        # println(i)
        # Plots.density!(data_NF_his[:,i], xaxis=:log, label = false, color = "forestgreen")
        # Plots.density!(data_RN_his[:,i], xaxis=:log, label = false, color = "tomato")
        # Plots.density!(data_RS_his[:,i], label = false, color = "yellowgreen")
        # Plots.density!(data_RA_his[:,i], label = false, color = "goldenrod1")
        # plot!([mean(log10.(data_NF_his[:,i]))], label = false, seriestype="vline", color = "skyblue")
        # plot!([mean(log10.(data_RN_his[:,i]))], label = false, seriestype="vline", color = "orchid")
        # plot!([mean(log10.(data_RS_his[:,i]))], label = false, seriestype="vline", color = "yellowgreen")
        # plot!([mean(log10.(data_RA_his[:,i]))], label = false, seriestype="vline", color = "goldenrod1")
    else
        # println(i)
        # Plots.density!(data_NF_his[:,i], xaxis=:log, label = false, color = "forestgreen")
        # Plots.density!(data_RN_his[:,i], xaxis=:log, label = false, color = "tomato")
        # Plots.density!(data_RS_his[:,i], label = false, color = "yellowgreen")
        # Plots.density!(data_RA_his[:,i], label = false, color = "goldenrod1")
        # plot!([mean(log10.(data_NF_his[:,i]))], label = false, seriestype="vline", color = "skyblue")
        # plot!([mean(log10.(data_RN_his[:,i]))], label = false, seriestype="vline", color = "orchid")
        # plot!([mean(log10.(data_RS_his[:,i]))], label = false, seriestype="vline", color = "yellowgreen")
        # plot!([mean(log10.(data_RA_his[:,i]))], label = false, seriestype="vline", color = "goldenrod1")
    end 
end

Plots.savefig(plt_density,  "err_density_NFNA.png")
