# setup callback functions for Ipopt
function eval_f_par(pcof::Vector{Float64}, params:: Juqbox.objparams, wa::Working_Arrays,
                    nodes::AbstractArray=[0.0],weights::AbstractArray=[1.0],
                    switch_loss::Int64=1, hyper::Float64=0.0)
    
    # Loop over specified nodes and compute risk-aware objective value. Default is usual optimization.
    # the noise is in the diagonal, 4-dimensional, and uniformly distributed.

    nquad = size(nodes, 2) # (dim of noise)*nquad
    # exp_v = zeros(Float64, nquad)
    exp_v = wa.exp_v
    exp_v .= 0.0

    # H0_old = copy(params.Hconst)
    for i = 1:nquad 
        ep = nodes[:,i] # noise vector

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

        E = Juqbox.traceobjgrad(pcof,params,wa,false,false)
        exp_v[i] = E[1]

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
    # copy!(params.Hconst,H0_old)
    loss = reweightf(switch_loss, exp_v, weights, hyper)
    return loss
  end


function reweightf(switch_loss::Int64,exp_v::Vector{Float64},weights::Vector{Float64}, hyper::Float64=0.0)
    if switch_loss == 1
        # do RN 
        return reweightf_RN(exp_v, weights)
    elseif switch_loss == 2
        # do RA
        return reweightf_RS(exp_v, weights, hyper)
    elseif switch_loss == 3
        # do RS
        return reweightf_RA(exp_v, weights, hyper)
    else 
        println("PICK SOMETHING ELSE!")
    end
end

function eval_g(pcof)
    return 0.0
end


function eval_grad_f_par(pcof::Vector{Float64}, grad_f::Vector{Float64}, params:: Juqbox.objparams, wa::Working_Arrays,
                        nodes::AbstractArray=[0.0],weights::AbstractArray=[1.0],
                        switch_loss::Int64=1, hyper::Float64=0.0)
    grad_f .= 0.0
    nquad = size(nodes, 2)

    grad_f_his = wa.grad_f_his
    objfv_his = wa.objfv_his
    grad_f_his .= 0.0
    objfv_his .= 0.0

    # pre-allocate memory for sample gradient
    # grad_f_his = zeros(Float64, (length(grad_f), nquad))
    # objfv_his = zeros(wa., nquad)

    # H0_old = copy(params.Hconst)
    exp_inf = 0.0
    exp_sec = 0.0
    for i = 1:nquad 
        ep = nodes[:,i]

        # Additive noise
        for j = 1:size(params.Hconst,2)
            params.Hconst[j,j] += 0.0001*ep[j]*(10.0^(j))
        end
        # params.Hconst[1,2] += 0.001*ep[5]
        # params.Hconst[2,1] += 0.001*ep[5]
        # params.Hconst[2,3] += 0.001*ep[6]
        # params.Hconst[3,2] += 0.001*ep[6]
        # params.Hconst[3,4] += 0.001*ep[7]
        # params.Hconst[4,3] += 0.001*ep[7]

        objfv, Gtemp, _, secondaryobjf, traceinfid = Juqbox.traceobjgrad(pcof,params,wa,false, true)

        Gtemp = vcat(Gtemp...) 
        for j in 1:length(Gtemp)
            grad_f_his[j,i] += Gtemp[j]
            #grad_f[j] += Gtemp[j]*weights[i]
        end
        objfv_his[i] += objfv

        # Accumulate expected value of infidelity and guard level penalty
        exp_inf += weights[i]*traceinfid
        exp_sec += weights[i]*secondaryobjf

        # Reset
        for j = 1:size(params.Hconst,2)
            params.Hconst[j,j] -= 0.0001*ep[j]*(10.0^(j))
        end
        # params.Hconst[1,2] -= 0.001*ep[5]
        # params.Hconst[2,1] -= 0.001*ep[5]
        # params.Hconst[2,3] -= 0.001*ep[6]
        # params.Hconst[3,2] -= 0.001*ep[6]
        # params.Hconst[3,4] -= 0.001*ep[7]
        # params.Hconst[4,3] -= 0.001*ep[7]
    end
    # copy!(params.Hconst,H0_old)

    # compute weight for gradients
    weights2 = reweightg(switch_loss, objfv_his, weights, hyper)

    # update reweighted gradients

    grad_f2 = grad_f_his * weights2
    grad_f .= grad_f2

    # remember the value of the primary obj func (for termination in intermediate_par)
    params.lastTraceInfidelity = exp_inf
    params.lastLeakIntegral = exp_sec

    # Save intermediate parameter vectors
    if params.save_pcof_hist
        push!(params.pcof_hist, copy(pcof)) #pcof_hist is and Array of Vector{Float64}
    end
end

function reweightg(switch_loss::Int64,objfv_his::Vector{Float64},weights::Vector{Float64}, hyper::Float64=0.0)
    if switch_loss == 1
        # do RN 
        return reweightg_RN(objfv_his, weights)
    elseif switch_loss == 2
        # do RA
        return reweightg_RS(objfv_his, weights, hyper)
    elseif switch_loss == 3
        # do RS
        return reweightg_RA(objfv_his, weights, hyper)
    else 
        println("PICK SOMETHING ELSE!")
    end
end

# delete _R* for the two reweight functions you want to use
# reweightg function returns adjusted weight for each node's gradient, a vector
# reweightf function returns reweighted loss, which is a real number
# reweightg_RN
function reweightg_RN(loss::Vector{Float64}, weights::AbstractArray=[1.0])
    return weights
end

# reweightf_RN
function reweightf_RN(loss::Vector{Float64}, weights::AbstractArray=[1.0])
    return sum(loss .* weights)
end

# reweightg_RS, without normalization
function reweightg_RS(loss::Vector{Float64}, weights::AbstractArray=[1.0], hyper::Float64=0.0)
    # hyper parameter mu>0
    mu = hyper

    # allocate space
    weight2 = zeros(Float64, length(weights))

    # compute adjusted weight
    weight2 .= mu.*weights .* exp.(mu.*loss)

    return weight2
end

function reweightf_RS(loss::Vector{Float64}, weights::AbstractArray=[1.0], hyper::Float64=0.0)
    mu = hyper
    return sum(exp.(mu.*loss) .* weights) - 1
end

# reweightg_RA
function reweightg_RA(loss::Vector{Float64}, weights::AbstractArray=[1.0], hyper::Float64=0.0)
    # hyper parameter theta
    theta = hyper

    # allocate space
    weights2 = zeros(Float64, length(weights))

    # mean of loss
    mloss = sum(weights .* loss)

    # compute adjusted weight
    weights2 .= weights .* (1 .+ theta .*(loss .- mloss))
    return weights2
end

function reweightf_RA(loss::Vector{Float64}, weights::AbstractArray=[1.0], hyper::Float64=0.0)
    # hyper parameter theta
    theta = hyper / 2

    mloss = sum(loss.*weights)
    return sum(weights .* (loss .+ theta .* (loss .- mloss).^2))
end


function eval_jac_g(
    x::Vector{Float64},         # Current solution
    mode,                       # Either :Structure or :Values
    rows::Vector{Int32},        # Sparsity structure - row indices
    cols::Vector{Int32},        # Sparsity structure - column indices
    values::Vector{Float64})    # The values of the Hessian

    if mode == :Structure
        # rows[...] = ...
        # ...
        # cols[...] = ...
    else
        # values[...] = ...
    end
end

function eval_jac_g(
    x::Vector{Float64},         # Current solution
    rows::Vector{Int32},        # Sparsity structure - row indices
    cols::Vector{Int32},        # Sparsity structure - column indices
    values::Union{Nothing,Vector{Float64}})    # The values of the Hessian
    return nothing
end

function eval_h(
    x::Vector{Float64},
    rows::Vector{Int32},
    cols::Vector{Int32},
    obj_factor::Float64,
    lambda::Vector{Float64},
    values::Union{Nothing,Vector{Float64}})

    if mode == :Structure
        # rows[...] = ...
        # ...
        # cols[...] = ...
    else
        # values[...] = ...
    end
end

function intermediate_par(
    alg_mod::Union{Int32,Int64},
    iter_count::Union{Int32,Int64},
    obj_value::Float64,
    inf_pr::Float64, inf_du::Float64,
    mu::Float64, d_norm::Float64,
    regularization_size::Float64,
    alpha_du::Float64, alpha_pr::Float64,
    ls_trials::Union{Int32,Int64},
    params:: Juqbox.objparams)
  # ...
    if params.saveConvHist 
        push!(params.objHist, obj_value)
        push!(params.dualInfidelityHist, inf_du)
        push!(params.primaryHist, params.lastTraceInfidelity)
        push!(params.secondaryHist,  params.lastLeakIntegral)
    end
    if params.lastTraceInfidelity < params.traceInfidelityThreshold
        println("Stopping because trace infidelity = ", params.lastTraceInfidelity,
                " < threshold = ", params.traceInfidelityThreshold)
        return false
    else
        return true  # Keep going
    end
end

"""
    prob = setup_ipopt_problem(params, wa, nCoeff, minCoeff, maxCoeff; maxIter=50, 
                            lbfgsMax=10, startFromScratch=true, ipTol=1e-5, acceptTol=1e-5, acceptIter=15,
                            nodes=[0.0], weights=[1.0])

Setup structure containing callback functions and convergence criteria for 
optimization via IPOPT. Note the last two inputs, `nodes', and 
`weights', are to be used when performing a simple risk-neutral optimization
where the fundamental frequency is random.

# Arguments
- `params:: objparams`: Struct with problem definition
- `wa::Working_Arrays`: Struct containing preallocated working arrays
- `nCoeff:: Int64`: Number of parameters in optimization
- `minCoeff:: Array{Float64, 1}`: Minimum allowable value for each parameter
- `maxCoeff:: Array{Float64, 1}`: Maximum allowable value for each parameter
- `maxIter:: Int64`: Maximum number of iterations to be taken by optimizer (keyword arg)
- `lbfgsMax:: Int64`: Maximum number of past iterates for Hessian approximation by L-BFGS (keyword arg)
- `startFromScratch:: Bool`: Specify whether the optimization is starting from file or not (keyword arg)
- `ipTol:: Float64`: Desired convergence tolerance (relative) (keyword arg)
- `acceptTol:: Float64`: Acceptable convergence tolerance (relative) (keyword arg)
- `acceptIter:: Int64`: Number of acceptable iterates before triggering termination (keyword arg)
- `nodes:: Array{Float64, 1}`: Risk-neutral opt: User specified quadrature nodes on the interval [-ϵ,ϵ] for some ϵ (keyword arg)
- `weights:: Array{Float64, 1}`: Risk-neutral opt: User specified quadrature weights on the interval [-ϵ,ϵ] for some ϵ (keyword arg)
"""
function setup_ipopt_problem(params:: Juqbox.objparams, wa::Working_Arrays, nCoeff:: Int64, minCoeff:: Array{Float64, 1}, maxCoeff:: Array{Float64, 1};
                             maxIter:: Int64=50, lbfgsMax:: Int64=10, 
                             startFromScratch:: Bool=true, ipTol:: Float64=1e-5, 
                             acceptTol:: Float64=1e-5, acceptIter:: Int64=15,
                             nodes::AbstractArray=[0.0], 
                             weights::AbstractArray=[1.0],
                             switch_loss::Int64=1,
                             hyper::Float64=0.0)
    # callback functions need access to the params object
    eval_f(pcof) = eval_f_par(pcof, params, wa, nodes, weights, switch_loss, hyper)
    eval_grad_f(pcof, grad_f) = eval_grad_f_par(pcof, grad_f, params, wa, nodes, weights, switch_loss, hyper)
    intermediate(alg_mod, iter_count, obj_value, inf_pr, inf_du, mu,
                 d_norm, regularization_size, alpha_du, alpha_pr, ls_trials) =
                     intermediate_par(alg_mod, iter_count, obj_value, inf_pr, inf_du, mu,
                                      d_norm, regularization_size, alpha_du, alpha_pr, ls_trials, params)
    # setup the Ipopt data structure
    mNLconstraints = 0;
    nEleJac = 0;
    nEleHess = 0;
    dum0 = zeros(0);
    dum1 = zeros(0);
    # tmp
    #    println("setup_ipopt_problem: nCoeff = ", nCoeff, " length(minCoeff) = ", length(minCoeff))
    if @isdefined createProblem
        prob = createProblem( nCoeff, minCoeff, maxCoeff, mNLconstraints, dum0, dum1, nEleJac, nEleHess, eval_f, eval_g, eval_grad_f, eval_jac_g);
    else
        prob = CreateIpoptProblem( nCoeff, minCoeff, maxCoeff, mNLconstraints, dum0, dum1, nEleJac, nEleHess, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h);
    end

    if @isdefined addOption
        addOption( prob, "hessian_approximation", "limited-memory");
        addOption( prob, "limited_memory_max_history", lbfgsMax);
        addOption( prob, "max_iter", maxIter);
        addOption( prob, "tol", ipTol);
        addOption( prob, "acceptable_tol", acceptTol);
        addOption( prob, "acceptable_iter", acceptIter);
        if !startFromScratch # enable warm start of Ipopt
            addOption( prob, "warm_start_init_point", "yes")
            # addOption( prob, "mu_init", 1e-6) # not sure how to set this parameter
            # addOption( prob, "nlp_scaling_method", "none") # what about scaling?
            #
            # the following settings prevent the initial parameters to be pushed away from their bounds
            addOption( prob, "warm_start_bound_push", 1e-16)
            addOption( prob, "warm_start_bound_frac", 1e-16)
            addOption( prob, "warm_start_slack_bound_frac", 1e-16)
            addOption( prob, "warm_start_slack_bound_push", 1e-16)

            if !params.quiet
                println("Ipopt: Enabling warm start option")
            end
        end
    else
        AddIpoptStrOption( prob, "hessian_approximation", "limited-memory");
        AddIpoptIntOption( prob, "limited_memory_max_history", lbfgsMax);
        AddIpoptIntOption( prob, "max_iter", maxIter);
        AddIpoptNumOption( prob, "tol", ipTol);
        AddIpoptNumOption( prob, "acceptable_tol", acceptTol);
        AddIpoptIntOption( prob, "acceptable_iter", acceptIter);
        if !startFromScratch # enable warm start of Ipopt
            AddIpoptStrOption( prob, "warm_start_init_point", "yes")
            # AddIpoptNumOption( prob, "mu_init", 1e-6) # not sure how to set this parameter
            # AddIpoptStrOption( prob, "nlp_scaling_method", "none") # what about scaling?
            #
            # the following settings prevent the initial parameters to be pushed away from their bounds
            AddIpoptNumOption( prob, "warm_start_bound_push", 1e-16)
            AddIpoptNumOption( prob, "warm_start_bound_frac", 1e-16)
            AddIpoptNumOption( prob, "warm_start_slack_bound_frac", 1e-16)
            AddIpoptNumOption( prob, "warm_start_slack_bound_push", 1e-16)

            if !params.quiet
                println("Ipopt: Enabling warm start option")
            end
        end        
    end

    # intermediate callback function
    if @isdefined setIntermediateCallback
        setIntermediateCallback(prob, intermediate)
    else 
        SetIntermediateCallback(prob,intermediate)
    end

# output some of the settings
    if !params.quiet
        println("Ipopt parameters: max # iterations = ", maxIter)
        println("Ipopt parameters: max history L-BFGS = ", lbfgsMax)
        println("Ipopt parameters: tol = ", ipTol)
        println("Ipopt parameters: atol = ", acceptTol)
        println("Ipopt parameters: accept # iter = ", acceptIter)
    end
    
    return prob
end

"""
    pcof = run_optimizer(prob, pcof0 [, baseName:: String=""])

Call IPOPT to  optimizize the control functions.

# Arguments
- `prob:: IpoptProblem`: Struct containing Ipopt parameters callback functions
- `pcof0:: Array{Float64, 1}`: Initial guess for the parameter values
- `baseName:: String`: Name of file for saving the optimized parameters; extension ".jld2" is appended
"""
function run_optimizer(prob:: IpoptProblem, pcof0:: Array{Float64, 1}, baseName:: String="")
    # takes at most max_iter >= 0 iterations. Set with addOption(prob, "max_iter", nIter)

    # initial guess for IPOPT
    prob.x = pcof0;

    # Ipopt solver
    println("*** Starting the optimization ***")
    if @isdefined solveProblem
        @time solveProblem(prob);
    else 
        @time IpoptSolve(prob);
    end
    pcof = prob.x;

    #save the b-spline coeffs on a JLD2 file
    if length(baseName)>0
        fileName = baseName * ".jld2"
        save_pcof(fileName, pcof)
        println("Saved B-spline parameters on binary jld2-file '", fileName, "'");
    end

    return pcof

end # run_optimizer
