module QPModelGurobi

using Gurobi
using QuadraticModels
using NLPModels
using SolverTools
using LinearAlgebra
using SparseArrays

function createQuadraticModel(qpdata)
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
                        Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
                        lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
                        c0=qpdata.c0)

end



function optimizeGurobi(QM; method=2, presolve=true, scaling=true, crossover=true)

    SM = SlackModel(QM)
    env = Gurobi.Env()

    # -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier,
    # 3=concurrent, 4=deterministic concurrent, 5=deterministic concurrent simplex.
    # default to barrier
    setparam!(env, "Method", 2)
    setparam!(env, "Presolve", 0) # set presolve to 0
    setparam!(env, "ScaleFlag", 0) # no scaling
    setparam!(env, "Crossover", 0)
    setparam!(env, "BarIterLimit", 0)

    Aeq = jac(SM, SM.meta.x0)
    beq = SM.meta.lcon
    f = grad(SM, zeros(length(SM.meta.x0)))
    H = hess(SM, zeros(length(SM.meta.x0)))
    H = Matrix(Symmetric(H, :L))
    n,m = size(Aeq)
    model = gurobi_model(env; f = f, H = H,
                        Aeq = Aeq, beq = beq,
                        lb = SM.meta.lvar, ub = SM.meta.uvar)
     # run optimization
    optimize(model)

    # y with dual: b'*y   s.t. A'*y <= c and y >= 0
    y = zeros(n)
    for i=1:n
        y[i] = Gurobi.get_dblattrelement(model, "Pi", i)
    end

    s = zeros(m) # s_l - s_u
    for i=1:m
        s[i] = Gurobi.get_dblattrelement(model, "RC", i)
    end

    # outputs
    optim_info = get_optiminfo(model)
    if optim_info.status == :optimal
        status = :acceptable
    elseif optim_info.status == :iteration_limit
        status = :max_iter
    elseif optim_info.status == :unbounded
        status = :unbounded
    else
        status = :unknown
    end

    x = get_solution(model)
    stats = GenericExecutionStats(status, SM, solution = x[1:SM.meta.nvar],
                                  objective = get_objval(model),
                                  iter = Gurobi.get_intattr(model,"BarIterCount"),
                                  primal_feas = norm(Aeq * x - beq),
                                  dual_feas = norm(Aeq' * y - H*x + s - f),
                                  multipliers = y,
                                  elapsed_time = optim_info.runtime)
    return stats
end


function optimizeGurobi(qpdata::QPSData; method=2, presolve=true, scaling=true,
                        crossover=true)
    return optimizeGurobi(createQuadraticModel(qpdata); method=2, presolve=true,
                          scaling=true, crossover=true)
end


end
