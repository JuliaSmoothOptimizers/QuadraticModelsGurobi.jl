module QuadraticModelsGurobi

export gurobi

using Gurobi
using SolverTools
using SparseArrays

const gurobi_statuses = Dict(:loaded => :unknown,
                             :optimal => :acceptable,
                             :infeasible => :infeasible,
                             :inf_or_unbd => :infeasible,
                             :unbounded => :unbounded,
                             :cutoff => :exception,
                             :iteration_limit => :max_iter,
                             :node_limit => :exception,
                             :time_limit => :max_time,
                             :solution_limit => :exception,
                             :interrupted => :user,
                             :numeric => :exception,
                             :suboptimal => :exception,
                             :inprogress => :exception,
                             :user_obj_limit => :exception)

function sparse_csr(I, J, V, m=maximum(I), n=maximum(J))
    csrrowptr = zeros(Int, m+1)
    # Compute the CSR form's row counts and store them shifted forward by one in csrrowptr
    coolen = length(I)
    min(length(J), length(V)) >= coolen || throw(ArgumentError("J and V need length >= length(I) = $coolen"))
    @inbounds for k in 1:coolen
        Ik = I[k]
        if 1 > Ik || m < Ik
            throw(ArgumentError("row indices I[k] must satisfy 1 <= I[k] <= m"))
        end
        csrrowptr[Ik+1] += 1
    end

    # Compute the CSR form's rowptrs and store them shifted forward by one in csrrowptr
    countsum = 1
    csrrowptr[1] = 1
    @inbounds for i in 2:(m+1)
        overwritten = csrrowptr[i]
        csrrowptr[i] = countsum
        countsum += overwritten
    end

    # Counting-sort the column and nonzero values from J and V into csrcolval and csrnzval
    # Tracking write positions in csrrowptr corrects the row pointers
    csrcolval = zeros(Int, length(I))
    csrnzval = zeros(length(I))
    @inbounds for k in 1:coolen
        Ik, Jk = I[k], J[k]
        if 1 > Jk || n < Jk
            throw(ArgumentError("column indices J[k] must satisfy 1 <= J[k] <= n"))
        end
        csrk = csrrowptr[Ik+1]
        csrrowptr[Ik+1] = csrk + 1
        csrcolval[csrk] = Jk
        csrnzval[csrk] = V[k]
    end
    csrrowptr = csrrowptr[1:end-1]
    return csrrowptr, csrcolval, csrnzval
end

function gurobi(QM; method=2, kwargs...)
    env = Gurobi.Env()
    # -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier,
    # 3=concurrent, 4=deterministic concurrent, 5=deterministic concurrent simplex.
    # default to barrier
    setparam!(env, "Method", method)
    # use kwargs change to presolve, scaling and crossover mode
    # example: gurobi(QM, presolve=0) (see gurobi doc for other options)
    for (k, v) in kwargs
        if k==:presolve
            setparam!(env, "Presolve", v) # 0 = no presolve
        elseif k==:scaling
            setparam!(env, "ScaleFlag", v) # 0 = no scaling
        elseif k==:crossover
            setparam!(env, "Crossover", v) # 0 = no crossover
        end
    end

    model = Gurobi.Model(env, "")
    add_cvars!(model, QM.data.c, QM.meta.lvar, QM.meta.uvar)
    Gurobi.set_dblattr!(model, "ObjCon", QM.data.c0)
    update_model!(model)

    if QM.meta.nnzh > 0
        Hvals = zeros(eltype(QM.data.Hvals), length(QM.data.Hvals))
        for i=1:length(QM.data.Hvals)
            if QM.data.Hrows[i] == QM.data.Hcols[i]
                Hvals[i] = QM.data.Hvals[i] / 2
            else
                Hvals[i] = QM.data.Hvals[i]
            end
        end
        add_qpterms!(model, QM.data.Hrows, QM.data.Hcols, Hvals)
    end

    Acsrrowptr, Acsrcolval, Acsrnzval = sparse_csr(QM.data.Arows,QM.data.Acols,
                                                   QM.data.Avals, QM.meta.ncon,
                                                   QM.meta.nvar)

    # add_rangeconstrs!(model, sparse(QM.data.Arows, QM.data.Acols,
    # 								  QM.data.Avals, QM.meta.ncon,
    # 								  QM.meta.nvar),
    # 					QM.meta.lcon, QM.meta.ucon)
    add_rangeconstrs!(model, Acsrrowptr, Acsrcolval, Acsrnzval, QM.meta.lcon,
                      QM.meta.ucon)


    update_model!(model)

    optimize(model)

    n_constr = Gurobi.get_intattr(model, "NumConstrs")
    y = zeros(n_constr)
    for i=1:(n_constr)
        y[i] = Gurobi.get_dblattrelement(model, "Pi", i)
    end
    s = zeros(length(QM.data.c)) # s_l - s_u
    for i=1:length(QM.data.c)
        s[i] = Gurobi.get_dblattrelement(model, "RC", i)
    end

    optim_info = get_optiminfo(model)
    x = get_solution(model)
    stats = GenericExecutionStats(get(gurobi_statuses, optim_info.status, :unknown),
                                  QM, solution = x,
                                  objective = get_objval(model),
                                  iter = Gurobi.get_intattr(model,"BarIterCount"),
                                  primal_feas = Gurobi.get_dblattr(model, "ConstrResidual"),
                                  dual_feas = Gurobi.get_dblattr(model, "DualResidual"),
                                  solver_specific = Dict(:multipliers => y,
                                  :reduced_costs => s),
                                  elapsed_time = optim_info.runtime)
                                  return stats
end

end
