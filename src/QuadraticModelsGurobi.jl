module QuadraticModelsGurobi

export gurobi

using Gurobi
using QuadraticModels, SolverCore, SparseMatricesCOO
using LinearAlgebra, SparseArrays

const gurobi_statuses = Dict(1 => :unknown,
                             2 => :acceptable,
                             3 => :infeasible,
                             4 => :infeasible,
                             5 => :unbounded,
                             6 => :exception,
                             7 => :max_iter,
                             8 => :exception,
                             9 => :max_time,
                             10 => :exception,
                             11 => :user,
                             12 => :exception,
                             13 => :exception,
                             14 => :exception,
                             15 => :exception)

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

gurobi(QM::QuadraticModel{T, S}; kwargs...) where {T, S} = gurobi(
    convert(QuadraticModel{T, S, SparseMatrixCOO{T, Int}, SparseMatrixCOO{T, Int}}, QM);
    kwargs...
)

function gurobi(QM::QuadraticModel{T, S, M1, M2};
                method=2, kwargs...) where {T, S, M1 <: SparseMatrixCOO, M2 <: SparseMatrixCOO}
    env = Gurobi.Env()
    # -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier,
    # 3=concurrent, 4=deterministic concurrent, 5=deterministic concurrent simplex.
    # default to barrier
    GRBsetintparam(env, "Method", method)
    # use kwargs change to presolve, scaling and crossover mode
    # example: gurobi(QM, presolve=0) (see gurobi doc for other options)
    for (k, v) in kwargs
        if k==:presolve
            GRBsetintparam(env, "Presolve", v) # 0 = no presolve
        elseif k==:scaling
            GRBsetintparam(env, "ScaleFlag", v) # 0 = no scaling
        elseif k==:crossover
            GRBsetintparam(env, "Crossover", v) # 0 = no crossover
        elseif k==:display
            GRBsetintparam(env, "OutputFlag", v) # 0 = no display
        end
    end

    model = Ref{Ptr{Cvoid}}()
    GRBnewmodel(env, model, "", QM.meta.nvar, QM.data.c, QM.meta.lvar, QM.meta.uvar, C_NULL, C_NULL)
    GRBsetdblattr(model.x, "ObjCon", QM.data.c0)
    if QM.meta.nnzh > 0
        Hvals = zeros(eltype(QM.data.H.vals), length(QM.data.H.vals))
        for i=1:length(QM.data.H.vals)
            if QM.data.H.rows[i] == QM.data.H.cols[i]
                Hvals[i] = QM.data.H.vals[i] / 2
            else
                Hvals[i] = QM.data.H.vals[i]
            end
        end
        GRBaddqpterms(model.x, length(QM.data.H.cols), convert(Array{Cint,1}, QM.data.H.rows.-1),
                      convert(Array{Cint,1}, QM.data.H.cols.-1), Hvals)
    end

    Acsrrowptr, Acsrcolval, Acsrnzval = sparse_csr(QM.data.A.rows,QM.data.A.cols,
                                                   QM.data.A.vals, QM.meta.ncon,
                                                   QM.meta.nvar)
    GRBaddrangeconstrs(model.x, QM.meta.ncon, length(Acsrcolval), convert(Array{Cint,1}, Acsrrowptr.-1),
                       convert(Array{Cint,1}, Acsrcolval.-1), Acsrnzval, QM.meta.lcon, QM.meta.ucon, C_NULL)

    GRBoptimize(model.x)

    x = zeros(QM.meta.nvar)
    GRBgetdblattrarray(model.x, "X", 0, QM.meta.nvar, x)
    y = zeros(QM.meta.ncon)
    GRBgetdblattrarray(model.x, "Pi", 0, QM.meta.ncon, y)
    s = zeros(QM.meta.nvar)
    GRBgetdblattrarray(model.x, "RC", 0, QM.meta.nvar, s)
    status = Ref{Cint}()
    GRBgetintattr(model.x, "Status", status)
    baritcnt = Ref{Cint}()
    GRBgetintattr(model.x, "BarIterCount", baritcnt)
    objval = Ref{Float64}()
    GRBgetdblattr(model.x, "ObjVal", objval)
    p_feas = Ref{Float64}()
    GRBgetdblattr(model.x, "ConstrResidual", p_feas)
    d_feas = Ref{Float64}()
    GRBgetdblattr(model.x, "DualResidual", d_feas)
    elapsed_time = Ref{Float64}()
    GRBgetdblattr(model.x, "Runtime", elapsed_time)
    stats = GenericExecutionStats(get(gurobi_statuses, status[], :unknown),
                                  QM, solution = x,
                                  objective = objval[],
                                  iter = Int64(baritcnt[]),
                                  primal_feas = p_feas[],
                                  dual_feas = d_feas[],
                                  multipliers = y,
                                  elapsed_time = elapsed_time[])
    return stats
end

end
