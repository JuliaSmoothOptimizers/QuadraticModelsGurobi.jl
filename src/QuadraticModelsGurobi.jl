module QuadraticModelsGurobi

export gurobi

using Gurobi
using QuadraticModels
using QPSReader
using SolverTools
using LinearAlgebra
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


function createQuadraticModel(qpdata)
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
                        Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
                        lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
                        c0=qpdata.c0)
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

	# method 1
    # first_irow = 1
    # last_irow = 0
    # p = sortperm(QM.data.Arows)
    # for i=1:length(QM.meta.lcon)
    #     if last_irow < length(QM.data.Arows) && @views QM.data.Arows[p][last_irow+1] == i
    #         first_irow = last_irow + 1
    #         last_irow = @views first_irow-1+findlast(QM.data.Arows[p][first_irow:end] .== i)
    #         if QM.meta.lcon[i] == QM.meta.ucon[i]
	# 			add_constr!(model, QM.data.Acols[p][first_irow:last_irow],
	# 			 			QM.data.Avals[p][first_irow:last_irow],
	# 						'=', QM.meta.lcon[i])
    #         elseif QM.meta.lcon[i] == -Inf && QM.meta.ucon[i] != Inf
	# 			add_constr!(model, QM.data.Acols[p][first_irow:last_irow],
	# 					   QM.data.Avals[p][first_irow:last_irow],
	# 					   '<', QM.meta.ucon[i])
    #         elseif QM.meta.lcon[i] != -Inf && QM.meta.ucon[i] == Inf
	# 			add_constr!(model, QM.data.Acols[p][first_irow:last_irow],
	# 						.-QM.data.Avals[p][first_irow:last_irow],
	# 						'<', -QM.meta.lcon[i])
    #         elseif QM.meta.lcon[i] != -Inf && QM.meta.ucon[i] != Inf
	# 			add_constr!(model, QM.data.Acols[p][first_irow:last_irow],
	# 						QM.data.Avals[p][first_irow:last_irow],
	# 						'<', QM.meta.ucon[i])
	# 			add_constr!(model, QM.data.Acols[p][first_irow:last_irow],
	# 						.-QM.data.Avals[p][first_irow:last_irow],
	# 						'<', -QM.meta.lcon[i])
    #         end
    #     end
    # end

	eq = [j for j=1:length(QM.data.Arows) if QM.data.Arows[j] in QM.meta.jfix]
	add_constrs!(model, QM.data.Arows[eq], QM.data.Acols[eq], QM.data.Avals[eq],
				'=', QM.meta.lcon[QM.data.Arows[eq]])
	low = [j for j=1:length(QM.data.Arows) if QM.data.Arows[j] in QM.meta.jlow]
	add_constrs!(model, QM.data.Arows[low], QM.data.Acols[low], .-QM.data.Avals[low],
				'<', .-QM.meta.lcon[QM.data.Arows[low]])
	upp = [j for j=1:length(QM.data.Arows) if QM.data.Arows[j] in QM.meta.jupp]
	add_constrs!(model, QM.data.Arows[upp], QM.data.Acols[upp], QM.data.Avals[upp],
				'<', QM.meta.ucon[QM.data.Arows[upp]])
	rng = [j for j=1:length(QM.data.Arows) if QM.data.Arows[j] in QM.meta.jrng]
	add_constrs!(model, QM.data.Arows[rng], QM.data.Acols[rng], .-QM.data.Avals[rng],
				'<', .-QM.meta.lcon[QM.data.Arows[rng]])
	add_constrs!(model, QM.data.Arows[rng], QM.data.Acols[rng], QM.data.Avals[rng],
				'<', QM.meta.ucon[QM.data.Arows[rng]])

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


function gurobi(qpdata::QPSData; method=2, kwargs...)
    return gurobi(createQuadraticModel(qpdata); method=2, kwargs...)
end

end
