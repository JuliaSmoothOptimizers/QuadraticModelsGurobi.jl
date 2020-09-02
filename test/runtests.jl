using QPSReader, QuadraticModelsGurobi
using Test

@testset "QuadraticModelsGurobi.jl" begin
    qps1 = readqps("QAFIRO.SIF")
    stats1 = gurobi(qps1)
    @test abs(stats1.objective + 1.59078) < 1e-2
end
