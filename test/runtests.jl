using QPSReader, QuadraticModelsGurobi
using Test

@testset "QuadraticModelsGurobi.jl" begin
    qps1 = readqps("QAFIRO.SIF") #lower bounds
    stats1 = gurobi(qps1)
    @test abs(stats1.objective + 1.59078179) < 1e-2

    qps2 = readqps("HS21.SIF") # low/upp bounds
    stats2 = gurobi(qps2)
    @test abs(stats2.objective + 9.99599999e1) < 1e-2

    qps3 = readqps("HS52.SIF") # free bounds
    stats3 = gurobi(qps3)
    @test abs(stats3.objective - 5.32664756) < 1e-2
end
