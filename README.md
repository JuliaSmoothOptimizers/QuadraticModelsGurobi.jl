# QuadraticModelsGurobi

A package to use Gurobi to optimize linear and quadratic problems in QPSData
format (see QPSReader.jl)

# Usage

```julia
using QPSReader, QuadraticModelsGurobi
qps = readqps("AFIRO.SIF")
stats = gurobi(qps)
```
