# QPModelGurobi

A package to use Gurobi to optimize linear and quadratic problems in QPSData
format (see QPSReader.jl)

# Usage

```julia
using QPSReader
using QPModelGurobi
qpmodel = readqps("AFIRO.SIF")
stats = QPModelGurobi.optimizeGurobi(qpmodel)
```
