# QPModelGurobi

A package to use Gurobi to optimize linear and quadratic problems in QPSData
format (see QPSReader.jl)

# Usage

```julia
using QPSReader, QPModelGurobi
qps = readqps("AFIRO.SIF")
stats = gurobi(qps)
```
