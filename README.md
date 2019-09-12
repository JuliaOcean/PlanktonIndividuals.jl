
### To test the package

```
Pkg.develop(PackageSpec(path="AgentPhytModel_3D"))
Pkg.test("PhytoAgentModel")
```

_Note: this uses `AgentPhytModel_3D/samples/T_IR.csv`, `grid.jld`, & `uvw.jld`, and then compares results to `samples/testB1B2.csv`_

### To run the example

```
Pkg.develop(PackageSpec(path="AgentPhytModel_3D"))
using PhytoAgentModel
include("AgentPhytModel_3D/src/model_update.jl")
```

_Note: this runs the same example as `Pkg.test("PhytoAgentModel")` but interactively._

### To build and serve the docs

```
cd AgentPhytModel_3D/docs
julia make.jl
mkdocs build
mkdocs serve
```

_Note: this requires mkdocs since `format = Markdown()` is set in `docs/make.jl`_
