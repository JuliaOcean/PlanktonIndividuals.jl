
### To run the test case

```
;cd ./AgentPhytModel_3D/
include("src/model_update.jl")
```

(this uses `samples/T_IR.csv`, `grid.jld`, and `uvw.jld`)

### To test the package

```
]dev ./AgentPhytModel_3D/
]test PhytoAgentModel
```

(for now this just issues a print statement...)

### To build and serve the docs

```
cd AgentPhytModel_3D/docs
julia make.jl
mkdocs build
mkdocs serve
```

(this requires mkdocs since `format = Markdown()` is set in `docs/make.jl`)
