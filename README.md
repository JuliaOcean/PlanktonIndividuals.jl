
### To run the test case

```
;cd ./AgentPhytModel_3D/
include("src/model_update.jl")
```

### To test the package

```
]dev ./AgentPhytModel_3D/
]test PhytoAgentModel
```

### To build and serve the docs

(Requires mkdocs)

```
cd AgentPhytModel_3D/docs
julia make.jl
mkdocs build
mkdocs serve
```

