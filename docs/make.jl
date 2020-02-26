using Documenter
using DocumenterMarkdown
using PhytoAgentModel

makedocs(
    sitename = "PhytoAgentModel",
    format = Documenter.HTML(),
    modules = [PhytoAgentModel],
    pages = [
    "Home" => "index.md",
    "Various" => "various.md",
            ]
)

deploydocs(;
    repo="github.com/zhenwu0728/AgentPhytModel_3D/",
)
