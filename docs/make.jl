using Documenter
using PhytoAgentModel

makedocs(
    sitename = "PhytoAgentModel",
    format = Documenter.HTML(),
    modules = [PhytoAgentModel]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
