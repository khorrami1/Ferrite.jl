using JuAFEM
using Tensors
using Test
using ForwardDiff
import SHA
using Random
using LinearAlgebra

include("test_utils.jl")
include("test_interpolations.jl")
include("test_cellvalues.jl")
include("test_facevalues.jl")
include("test_quadrules.jl")
include("test_assemble.jl")
include("test_dofs.jl")
include("test_constraints.jl")
include("test_grid_dofhandler_vtk.jl")
# include("test_notebooks.jl")
include("test_examples.jl")

# Build the docs
if get(ENV, "TRAVIS_OS_NAME", "") == "linux" && get(ENV, "TRAVIS_JULIA_VERSION", "") == "1.0"
    include("../docs/make.jl")
end
