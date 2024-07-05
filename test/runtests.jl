using Ferrite
using Tensors
using Test
using Logging
using ForwardDiff
import SHA
using Random
using LinearAlgebra
using SparseArrays

const HAS_EXTENSIONS = isdefined(Base, :get_extension)

# https://github.com/JuliaLang/julia/pull/47749
const MODULE_CAN_BE_TYPE_PARAMETER = VERSION >= v"1.10.0-DEV.90"

if HAS_EXTENSIONS && MODULE_CAN_BE_TYPE_PARAMETER
    import Metis
end

const RUN_JET_TESTS = VERSION >= v"1.9"

if RUN_JET_TESTS
    using JET: @test_call
else
    # Just eat the macro on incompatible versions
    macro test_call(args...)
        nothing
    end
end

#
## Unit tests
#include("test_interfacevalues.jl")
#@test all(x -> isdefined(Ferrite, x), names(Ferrite))  # Test that all exported symbols are defined
#
## Integration tests
#include("integration/test_simple_scalar_convergence.jl")

include("test_p4est.jl")
include("test_p4est_example.jl")
