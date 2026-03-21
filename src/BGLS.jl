module BGLS

using LinearAlgebra

include("types.jl")
include("autofreq.jl")
include("bgls.jl")
include("sbgls.jl")

export bgls, bgls_auto, sbgls, sbgls_auto, autofrequencies
export BGLSResult, SBGLSResult
export frequencies, periods, power, best_frequency, best_period
export obs_times

end # module BGLS
