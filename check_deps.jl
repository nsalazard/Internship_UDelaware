for mod in [
    :BenchmarkTools, :LinearAlgebra, :Plots,
    :DifferentialEquations, :Tullio,
    :Base, :DelimitedFiles, :LaTeXStrings,
    :StaticArrays, :Serialization
]
    try
        @eval using $(mod)
        println("✅ loaded ", mod)
    catch e
        println("❌ failed to load ", mod, " — ", e)
    end
end

println("Threads: ", Threads.nthreads())
println("Julia Version: ", VERSION)