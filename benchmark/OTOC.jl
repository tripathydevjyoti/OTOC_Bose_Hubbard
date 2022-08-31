using .OTOC_example
using BenchmarkTools

function bench(
    L ::Int, 
    N ::Int
)

    t_vals = range(0,0.2,40)
    vals = []

    for t in t_vals
        push!(vals, OTOC_lattice([1,1,1,1,1,1], 6, 1, t, L, N))
    end
    
    return vals
end

L=6
N=6

@time vals = bench(L, N);
@time vals = bench(L, N);

otoc_bench = @benchmark bench(L,N)

nothing