const Instance = Union{String, Dict}

function unique_nodes(bh_tuples)
    sort(collect(Set(Iterators.flatten((i, j) for (i, j, _) ∈ bh_tuples))))
end

function lattice(instance::Instance)
    if instance isa String
        bhi = CSV.File(instance, types = [Int, Int, Float64], header=0, comment = "#")
    else
        bhi = [(i, j, J) for ((i, j), J) ∈ instance]
    end
    bhg = LabelledGraph{MetaGraph}(unique_nodes(bhi))

    set_prop!.(Ref(bhg), vertices(bhg), :U, 0)

    for (i, j, v) ∈ bhg
        if i == j
            set_prop!(bhg, i, :U, v)
        else
            add_edge!(bhg, i, j) || throw(ArgumentError("Duplicate Egde ($i, $j)"))
            set_prop!(bhg, i, j, :J, v)
        end
    end
    bhg
end
