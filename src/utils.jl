
## --. -. - .- -.- - -- -- - - ..- .-. -.-. -
export walkdict
function walkdict(f::Function, root_dict::AbstractDict; 
        path = []
    )
    for (key, dat) in root_dict
        push!(path, key)

        # function
        f(root_dict, path) === :break && return :break

        # recursive
        dat isa AbstractDict || continue
        walkdict(f, dat; path) === :break && return :break
    end
    empty!(path) # end of path
    return nothing
end

## --. -. - .- -.- - -- -- - - ..- .-. -.-. -
export findpath
function findpath(f::Function, root00::String = pkgdir(PaperSON))
    for (root, dirs, files) in walkdir(root0)
        for dir in dirs
            path = joinpath(root, dir)
            f(path) === true && return path
        end
        for file in files
            path = joinpath(root, file)
            f(path) === true && return path
        end
    end
    return nothing
end
findpath(bname::String, root0::String = pkgdir(PaperSON)) = 
    findpath(root0) do path
        basename(path) == bname
    end

## --. -. - .- -.- - -- -- - - ..- .-. -.-. -
nothing