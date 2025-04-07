src = """
x, y
14.898731469490958, 2.258117755502685
16.253273792388036, 5.564523918322426
17.41415105847497, 8.462771935231139
18.610735130946814, 11.701485518076968
20.564350574639732, 14.504870966125791
22.18610696096236, 16.959017690727585
23.19697090645429, 18.35948622064266
23.963254675138497, 17.557514893547634
0.7886666621126173, 2.775993146537611
3.225911631270316, 8.271738886671182
3.905565276282802, 9.934869681053854
4.314416787267658, 11.048969962279045
6.123461324896631, 13.505716451449409
8.877791494141377, 17.03889833283248
9.610508534451547, 16.304795181047126
"""

# extract
let
    # clear lines
    lines = split(src, "\n"; keepempty = false)
    filter!(lines) do line
        contains(line, "x, y") && return false
        contains(line, ",") || return false
        return true
    end

    xs = Float64[]
    ys = Float64[]
    for line in lines
        xstr, ystr = split(line, ",")
        x = parse(Float64, xstr)
        y = parse(Float64, ystr)
        push!(xs, x)
        push!(ys, y)
    end
    # comma_lines = String[]

    # Copy
    clipboard(string(
        join(repr.(xs), ", "), 
        "\n\n", 
        join(repr.(ys), ", "), 
    ))
end





















