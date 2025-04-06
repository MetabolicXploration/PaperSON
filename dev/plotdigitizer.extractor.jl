src = """
x, y
0, 0.4512440175446377
0, 0.4805679832646803
0, 0.4176415219856269
0, -0.023017573541920272
0, 0.24444367219443228
0, 0.4827374321274982
0, 0.4739283983027269
0, 0.4541417643618093
0, 0.4881628635559996
0, 0.244180224618845
0, 0.5843996892465954
0, 0.5744086274888792
0, 0.2993746540044142
0, 0.45582332773121115
0, 0.6014585596301675
0, 0.352024982949823
0, 0.47702269378014395
0, 0.46488521827425977
0, 0.36523866607269617
0, 0.47117071557348067
0, 0.5182626978311471
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





















