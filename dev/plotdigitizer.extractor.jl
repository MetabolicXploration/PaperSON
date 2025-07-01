src = """
x, y
0.1549715046481629, 1484.8630355137261
0.228600191219402, 1421.3781754583333
0.34709789043775374, 1253.4140076849862
0.4285512455493536, 1224.612537740435
0.528625675967805, 1005.1203130616829
0.6734279667850068, 855.2858028910896
0.7740451566022966, 736.6291295626114
0.924836169934932, 592.0504348182405
1.0401504025546107, 512.4674282116855
1.208996141967451, 400.78144353323114
1.45069693761496, 264.8965059946991
1.674735719807502, 183.43416075143745
2.066180681838836, 80.51550500939669
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





















