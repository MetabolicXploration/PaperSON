## ----------------------------------------------------------------------------
# BUNDLE from the ammonium-limited culture reported at
# Folsom, James Patrick, and Ross P. Carlson. “Physiological, Biomass Elemental Composition and Proteomic Analyses of Escherichia Coli Ammonium-Limited Chemostat Growth, and Comparison with Iron- and Glucose-Limited Chemostat Growth.” Microbiology (Reading, England) 161, no. 8 (August 2015): 1659–70. https://doi.org/10.1099/mic.0.000118.


## ----------------------------------------------------------------------------
# Sumpplementary file 00118s2.xlsx
const _S2 = Dict()
function _populate_s2()
    empty!(_S2)
    _S2["D"] = Dict(
        "unit" => "1/h",
        "val" => [0.1, 0.2, 0.3, 0.4],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    _S2["YX_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.140, 0.196, 0.186, 0.260],
        "err" => [0.032, 0.015, 0.063, 0.049],
    )
    _S2["Ypyr_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.003, 0.0, 0.0, 0.0],
        "err" => [0.002, 0.0, 0.0, 0.0],
    )
    _S2["Ysucc_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.003, 0.000, 0.000, 0.000],
        "err" => [0.003, 0.00, 0.0, 0.0],
    )
    _S2["Ylac_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.078, 0.047, 0.026, 0.026],
        "err" => [0.02, 0.005, 0.01, 0.005],
    )
    _S2["Yform_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.007, 0.0002, 0.0, 0.002],
        "err" => [0.004, 0.00009, 0.000, 0.001],
    )
    _S2["Yac_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.310, 0.262, 0.183, 0.224],
        "err" => [0.073, 0.019, 0.065, 0.043],
    )
    _S2["Yco2_glc"] = Dict(
        "unit" => "mol Co2/Cmmol glc",
        "val" => [0.460, 0.495, 0.605, 0.488],
        "err" => [0.134, 0.04, 0.138, 0.098],
    )
    _S2["Yo2_glc"] = Dict(
        "unit" => "mol o2/Cmmol glc",
        "val" => [0.459, 0.487, 0.598, 0.479],
        "err" => [0.133, 0.039, 0.136, 0.097],
    )

    _S2["qglc"] = Dict(
        "unit" => "Cmmol/gCDW-h",
        "val" => [27.8, 39.7, 62.8, 59.9],
        "err" => [6.41, 3.1, 21.3, 11.2],
    )
    _S2["qO2"] = Dict(
        "unit" => "mmol o2/gCDW-h",
        "val" => [12.77, 19.33, 37.58, 28.68],
        "err" => [6.65, 3.05, 21.28, 11.18],
    )
    return _S2
end

## ----------------------------------------------------------------------------
# Sumpplementary file 000118s4.xlsx
const _S4 = Dict()
function _populate_s4()
    empty!(_S4)
    _S4["D"] = Dict(
        "unit" => "1/h",
        "val" => [0.1, 0.2, 0.3, 0.4],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    _S4["Y_X/N"] = Dict(
        "unit" => "g/g_N",
        "val" => [9.22, 9.84, 8.82, 8.81],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    _S4["Y_X/glc"] = Dict(
        "unit" => "g/g_glc",
        "val" => [0.12, 0.17, 0.16, 0.22],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    _S4["qglc"] = Dict(
        "unit" => "mmol/ gCDW h",
        "val" => [4.6, 6.6, 10.5, 10.0],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    _S4["qac"] = Dict(
        "unit" => "mmol/ gCDW h",
        "val" => [4.3, 5.2, 5.7, 6.71],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    _S4["qpyr"] = Dict(
        "unit" => "mmol/ gCDW h",
        "val" => [0.03, 0.0, 0.0, 0.0],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    _S4["qo2"] = Dict(
        "unit" => "mmol/ gCDW h",
        "val" => [12.8, 19.33, 37.6, 28.7],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    _S4["qco2"] = Dict(
        "unit" => "mmol/ gCDW h",
        "val" => [12.8, 19.6, 38.0, 29.2],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    return _S4
end

## ----------------------------------------------------------------------------
#  M9 minima medium, all in mM (Ammonium-limited)
MEDIUM = Dict()
function _populate_medium()
    empty!(MEDIUM)
    MEDIUM["nh4"] = 1.31 # mM
    MEDIUM["fe"] = 18.7 # mM
    MEDIUM["glc"] = 18.7 # mM
    return MEDIUM
end

## ----------------------------------------------------------------------------
# Bundle, converting
const BUNDLE = Dict()
function _populate_bundle()

    empty!(BUNDLE)
    _populate_s2()
    _populate_s4()
    _populate_medium()

    # flxs
    BUNDLE["D"] = deepcopy(_S2["D"])
    BUNDLE["μ"] = deepcopy(_S2["D"])
    BUNDLE["uGLC"] = Dict(
        # q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "glucose uptake rate",
        "unit" => "mmol/gCDW h",
        "val" => -_S2["uglc"]["val"] ./ 6,
        "err" => _S2["uglc"]["err"] ./ 6,
    )
    BUNDLE["uO2"] = Dict(
        # q[mmol/gCDW-h]
        "name" => "oxigen uptake rate",
        "unit" => "mmol/gCDW h",
        "val" => -_S2["qO2"]["val"],
        "err" => _S2["qO2"]["err"]
    )
    
    # yields
    # qglc[Cmmol/gCDW-h]
    qglc = _S2["uglc"]["val"]

    BUNDLE["uPYR"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "pyruvate production rate",
        "unit" => "mmol/gCDW h",
        "val" => _S2["Ypyr_glc"]["val"] .* qglc ./ 3,
        "err" => _S2["Ypyr_glc"]["err"] .* qglc ./ 3,
    )
    BUNDLE["uSUCC"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "succinate production rate",
        "unit" => "mmol/gCDW h",
        "val" => _S2["Ysucc_glc"]["val"] .* qglc ./ 4,
        "err" => _S2["Ysucc_glc"]["err"] .* qglc ./ 4,
    )
    BUNDLE["uLAC"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "lactate production rate",
        "unit" => "mmol/gCDW h",
        "val" => _S2["Ylac_glc"]["val"] .* qglc ./ 3,
        "err" => _S2["Ylac_glc"]["err"] .* qglc ./ 3,
    )
    BUNDLE["uFORM"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "formate production rate",
        "unit" => "mmol/gCDW h",
        "val" => _S2["Yform_glc"]["val"] .* qglc ./ 1,
        "err" => _S2["Yform_glc"]["err"] .* qglc ./ 1,
    )
    BUNDLE["uAC"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "acetate production rate",
        "unit" => "mmol/gCDW h",
        "val" => _S2["Yac_glc"]["val"] .* qglc ./ 2,
        "err" => _S2["Yac_glc"]["err"] .* qglc ./ 2,
    )
    BUNDLE["uCO2"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "acetate production rate",
        "unit" => "mmol/gCDW h",
        "val" => _S2["Yco2_glc"]["val"] .* qglc ./ 1,
        "err" => _S2["Yco2_glc"]["err"] .* qglc ./ 1,
    )

    # There negigeble glucose present in the medium at steady state 
    # so we can compute X = c * D/qglc
    BUNDLE["X"] = Dict(
        # c[mmol/L] * D[1/h] / q[mmol/gCDW h] = X[gCDW/L]
        "name" => "cell concentration",
        "unit" => "gCDW/ L",
        "val" => abs.(MEDIUM["glc"] .* BUNDLE["D"]["val"] ./ BUNDLE["uGLC"]["val"]),
        "err" => zero(BUNDLE["uGLC"]["err"]) #TODO: fix this
    )

    # xi
    BUNDLE["xi"] = Dict(
        "name" => "dilution specific cell concentration",
        "unit" => "gCDW/ L h",
        "val" => BUNDLE["X"]["val"] ./ BUNDLE["D"]["val"],
        "err" => BUNDLE["X"]["err"] #TODO: fix this
    )

    # medium
    BUNDLE["cGLC"] = Dict(
        "name" => "glucose feed concentration",
        "unit" => "gCDW/ L h",
        "val" => fill(MEDIUM["glc"], 4),
        "err" => zeros(4)
    )

    BUNDLE
end

## ------------------------------------------------------------------
# API
const EXPS = 1:4
const msd_mets = ["GLC", "PYR", "SUCC", "LAC", "FORM", "AC"]

_get_val(id, dk) = BUNDLE[string(id)][dk]
_get_val(id, dk, exp::Int) = BUNDLE[string(id)][dk][exp]
function _get_val(id, dk, D::Float64)
    exp = findfirst(BUNDLE["D"][dk] .== D)
    isnothing(exp) && error("No experiment with D = $D")
    _get_val(id, dk, exp)
end
_get_val(id, dk, ref, dflt) =
    try; _get_val(id, dk, ref); catch err; dflt end

for fun in [:val, :err]
    dk = string(fun)
    @eval $fun(id) = _get_val(id, $dk)
    @eval $fun(id, ref) = _get_val(id, $dk, ref)
    @eval $fun(id, ref, dflt) = _get_val(id, $dk, ref, dflt)

    for p in [:u, :c]
        pstr = string(p)
        pfun = Symbol(p, fun)
        @eval $pfun(id, args...) = $fun(string($pstr, id), args...)
    end
end

name(id) = BUNDLE[string(id)]["name"]
unit(id) = BUNDLE[string(id)]["unit"]

ciD_X(id) = [cval(id, exp, 0.0) * val(:D, exp) / val(:X, exp) for exp in EXPS]
ciD_X(id, exp) = cval(id, exp, 0.0) * val(:D, exp) / val(:X, exp)