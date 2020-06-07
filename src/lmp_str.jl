module lmp_str

import Dates
using Printf
using Statistics
import LinearAlgebra.norm

abstract type Str end
abstract type Data end
abstract type Unit end

include("Constants.jl")
#=
export Const_k_b, Const_n_a, Const_density_wat,
    Const_kcal2j, Const_kcalm2j, Const_kcalm2t,
    Const_g2kg, Const_kg2g, Const_gm2g, Const_gm2kg,
    Const_an2m, Const_an2nm, Const_nm2m, Const_cm2an, Const_dm2m, Const_m2dm,
    Const_fs2s, Const_ps2s, Const_ns2s
=#

include("./Family/Family_Basic.jl")
export Family_Basic, transform, SC, BCC, FCC, DC

include("./Family/Family_Wat.jl")
export Family_Wat, Tip3p, SPC, SPCE, Tip4p, Tip4p_2005, Tip4p_Ew, Tip4p_FQ, Tip4p_Ice, Tip5p, Tip5p_2018

include("./Family/Family_Si.jl")
export Family_Si, Si, Si3N4, Si3N4_Ort, SiO2

include("./Family/Family_C.jl")
export Family_C, Graphene, Graphene_Ort

include("./Family/Family_Ti.jl")
export Family_Ti, TiO2_Anatase, TiO2_Rutile

include("./Family/Family_Ca.jl")
export Family_Ca, CaCO3

include("./Family/Family_Ion.jl")
export Family_Ion

include("Data.jl")
export Atom, Bond, Angle, Data_Cell, Data_Basic, Data_Unit, Data_Sum

include("useful_funcs.jl")
export max, min, get_data, add!, change!, diag, conv, dist, norm_vec, rot_mat, central_point_atom, central_point_box, copy_array

include("genr.jl")
export genr, genr_cell, genr_atom, genr_bond, genr_angle

include("modify.jl")
export move, find, select, remove!, remove_else!

include("genr_ions.jl")
export genr_ions

include("combine.jl")
export cat_data, sort_data!

include("write.jl")
export write_data, write_info, write_xyz

end # module
