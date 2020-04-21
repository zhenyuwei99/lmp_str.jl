module lmp_str

import Dates
using Printf

abstract type Str end
abstract type Data end
abstract type Unit end

include("Constants.jl")

include("./Family/Family_Basic.jl")
export Family_Basic, transform, SC, BCC, FCC, DC

include("./Family/Family_Wat.jl")
export Family_Wat, Tip3p, SPC, SPCE

include("./Family/Family_Si.jl")
export Family_Si, Si, Si3N4, Si3N4_Ort, SiO2

include("./Family/Family_C.jl")
export Family_C, Graphene, Graphene_Ort

include("./Family/Family_Ti.jl")
export Family_Ti, TiO2

include("./Family/Family_Ion.jl")
export Family_Ion

include("Data.jl")
export Atom, Bond, Angle, Data_Cell, Data_Basic, Data_Unit, Data_Sum

include("useful_funcs.jl")
export max, min, add!, change!, diag, conv, dist


include("genr.jl")
export genr, genr_cell, genr_atom, genr_bond, genr_angle

include("modify.jl")
export move, find, select, remove!

include("addions.jl")
export addions

include("combine.jl")
export cat_data, sort_data!

include("write.jl")
export write_data, write_info, write_xyz

end # module
