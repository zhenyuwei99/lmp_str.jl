module lmp_str

    import Dates
    using Printf
    using Statistics
    import LinearAlgebra.norm

    abstract type Str end
    abstract type Potential end
    abstract type Potential_Unit end
    abstract type Data end
    abstract type Unit end
    
    # Main Files
    include("Data.jl")
    export Atom, Bond, Angle, Dihedral, Improper
    export Data_Cell, Data_Basic, Data_Unit, Data_Sum

    include("useful_funcs.jl")
    export max, min, get_data, diag, conv, dist, norm_vec, rot_mat
    export add!, change!, central_point_atom, central_point_box, copy_array

    include("Constants.jl")
    export const_list
    #=
    export Const_k_b, Const_n_a, Const_density_wat,
        Const_kcal2j, Const_kcalm2j, Const_kcalm2t,
        Const_g2kg, Const_kg2g, Const_gm2g, Const_gm2kg,
        Const_an2m, Const_an2nm, Const_nm2m, Const_cm2an, Const_dm2m, Const_m2dm,
        Const_fs2s, Const_ps2s, Const_ns2s
    =#

    # Structure Files
    include("./Structure/Structure_Basic.jl")
    export Structure_Basic, transform
    export structure_sc, structure_bcc, structure_fcc, structure_dc

    include("./Structure/Structure_Wat.jl")
    export Structure_Wat
    export structure_tip3p, structure_spc, structure_spce
    export structure_tip4p_cut, structure_tip4p_long
    export structure_tip4p_2005, structure_tip4p_ice
    export structure_tip4p_ew, structure_tip4p_fq
    export structure_tip5p, structure_tip4p_2005, structure_tip5p_2018
    export structure_tip7p

    include("./Structure/Structure_Si.jl")
    export Structure_Si
    export structure_si, structure_sio2
    export structure_si3n4, structure_si3n4_ort

    include("./Structure/Structure_C.jl")
    export Structure_C
    export structure_graphene, structure_graphene_ort

    include("./Structure/Structure_Ti.jl")
    export Structure_Ti
    export structure_tio2_anatase, structure_tio2_rutile

    include("./Structure/Structure_Ca.jl")
    export Structure_Ca
    export structure_caco3

    include("./Structure/Structure_Ion.jl")
    export Structure_Ion

    include("./Structure/Structure_VMD.jl")
    export Structure_VMD

    # generator Files
    include("./generator/genr.jl")
    export genr, genr_cell, genr_atom
    export genr_bond, genr_angle, genr_dihedral, genr_improper

    include("./generator/modify.jl")
    export move, find, select, remove!, remove_else!

    include("./generator/genr_ions.jl")
    export genr_ions

    include("./generator/combine.jl")
    export cat_data, sort_data!

    # Potential Files
    include("./Potential/Potential_Charmm36.jl")
    export Potential_Charmm36, Potential_Charmm36_Unit, potential_charmm36

    # converter Files
    include("./converter/useful_funcs.jl")

    include("converter/converter_vmd.jl")
    export converter_vmd

    # writeer Files
    include("./Writer/write_data.jl")
    export write_data, write_info

    include("./Writer/write_xyz.jl")
    export write_xyz

    include("./Writer/write_pdb.jl")

end # module
