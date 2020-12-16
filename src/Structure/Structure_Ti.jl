"""
    mutable struct Structure_Ti  <: Str

This contains information for model constructing of Ti and its compounds.

# Supported list:
- structure_tio2_anatase()
- structure_tio2_rutile()
"""
mutable struct Structure_Ti <: Str
    atom_vec::Matrix{Float64}
    cell_vec::Matrix{Float64}
    atom_type
    atom_name
    atom_charge
    para_mass
    num_atoms
    num_atom_types
    vec_type_id
    res_name 
end

"""
    function structure_tio2_anatase()
Do this will generate a Structure_Ti type which contains all information needed to build a TiO2 moedel in structure of Anatase.
"""
function structure_tio2_anatase()
    atom_vec = [
        0.500	0.000	0.000
        0.500	0.500	0.250
        0.000	0.500	0.500
        0.000	0.000	0.750
        0.500	0.000	0.208
        0.500	0.500	0.458
        0.500	0.500	0.042
        0.500	0.000	0.792
        0.000	0.500	0.292
        0.000	0.000	0.958
        0.000	0.500	0.708
        0.000	0.000	0.542
    ]
    cell_vec = [
        3.7842      0           0
        0           3.7842      0
        0           0           9.5146
    ] 
    atom_type = vcat([1 for n = 1:4], [2 for n = 5:12])
    atom_charge = 0 .* atom_type
    para_mass = [47.867, 15.999] # Ti: 47.867, O: 15.999. Unit: g/mol
    atom_name = split("Ti O")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    res_name = "TIO"
    Structure_Ti(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, num_atom_types, [1, 2], res_name)
end

"""
    function structure_tio2_rutile()
Do this will generate a Structure_Ti type which contains all information needed to build a TiO2 moedel in structure of Rutile.
"""
function structure_tio2_rutile()
    atom_vec = [
        0.804921  0.804921  0.0     
        0.305029  0.305029  0.499831
        0.610059  0.610059  0.499831
        0.110168  0.499891  0.0     
        0.499891  0.110168  0.0     
        0.0       0.0       0.499831
    ]
    cell_vec = [
        4.593  0.0    0.0  
        0.0    4.593  0.0  
        0.0    0.0    2.959
    ] 
    atom_type = vcat([1 for n = 1:2], [2 for n = 1:4])
    atom_charge = 0 .* atom_type
    para_mass = [47.867, 15.999] # Ti: 47.867, O: 15.999. Unit: g/mol
    atom_name = split("Ti O")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    res_name = "TIO"
    Structure_Ti(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, num_atom_types, [1, 2], res_name)
end