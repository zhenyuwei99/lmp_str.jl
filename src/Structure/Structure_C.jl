"""
    mutable struct Structure_C  <: Str

This contains information for model constructing of C and its compounds.

# Supported list:
- structure_graphene()
- structure_graphene_ort()
"""
mutable struct Structure_C <: Str
    atom_vec::Matrix{Float64}
    cell_vec::Matrix{Float64}
    atom_type
    atom_name
    atom_charge
    para_mass
    num_atoms
    num_atom_types
    vec_type_id
end

"""
    structure_graphene()

Do this will generate a Structure_C type which contains all information needed to build a Graphene moedel.
"""
function structure_graphene()
    atom_vec = [
        1/3 1/3 0
        2/3 2/3 0
    ]
    cell_vec = [
        2.46 0                0
        1.23 2.46 * sin(pi/3) 0
        0    0                3.35
    ]
    atom_type = [1 1]
    atom_charge = 0 .* atom_type
    para_mass = [12]
    atom_name = split("C")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Structure_C(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, num_atom_types, [1])
end

"""
    structure_graphene_ort()

Do this will generate a Famliy_C type which contains all information needed to build a Graphene moedel with orthogonal unit cell.
"""
function structure_graphene_ort()
    atom_vec = [
        0   0   0
        0   1/3 0
        1/2 1/2 0
        1/2 5/6 0
    ]
    cell_vec = [
        2.46 0                0
        0    3.69 / cos(pi/6) 0
        0    0                0.35
    ]
    atom_type = [1 1 1 1]
    atom_charge = 0 .* atom_type
    para_mass = [12]
    atom_name = split("C")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Structure_C(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, num_atom_types, [1])
end