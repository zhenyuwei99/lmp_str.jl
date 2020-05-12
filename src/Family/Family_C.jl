"""
    mutable struct Family_C  <: Str

This contains information for model constructing of C and its compounds.

# Supported list:
- Graphene
- Graphene_Ort
Notice:Case sensetive
"""
mutable struct Family_C <: Str
    atom_vec::Matrix{Float64}
    cell_vec::Matrix{Float64}
    atom_type
    atom_name
    atom_charge
    atom_mass
    num_atoms
    num_atom_types
    vec_type_id
end

"""
    Graphene()

Do this will generate a Family_C type which contains all information needed to build a Graphene moedel.
"""
function Graphene()
    atom_vec = [
        1/3 1/3 0
        2/3 2/3 0
    ]
    cell_vec = [
        2.46 0                0
        1.23 2.46 * sin(pi/3) 0
        0    0                0.35
    ]
    atom_type = [1 1]
    atom_charge = 0 .* atom_type
    atom_mass = [12]
    atom_name = split("C")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_C(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, [1])
end

"""
    Graphene_Ort()

Do this will generate a Famliy_C type which contains all information needed to build a Graphene moedel with orthogonal unit cell.
"""
function Graphene_Ort()
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
    atom_mass = [12]
    atom_name = split("C")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_C(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, [1])
end