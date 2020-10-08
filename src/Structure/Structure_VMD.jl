"""
    mutable struct Structure_VMD  <: Str

This contains information for `converter_vmd()` to assign all potential parameters.

"""
mutable struct Structure_VMD <: Str
    atom_name
    num_atom_types::Int64
    num_bond_types::Int64
    num_angle_types::Int64
    num_dihedral_types::Int64
    num_improper_types::Int64
    para_mass
    para_pair
    para_bond
    para_angle
    para_dihedral
    para_improper
end
