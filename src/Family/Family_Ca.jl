"""
    mutable struct Family_Ca  <: Str

This contains information for model constructing of C and its compounds.

# Supported list:
- CaCO3
Notice:Case sensetive
"""
mutable struct Family_Ca <: Str
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
    CaCO3()

Do this will generate a Family_Ca type which contains all information needed to build a CaCO3 moedel.
"""
function CaCO3()
    atom_vec = [
		0.772116  0.360801   0.0909051
		0.32908   0.721614   0.454547 
		0.329085  0.0        0.818179 
		0.329085  0.0        0.272726 
		0.772116  0.360801   0.636358 
		0.32908   0.721614   1.0      
		0.329085  0.0        0.0      
		0.772116  0.360801   0.363631 
		0.32908   0.721614   0.727274 
		0.32908   0.721614   0.181821 
		0.329085  0.0        0.545453 
		0.772116  0.360801   0.909084 
		1.0       0.360801   0.0909051
		0.556973  0.721614   0.454547 
		0.556969  0.0        0.818179 
		0.658165  0.0824151  0.0909051
		0.215143  0.443217   0.454547 
		0.658169  0.804029   0.818179 
		0.658169  0.639199   0.0909051
		0.215138  1.0        0.454547 
		0.215143  0.278386   0.818179 
		0.987254  0.0        0.272726 
		0.544223  0.360801   0.636358 
		0.101196  0.721614   1.0      
		0.443027  0.278386   0.272726 
		0.886054  0.639199   0.636358 
		0.443031  1.0        1.0      
		0.0       0.804029   0.272726 
		0.886058  0.0824151  0.636358 
		0.443027  0.443217   1.0     
    ]
    cell_vec = [
        5.69943  0.0       0.0   
		0.0      4.04044   0.0   
		0.0      0.0      15.8088
    ]
    atom_type = vcat([1 for n=1:6], [2 for n=1:6], [3 for n=1:18])
    atom_charge = 0 .* atom_type
    atom_mass = [40.078, 12.011, 15.999]
    atom_name = split("Ca C O")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Ca(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, [1])
end