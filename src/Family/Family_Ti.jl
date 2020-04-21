mutable struct Family_Ti <: Str
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

function TiO2()
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
    atom_mass = [47.867, 15.999] # Ti: 47.867, O: 15.999. Unit: g/mol
    atom_name = split("Ti O")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Ti(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, [1, 2])
end