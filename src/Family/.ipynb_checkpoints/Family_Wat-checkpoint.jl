"""
    mutable struct Family_Wat  <: Str

This contains information for model constructing of the famliy of water models.

# Supported list:
 - Tip3p
 - SPC
 - SPCE
Notice:Case sensetive
"""
mutable struct Family_Wat <: Str
    atom_vec::Matrix{Float64}
    cell_vec::Matrix{Float64}
    atom_type
    atom_name
    atom_charge
    atom_mass
    num_atoms
    num_atom_types
    bond_mode
    num_bonds
    num_bond_types
    angle_mode
    num_angles
    num_angle_types
    vec_type_id
end

"""
    Tip3p()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip3p moedel.
"""
function Tip3p()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.25;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572;          # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details

    atom_type = [1, 2, 1]
    atom_charge = [0.41, -0.82, 0.41]
    atom_mass = [1.00784, 15.9994]
    atom_name = split("H O")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    atom_vec = [
        0          0          0
        1/ratio[1] 1/ratio[2] 0
        2/ratio[1] 0          0
    ]
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 2)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    num_angles = length(angle_mode)
    num_angle_types = 1

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end

"""
    SPC()

Do this will generate a Famliy_Wat type which contains all information needed to build a SPC moedel.
"""
function SPC()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 109.47;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 1.0;          # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details

    atom_type = [1, 2, 1]
    atom_charge = [0.41, -0.82, 0.41]
    atom_mass = [1.00784, 15.9994]
    atom_name = split("H O")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    atom_vec = [
        0          0          0
        1/ratio[1] 1/ratio[2] 0
        2/ratio[1] 0          0
    ]
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 2)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    num_angles = length(angle_mode)
    num_angle_types = 1

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end

"""
    SPCE()

Do this will generate a Famliy_Wat type which contains all information needed to build a SPCE moedel.
"""
function SPCE()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 109.47;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 1.0;          # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details

    atom_type = [1, 2, 1]
    atom_charge = [0.4238, -0.8476, 0.4238]
    atom_mass = [1.00784, 15.9994]
    atom_name = split("H O")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    atom_vec = [
        0          0          0
        1/ratio[1] 1/ratio[2] 0
        2/ratio[1] 0          0
    ]
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 2)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    num_angles = length(angle_mode)
    num_angle_types = 1

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end