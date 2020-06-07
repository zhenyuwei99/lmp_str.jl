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
    charge_vec = [0.4170, -0.8340]

    atom_type = [1, 2, 1]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
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
    charge_vec = [0.41, -0.82]

    atom_type = [1, 2, 1]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
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
    charge_vec = [0.4238, -0.8476]

    atom_type = [1, 2, 1]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
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
Tip4p()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip4p moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function Tip4p()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.15            # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    charge_vec = [0.52, 0.0, -1.04]

    atom_type = [1, 2, 1, 3]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    atom_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        0          0          0
        1/ratio[1] 1/ratio[2] 0
        2/ratio[1] 0          0
        1/ratio[1] 1/cell_vec[2]*dummy_len 0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 3)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    bond_mode[3] = Bond(0, 1, [2, 4])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 2)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    angle_mode[2] = Angle(2, 1, [1, 2, 4])
    num_angles = length(angle_mode)
    num_angle_types = 2

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end

"""
Tip4p_Ew()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip4p_Ew moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function Tip4p_Ew()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.125           # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    charge_vec = [0.52422, 0.0, -1.04844]

    atom_type = [1, 2, 1, 3]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    atom_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        0          0          0
        1/ratio[1] 1/ratio[2] 0
        2/ratio[1] 0          0
        1/ratio[1] 1/cell_vec[2]*dummy_len 0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 3)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    bond_mode[3] = Bond(0, 1, [2, 4])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 2)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    angle_mode[2] = Angle(2, 1, [1, 2, 4])
    num_angles = length(angle_mode)
    num_angle_types = 2

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end

"""
Tip4p_FQ()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip4p_FQ moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function Tip4p_FQ()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.15            # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    charge_vec = [0.63, 0.0, -1.26]

    atom_type = [1, 2, 1, 3]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    atom_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        0          0          0
        1/ratio[1] 1/ratio[2] 0
        2/ratio[1] 0          0
        1/ratio[1] 1/cell_vec[2]*dummy_len 0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 3)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    bond_mode[3] = Bond(0, 1, [2, 4])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 2)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    angle_mode[2] = Angle(2, 1, [1, 2, 4])
    num_angles = length(angle_mode)
    num_angle_types = 2

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end

"""
Tip4p_2005()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip4p_2005 moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function Tip4p_2005()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.1546          # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    charge_vec = [0.5564, 0.0, -1.1128]

    atom_type = [1, 2, 1, 3]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    atom_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        0          0          0
        1/ratio[1] 1/ratio[2] 0
        2/ratio[1] 0          0
        1/ratio[1] 1/cell_vec[2]*dummy_len 0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 3)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    bond_mode[3] = Bond(0, 1, [2, 4])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 2)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    angle_mode[2] = Angle(2, 1, [1, 2, 4])
    num_angles = length(angle_mode)
    num_angle_types = 2

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end


"""
Tip4p_Ice()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip4p_Ice moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function Tip4p_Ice()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.1577          # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    charge_vec = [0.5897, 0.0, -1.1794]

    atom_type = [1, 2, 1, 3]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    atom_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        0          0          0
        1/ratio[1] 1/ratio[2] 0
        2/ratio[1] 0          0
        1/ratio[1] 1/cell_vec[2]*dummy_len 0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 3)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    bond_mode[3] = Bond(0, 1, [2, 4])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 2)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    angle_mode[2] = Angle(2, 1, [1, 2, 4])
    num_angles = length(angle_mode)
    num_angle_types = 2

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end

"""
Tip5p()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip5p moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function Tip5p()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.70            # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    dummy_ang = 109.47 / 180 * pi
    charge_vec = [0.2410, 0.0, -0.2410]

    atom_type = [1, 2, 1, 3, 3]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    atom_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        0          0          0.5
        1/ratio[1] 1/ratio[2] 0.5
        2/ratio[1] 0          0.5
        1/ratio[1] 1/ratio[2]+1/cell_vec[2]*dummy_len*cos(dummy_ang/2) 0.5 + 1/cell_vec[3]*dummy_len*sin(dummy_ang/2)
        1/ratio[1] 1/ratio[2]+1/cell_vec[2]*dummy_len*cos(dummy_ang/2) 0.5 - 1/cell_vec[3]*dummy_len*sin(dummy_ang/2)
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 4)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    bond_mode[3] = Bond(0, 1, [2, 4])
    bond_mode[4] = Bond(0, 1, [2, 5])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 2)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    angle_mode[2] = Angle(2, 1, [2, 4, 5])
    num_angles = length(angle_mode)
    num_angle_types = 2

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end

"""
Tip5p_2018()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip5p_2018 moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function Tip5p_2018()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.70            # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    dummy_ang = 109.47 / 180 * pi
    charge_vec = [0.394137, -0.641114, -0.07358]

    atom_type = [1, 2, 1, 3, 3]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    atom_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        0          0          0.5
        1/ratio[1] 1/ratio[2] 0.5
        2/ratio[1] 0          0.5
        1/ratio[1] 1/ratio[2]+1/cell_vec[2]*dummy_len*cos(dummy_ang/2) 0.5 + 1/cell_vec[3]*dummy_len*sin(dummy_ang/2)
        1/ratio[1] 1/ratio[2]+1/cell_vec[2]*dummy_len*cos(dummy_ang/2) 0.5 - 1/cell_vec[3]*dummy_len*sin(dummy_ang/2)
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 4)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    bond_mode[3] = Bond(0, 1, [2, 4])
    bond_mode[4] = Bond(0, 1, [2, 5])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 2)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    angle_mode[2] = Angle(2, 1, [2, 4, 5])
    num_angles = length(angle_mode)
    num_angle_types = 2

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end

"""
Tip7p()

Do this will generate a Famliy_Wat type which contains all information needed to build a Tip7p moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function Tip7p()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.41            # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    dummy_ang = 109.47 / 180 * pi
    charge_vec = [0.58014, 0.11094, -0.17724, 0.45837]

    atom_type = [1, 2, 1, 3, 3, 4, 4]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    atom_mass = [1.00784, 15.9994, 0, 0]
    atom_name = split("H O Dummy01 Dummy02")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (atom_mass[1]*2+atom_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        0          0          0.5
        1/ratio[1] 1/ratio[2] 0.5
        2/ratio[1] 0          0.5
        1/ratio[1] 1/ratio[2]+1/cell_vec[2]*dummy_len*cos(dummy_ang/2) 0.5 + 1/cell_vec[3]*dummy_len*sin(dummy_ang/2)
        1/ratio[1] 1/ratio[2]+1/cell_vec[2]*dummy_len*cos(dummy_ang/2) 0.5 - 1/cell_vec[3]*dummy_len*sin(dummy_ang/2)
        1/ratio[1]/2 1/ratio[2]/2 0.5
        3/ratio[1]/2 1/ratio[2]/2 0.5
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 6)
    bond_mode[1] = Bond(0, 1, [1, 6])
    bond_mode[2] = Bond(0, 1, [3, 7])
    bond_mode[3] = Bond(0, 1, [2, 4])
    bond_mode[4] = Bond(0, 1, [2, 5])
    bond_mode[5] = Bond(0, 1, [2, 6])
    bond_mode[6] = Bond(0, 1, [2, 7])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 2)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    angle_mode[2] = Angle(2, 1, [2, 4, 5])
    num_angles = length(angle_mode)
    num_angle_types = 2

    Family_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2])
end