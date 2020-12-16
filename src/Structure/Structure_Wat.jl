"""
    mutable struct Structure_Wat  <: Str

This contains information for model constructing of a series of water models.

# Supported list:
 - structure_tip3p()
 - structure_spc()
 - structure_spce()
 - structure_tip4p_cut()
 - structure_tip4p_long()
 - structure_tip4p_ew()
 - structure_tip4p_fq()
 - structure_tip4p_2005()
 - structure_tip4p_ice()
 - structure_tip5p()
 - structure_tip4p_2005()
 - structure_tip5p_2018()
 - structure_tip7p()

# Parameters source:
- https://en.wikipedia.org/wiki/Water_model
- http://www1.lsbu.ac.uk/water/water_models.html
- https://lammps.sandia.gov/doc/Howto_tip3p.html
- https://lammps.sandia.gov/doc/Howto_tip4p.html
- https://lammps.sandia.gov/doc/Howto_spc.html
"""
mutable struct Structure_Wat <: Str
    atom_vec::Matrix{Float64}
    cell_vec::Matrix{Float64}
    atom_type
    atom_name
    atom_charge
    para_mass
    num_atoms
    num_atom_types
    bond_mode
    num_bonds
    num_bond_types
    angle_mode
    num_angles
    num_angle_types
    vec_type_id
    res_name
end

"""
    structure_tip3p()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip3p moedel.
"""
function structure_tip3p()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.25;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572;          # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    charge_vec = [0.4170, -0.8340]

    atom_type = [1, 2, 1]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    para_mass = [1.00784, 15.9994]
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
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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
    res_name = "TIP3"

    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_spc()

This will generate a `Structure_Wat` type which contains all information needed to build a SPC moedel.
"""
function structure_spc()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 109.47;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 1.0;          # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    charge_vec = [0.41, -0.82]

    atom_type = [1, 2, 1]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    para_mass = [1.00784, 15.9994]
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
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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
    res_name="SPC"

    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_spce()

This will generate a `Structure_Wat` type which contains all information needed to build a SPCE moedel.
"""
function structure_spce()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 109.47;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 1.0;          # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    charge_vec = [0.4238, -0.8476]

    atom_type = [1, 2, 1]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    para_mass = [1.00784, 15.9994]
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
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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
    res_name = "SPCE"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_tip4p_cut()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip4p moedel used in lammps with lj/cut/tip4p/cut. Detail info: https://lammps.sandia.gov/doc/Howto_tip4p.html

# Model info
- O mass = 15.9994
- H mass = 1.008
- O charge = -1.040
- H charge = 0.520
- ``r_0`` of OH bond = 0.9572
- ``θ`` of HOH angle = 104.52
- OM distance = 0.15
- LJ ``ε`` of O-O = 0.1550
- LJ `` σ`` of O-O = 3.1536
- LJ ``ε``, `` σ`` of OH, HH = 0.0
- Coulomb cutoff = 8.5
"""
function structure_tip4p_cut()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    charge_vec = [-1.040, 0.520]

    atom_type = [1, 2, 2]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    para_mass = [15.9994, 1.008]
    atom_name = split("O H")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        1/ratio[1] 1/ratio[2] 0
        0          0          0
        2/ratio[1] 0          0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 2)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [1, 3])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, [2, 1, 3])
    num_angles = length(angle_mode)
    num_angle_types = 1
    res_name = "TIP4P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_tip4p_long()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip4p moedel used in lammps with lj/cut/tip4p/long. Detail info: https://lammps.sandia.gov/doc/Howto_tip4p.html

# Model info
- O mass = 15.9994
- H mass = 1.008
- O charge = -1.0484
- H charge = 0.5242
- ``r_0`` of OH bond = 0.9572
- ``θ`` of HOH angle = 104.52
- OM distance = 0.1250
- LJ ``ε`` of O-O = 0.16275
- LJ `` σ`` of O-O = 3.16435
- LJ ``ε``, `` σ`` of OH, HH = 0.0
- Coulomb cutoff = 8.5
"""
function structure_tip4p_long()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    charge_vec = [-1.0484, 0.5242]

    atom_type = [1, 2, 2]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    para_mass = [15.9994, 1.008]
    atom_name = split("O H")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        1/ratio[1] 1/ratio[2] 0
        0          0          0
        2/ratio[1] 0          0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 2)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [1, 3])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, [2, 1, 3])
    num_angles = length(angle_mode)
    num_angle_types = 1

    res_name = "TIP4P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_tip4p_ew()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip4p_Ew moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function structure_tip4p_ew()
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
    para_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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

    res_name = "TIP4P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_tip4p_fq()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip4p_FQ moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function structure_tip4p_fq()
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
    para_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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

    res_name = "TIP4P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_tip4p_2005()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip4p_2005 moedel used in lammps with lj/cut/tip4p/cut. Detail info: https://lammps.sandia.gov/doc/Howto_tip4p.html

# Model info
- O mass = 15.9994
- H mass = 1.008
- O charge = -1.1128
- H charge = 0.5564
- ``r_0`` of OH bond = 0.9572
- ``θ`` of HOH angle = 104.52
- OM distance = 0.1546
- LJ ``ε`` of O-O = 0.1852
- LJ `` σ`` of O-O = 3.1589
- LJ ``ε``, `` σ`` of OH, HH = 0.0
- Coulomb cutoff = 8.5
"""
function structure_tip4p_2005()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.15            # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    charge_vec = [0.5564, -1.1128]

    atom_type = [1, 2, 2]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    para_mass = [15.9994, 1.008]
    atom_name = split("O H")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        1/ratio[1] 1/ratio[2] 0
        0          0          0
        2/ratio[1] 0          0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 2)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [1, 3])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, [2, 1, 3])
    num_angles = length(angle_mode)
    num_angle_types = 1

    res_name = "TIP4P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end


"""
    structure_tip4p_ice()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip4p_Ice moedel used in lammps with lj/cut/tip4p/cut. Detail info: https://lammps.sandia.gov/doc/Howto_tip4p.html

# Model info
- O mass = 15.9994
- H mass = 1.008
- O charge = -1.1794
- H charge = 0.5897
- ``r_0`` of OH bond = 0.9572
- ``θ`` of HOH angle = 104.52
- OM distance = 0.1577
- LJ ``ε`` of O-O = 0.21084
- LJ `` σ`` of O-O = 3.1668
- LJ ``ε``, `` σ`` of OH, HH = 0.0
- Coulomb cutoff = 8.5"""
function structure_tip4p_ice()
    # Parameters of water
    density = 1 / Const_cm2an^3;      # Unit: g/A^3
    angle = 104.52;
    angle = angle / 180 * pi;   # Convert int to rad
    ratio = [4, 5]              # Ratio to determing shape of box
    bond_len = 0.9572           # Unit: Anstrom. https://en.wikipedia.org/wiki/Water_model for more details
    dummy_len = 0.15            # Unit: Anstrom. http://www1.lsbu.ac.uk/water/water_models.html for more details
    charge_vec = [0.5897, -1.1794]

    atom_type = [1, 2, 2]
    atom_charge = [charge_vec[atom_type[n]] for n=1:length(atom_type)]
    para_mass = [15.9994, 1.008]
    atom_name = split("O H")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    atom_vec = [
        1/ratio[1] 1/ratio[2] 0
        0          0          0
        2/ratio[1] 0          0
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 2)
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [1, 3])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, [2, 1, 3])
    num_angles = length(angle_mode)
    num_angle_types = 1

    res_name = "TIP4P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_tip5p()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip5p moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function structure_tip5p()
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
    para_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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

    res_name = "TIP5P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_tip5p_2018()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip5p_2018 moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function structure_tip5p_2018()
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
    para_mass = [1.00784, 15.9994, 0]
    atom_name = split("H O Dummy")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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

    res_name = "TIP5P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end

"""
    structure_tip7p()

This will generate a `Structure_Wat` type which contains all information needed to build a Tip7p moedel. Detail info: http://www1.lsbu.ac.uk/water/water_models.html
"""
function structure_tip7p()
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
    para_mass = [1.00784, 15.9994, 0, 0]
    atom_name = split("H O Dummy01 Dummy02")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)

    # Setting Atom Args
    cell_vec = [
        ratio[1] * bond_len * sin(angle/2)
        ratio[2] * bond_len * cos(angle/2)
        ( (para_mass[1]*2+para_mass[2]) * Const_gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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

    res_name = "TIP7P"
    Structure_Wat(atom_vec, cell_vec, atom_type, atom_name, atom_charge, para_mass, num_atoms, 
    num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types, [1, 2], res_name)
end
