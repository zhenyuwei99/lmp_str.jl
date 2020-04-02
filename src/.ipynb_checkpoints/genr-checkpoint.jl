# genr_cell
"""
    genr_cell(cell_vec)

Do this will return a variable in `Data_Cell` type

# Example
```julia-repl
data_cell = genr_cell([1 2 3])
```
"""
function genr_cell(cell_vec)
    num_cells = 1
    for i in cell_vec
        num_cells *= i
    end

    cell_mat = Matrix{Int64}(undef, (num_cells, 3))
    for x in 0 : cell_vec[1] - 1
        for y in 0 : cell_vec[2] - 1
            for z in 0 : cell_vec[3] - 1
                id = x*cell_vec[2]*cell_vec[3] + y*cell_vec[3] + z + 1
                cell_mat[id, :] = [x y z]
            end
        end
    end
    if typeof(cell_vec) == Array{Int64,2}
        cell_vec = conv(cell_vec, 1)
    end
    Data_Cell(cell_mat, cell_vec, num_cells)
end

# genr
"""
    genr(data_cell::Data_Cell, str::Str)

Do this will return a variable in `Data_Unit` type, containing all infomation needed for model of `str`.

# Example
```julia-repl
data_cell = genr_cell([5, 5, 5])
str = Tip3p()
data = genr(data_cell, str)
```
"""
function genr(data_cell::Data_Cell, str::Str)
    data = genr_atom(data_cell, str)
    str_fields = fieldnames(typeof(str))

    vec_bond = 0
    vec_angle = 0

    if in(:bond_mode, str_fields)
        data_bond = genr_bond(data_cell, data)
        vec_bond = data_bond.vec_bond
    end

    if in(:angle_mode, str_fields)
        data_angle = genr_angle(data_cell, data)
        vec_angle = data_angle.vec_angle
    end

    Data_Unit(data.data_basic, data.data_str, data.vec_atom, vec_bond, vec_angle)
end

# genr_atom
"""
    genr_atom(data_cell::Data_Cell, str::Str)

This will be called by `genr(...)`, genrating atomic information
"""
function genr_atom(data_cell::Data_Cell, str::Str)
    # Variables Setting
    vec_atom = Vector{Atom}(undef, data_cell.num_cells*str.num_atoms)

    # Generate Atom Data
    for cell in 1 : data_cell.num_cells
        for atom in 1 : str.num_atoms
            id_now = (cell-1) * str.num_atoms + atom
            coord = conv((data_cell.cell_mat[cell, :]+str.atom_vec[atom, :]), 1) * str.cell_vec
            coord = conv(coord, 1) # Reduce Matrix to Vector
            vec_atom[id_now] = Atom(id_now, cell, str.atom_type[atom], str.atom_charge[atom], coord)
        end
    end

    # Generate Basic Data
    num_atoms = data_cell.num_cells * str.num_atoms
    num_bonds = 0
    num_angles = 0
    num_dihedrals = 0
    num_impropers = 0
    num_atom_types = str.num_atom_types
    num_bond_types = 0
    num_angle_types = 0
    num_dihedral_types = 0
    num_improper_types = 0
    box_size = Matrix{Float64}(undef, (3, 2))
    str_vec = Matrix{Float64}(undef, (3, 3))
    str_vec[:, :] = str.cell_vec[:, :]
    str_vec[2,1] = 0
    str_vec[3,1] = 0
    str_vec[3,2] = 0
    box_size[:, 1] = min(vec_atom, "coord")
    box_size[:, 2] = conv(conv(min(vec_atom, "coord"), 1) + conv(data_cell.cell_vec, 1)*str_vec, 1)
    box_tilt = Vector{Float64}(undef, 3)
    box_tilt[1] = data_cell.cell_vec[2] * str.cell_vec[2,1]
    box_tilt[2] = data_cell.cell_vec[3] * str.cell_vec[3,1]
    box_tilt[3] = data_cell.cell_vec[2] * str.cell_vec[3,2]

    data_basic = Data_Basic(num_atoms, num_bonds, num_angles, num_dihedrals, num_impropers, num_atom_types, num_bond_types, num_angle_types, num_dihedral_types, num_improper_types, box_size, box_tilt)

    # Output
    Data_Unit(data_basic, str, vec_atom, 0, 0)
end

# genr_bond
"""
    genr_bond(data_cell::Data_Cell, data::Data)

This will be called by `genr(...)` if `bond_mode` is contained in `str`, genrating information of bonds
"""
function genr_bond(data_cell::Data_Cell, data::Data)
    # Reading Input
    data_basic = data.data_basic
    data_str = data.data_str
    num_cells = data_cell.num_cells
    num_cell_bonds = data_str.num_bonds
    num_bonds = num_cells * num_cell_bonds
    num_bond_types = data_str.num_bond_types

    data_basic.num_bonds = num_bonds
    data_basic.num_bond_types = num_bond_types

    # Generate Bonds
    vec_bond = Vector{Bond}(undef, num_bonds)
    for cell in 1 : num_cells
        for bond = 1 : num_cell_bonds
            id_now = (cell-1) * num_cell_bonds + bond
            bond_now = Bond(data_str.bond_mode[bond])
            bond_now.id = id_now
            atom_tilt = (cell-1) * data_str.num_atoms
            bond_now.atom .+= atom_tilt
            vec_bond[id_now] = bond_now
        end
    end

    # Output
    Data_Unit(data_basic, data_str, data.vec_atom, vec_bond, 0)
end

# genr_angle
"""
    genr_bond(data_cell::Data_Cell, data::Data)

This will be called by `genr(...)` if `angle_mode` is contained in `str`,
genrating information of angles
"""
function genr_angle(data_cell::Data_Cell, data::Data)
    # Reading Input
    data_basic = data.data_basic
    data_str = data.data_str
    num_cells = data_cell.num_cells
    num_cell_angles = data_str.num_angles
    num_angles = num_cells * num_cell_angles
    num_angle_types = data_str.num_angle_types

    data_basic.num_angles = num_angles
    data_basic.num_angle_types = num_angle_types

    # Generate Angle
    vec_angle = Vector{Angle}(undef, num_angles)
    for cell in 1 : num_cells
        for angle in 1 : num_cell_angles
            id_now = (cell-1) * num_cell_angles + angle
            angle_now = Angle(data_str.angle_mode[angle])
            angle_now.id = id_now
            atom_tilt = (cell-1) * data_str.num_atoms
            angle_now.atom .+= atom_tilt
            vec_angle[id_now] = angle_now
        end
    end

    # Output
    Data_Unit(data_basic, data_str, data.vec_atom, data.vec_bond, vec_angle)
end