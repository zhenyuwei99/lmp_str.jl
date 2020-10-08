"""
    function read_pdb(file_pdb)

This will read a .pdb file, getting all the coordinate information.
- `file_pdb`: path of target .pdb file.

# Output
- `id_vec`
- ``

# Example
```julia-repl
file_pdb = "gb1.pdb"
coord = read_pdb(file_pdb)
```
"""
function read_pdb(file_pdb)

    # Setting variabls
    col_range_id = 7:11
    col_range_atom_name = 13:16
    col_range_mol = 23:26
    col_range_x = 31:38
    col_range_y = 39:46
    col_range_z = 47:54

    # Reading file
    io = open(file_pdb)
    str = read(io, String)
    reg = r"ATOM..*[A-Z]"
    res = findall(reg, str)
    res = [str[n] for n in res]

    # Ansys
    num_atoms = length(res)
    atom_vec = [String(strip(res[atom][col_range_atom_name], ' ')) for atom=1:num_atoms]
    mol_vec = zeros(Int, num_atoms, 1)
    id_vec = zeros(Int, num_atoms, 1)
    coord_mat = zeros(num_atoms, 3)

    for atom = 1:num_atoms
        id_vec[atom, 1] = parse(Int, res[atom][col_range_id])
        mol_vec[atom, 1] = parse(Int, res[atom][col_range_mol])
        coord_mat[atom, 1] = parse(Float64, res[atom][col_range_x])
        coord_mat[atom, 2] = parse(Float64, res[atom][col_range_y])
        coord_mat[atom, 3] = parse(Float64, res[atom][col_range_z])
    end
    
    return id_vec, atom_vec, mol_vec, coord_mat
end

"""
    function read_pdb()

This will read a .pdb file, getting all the coordinate information.
- `file_pdb`: path of target .pdb file.

# Example
```julia-repl
file_pdb = "gb1.pdb"
coord = read_pdb(file_pdb)
```
"""
function read_psf(file_psf)
    # Read psf file
    io = open(file_psf)
    raw_info = split(read(io, String), "\n")
    raw_start = findall(x->occursin("!N", x), raw_info)[2:7] .+ 1 # Ignore !NTITLE
    raw_end = raw_start .- 3 

    # Atom info
    col_type = 6
    col_charge = 7
    num_atoms = raw_end[2] - raw_start[1] + 1
    type_vec = [String(split(raw_info[i])[col_type]) for i=raw_start[1]:raw_end[2]]
    charge_vec = [parse(Float64, String(split(raw_info[i])[col_charge])) for i=raw_start[1]:raw_end[2]]

    # Bond info
    bond_topo = [strvec2intvec(split(raw_info[i])) for i=raw_start[2]:raw_end[3]]
    bond_topo = Array(reshape(vcat(bond_topo...), (2, :))')

    # Angle info
    angle_topo = [strvec2intvec(split(raw_info[i])) for i=raw_start[3]:raw_end[4]]
    angle_topo = Array(reshape(vcat(angle_topo...), (3, :))')

    # Dihedral info
    dihedral_topo = [strvec2intvec(split(raw_info[i])) for i=raw_start[4]:raw_end[5]]
    dihedral_topo = Array(reshape(vcat(dihedral_topo...), (4, :))')

    # Improper info
    improper_topo = [strvec2intvec(split(raw_info[i])) for i=raw_start[5]:raw_end[6]]
    improper_topo = Array(reshape(vcat(improper_topo...), (4, :))')
    
    return type_vec, charge_vec, bond_topo, angle_topo, dihedral_topo, improper_topo
end

function assign_para_atom(id_vec, mol_vec, type_vec, charge_vec, coord_mat, potential::T) where T <:Potential
    num_atoms = size(id_vec, 1)
    type_list = unique(type_vec)
    mol_vec = mol_vec .- minimum(mol_vec) .+ 1  # Mol id start from 1
    num_atom_types = size(type_list, 1)

    vec_atom = Vector{Atom}(undef, num_atoms)
    
    for atom = 1 : num_atoms
        vec_atom[atom] = Atom(id_vec[atom], mol_vec[atom], findfirst(x->x==type_vec[atom], type_list),
                        charge_vec[atom], coord_mat[atom, :])
    end

    para_pair = [findfirst(x->judgeequal(x, [type_list[i] for i in id_vec[atom]]),
                             potential.para_pair) for atom = 1:num_atom_types]
    para_pair = vcat([pair.para for pair in potential.para_pair[para_pair]]...)
    para_pair = reshape(para_pair, (4, :))'

    para_mass = [findfirst(x->judgeequal(x, [type_list[i] for i in id_vec[atom]]),
    potential.para_mass) for atom = 1:num_atom_types]
    para_mass = vcat([pair.para for pair in potential.para_mass[para_mass]]...)


    return num_atoms, num_atom_types, vec_atom, para_pair, para_mass, type_list
end

function assign_para_bond(type_vec, bond_topo, potential)
    num_bonds = size(bond_topo, 1)
    num_bond_types = 0
    vec_bond = Vector{Bond}(undef, num_bonds)
    type_list = [0]
    for id = 1:num_bonds
        para_type_id = findfirst(x->judgeequal(x, [type_vec[i] for i in bond_topo[id, :]]), potential.para_bond)
        if para_type_id == nothing
            print([type_vec[i] for i in bond_topo[id, :]])
            error("No matching bond type is found in Charmm36 force filed.")
        end
        if !in(para_type_id, type_list)
            num_bond_types += 1
            append!(type_list, para_type_id)
        end
        vec_bond[id] = Bond(id, findfirst(x->x==para_type_id, type_list) - 1, bond_topo[id, :])
    end
    para_bond = [bond.para for bond in potential.para_bond[type_list[2:end]]]
    para_bond = Array(hcat(para_bond...)')
    return num_bonds, num_bond_types, vec_bond, para_bond
end

function assign_para_angle(type_vec, angle_topo, potential)
    num_angles = size(angle_topo, 1)
    num_angle_types = 0
    vec_angle = Vector{Angle}(undef, num_angles)
    type_list = [0]
    for id = 1:num_angles
        para_type_id = findfirst(x->judgeequal(x, [type_vec[i] for i in angle_topo[id, :]]), potential.para_angle)
        if para_type_id == nothing
            print([type_vec[i] for i in angle_topo[id, :]])
            error("No matching angle type is found in Charmm36 force filed.")
        end
        if !in(para_type_id, type_list)
            num_angle_types += 1
            append!(type_list, para_type_id)
        end
        vec_angle[id] = Angle(id, findfirst(x->x==para_type_id, type_list)-1, angle_topo[id, :])
    end
    para_angle = [angle.para for angle in potential.para_angle[type_list[2:end]]]
    para_angle = Array(hcat(para_angle...)')
    return num_angles, num_angle_types, vec_angle, para_angle
end

function assign_para_dihedral(type_vec, dihedral_topo, potential)
    num_dihedrals = size(dihedral_topo, 1)
    num_dihedral_types = 0
    vec_dihedral = Vector{Dihedral}(undef, num_dihedrals)
    type_list = [0]
    for id = 1:num_dihedrals
        para_type_id = findfirst(x->judgeequal(x, [type_vec[i] for i in dihedral_topo[id, :]]), potential.para_dihedral)
        if para_type_id == nothing
            print([type_vec[i] for i in dihedral_topo[id, :]])
            error("No matching dihedral type is found in Charmm36 force filed.")
        end
        if !in(para_type_id, type_list)
            num_dihedral_types += 1
            append!(type_list, para_type_id)
        end
        vec_dihedral[id] = Dihedral(id, findfirst(x->x==para_type_id, type_list)-1, dihedral_topo[id, :])
    end
    para_dihedral = [dihedral.para for dihedral in potential.para_dihedral[type_list[2:end]]]
    para_dihedral = Array(hcat(para_dihedral...)')
    return num_dihedrals, num_dihedral_types, vec_dihedral, para_dihedral
end

function assign_para_improper(type_vec, improper_topo, potential)
    num_impropers = size(improper_topo, 1)
    num_improper_types = 0
    vec_improper = Vector{Improper}(undef, num_impropers)
    type_list = [0]
    for id = 1:num_impropers
        para_type_id = findfirst(x->judgeequal(x, [type_vec[i] for i in improper_topo[id, :]]), potential.para_improper)
        if para_type_id == nothing
            print([type_vec[i] for i in improper_topo[id, :]])
            error("No matching improper type shwon above is found in Charmm36 force filed.")
        end
        if !in(para_type_id, type_list)
            num_improper_types += 1
            append!(type_list, para_type_id)
        end
        vec_improper[id] = Improper(id, findfirst(x->x==para_type_id, type_list)-1, improper_topo[id, :])
    end
    para_improper = [improper.para for improper in potential.para_improper[type_list[2:end]]]
    para_improper = Array(hcat(para_improper...)')
    return num_impropers, num_improper_types, vec_improper, para_improper
end


"""
    function convert_vmd(file_pdb::String, file_psf::String, potential::Potential)
This will convert a set of files which created by VMD (.pdb and .psf) file into a `Data` variable that can be used to write a .data file.
- `file_pdb`: String of path of .pdb file.
- `file_psf`: String of path of .psf file.
- `poetntial`: Subclass of `Potential` which contains all information of forcefiled to automatically generate parameters, e.g. pair_coeff, bond_coeff etc.
"""
function convert_vmd(file_pdb::String, file_psf::String, potential::Potential, box_padding=2)
    # Reading Input
    id_vec, atom_vec, mol_vec, coord_mat = read_pdb(file_pdb) 
    type_vec, charge_vec, bond_topo, angle_topo, dihedral_topo, improper_topo = read_psf(file_psf)
    
    box_size = Array(cat(minimum(coord_mat, dims=1), maximum(coord_mat, dims=1), dims=1)')
    box_size[:, 1] .-= box_padding / 2
    box_size[:, 2] .+= box_padding / 2
    box_tilt = [0, 0, 0] 

    # Assign parameters
    num_atoms, num_atom_types, vec_atom, para_pair, para_mass, type_list = assign_para_atom(id_vec, 
            mol_vec, type_vec, charge_vec, coord_mat, potential_charmm36)
    num_bonds, num_bond_types, vec_bond, para_bond = assign_para_bond(type_vec, bond_topo, potential)
    num_angles, num_angle_types, vec_angle, para_angle = assign_para_angle(type_vec, angle_topo, potential)
    num_dihedrals, num_dihedral_types, vec_dihedral, para_dihedral = assign_para_dihedral(type_vec, dihedral_topo, potential)
    num_impropers, num_improper_types, vec_improper, para_improper = assign_para_improper(type_vec, improper_topo, potential)

    # Creat Data instance
    data_basic = Data_Basic(num_atoms, num_bonds, num_angles, num_dihedrals, num_impropers, 
    num_atom_types, num_bond_types, num_angle_types, num_dihedral_types, num_improper_types, box_size, box_tilt)
    data_str = Structure_VMD(type_list, num_atom_types, num_bond_types, num_angle_types, num_dihedral_types, num_improper_types,
        para_mass, para_pair, para_bond, para_angle, para_dihedral, para_improper)
    data = Data_Unit(data_basic, data_str,
        vec_atom, vec_bond, vec_angle, vec_dihedral, vec_improper)
    return data
end

