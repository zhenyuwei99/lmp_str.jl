module lmp_str

import Dates

## Global Variables

# Physical Constants
k_b = 1.38065e-23;
n_a = 6.02214e23;
density_wat = 1e3;      # Unit: kg/m^3
# Energy Converters
kcal2j = 4.184e3;
kcalm2j = kcal2j/n_a;
kcalm2t = kcalm2j/k_b;
# Mass Converters
g2kg = 1e-3;
kg2g = 1e3;
gm2g = 1/n_a;
gm2kg = g2kg/n_a;
# Length Converters
an2m = 1e-10;
an2nm = 1e-1;
nm2m = 1e-9;
cm2an = 1e8;
dm2m = 1e-1;
m2dm = 1/dm2m;
# Time Converters
fs2s = 1e-15;
ps2s = 1e-12;
ns2s = 1e-9;

## Definations of Types

# Type of Structure
abstract type Str end

mutable struct Tip3p <: Str
    atom_vec
    cell_vec
    atom_type
    atom_charge
    atom_mass
    atom_name
    num_atoms
    num_atom_types
    bond_mode
    num_bonds
    num_bond_types
    angle_mode
    num_angles
    num_angle_types
end

function Tip3p()
    # Parameters of water
    density = 1 / cm2an^3;      # Unit: g/A^3
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
        ( (atom_mass[1]*2+atom_mass[2]) * gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
    ]
    cell_vec = diag(cell_vec)

    # Setting Bond Modes
    bond_mode = Vector{Bond}(undef, 2)
    bond_mode[1] = Bond(0, 1, 1, 2)
    bond_mode[2] = Bond(0, 1, 2, 3)
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, 1, 2, 3)
    num_angles = length(angle_mode)
    num_angle_types = 1
    Tip3p(atom_vec, cell_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types)
end

mutable struct Si3N4 <: Str
    cell_vec
    atom_vec
    atom_type
    atom_charge
    atom_mass
    atom_name
    num_atoms
    num_atom_types
end

function Si3N4()
   cell_vec = [
        7.595  0         0
        3.7975 6.577463  0
        0      0         2.902
        ]
    atom_vec = [
        -0.0977580565637542 0.195516113127508   0
        0.0588368288671674  0.307869462739661   0
        -0.0535176961032620 0.576818144016926   0
        0.327654573699065   0.151274131074550   0
        0.403949967083608   0                   0.500000000000000
        0.438645801493259   0.309541840068124   0
        0.749567567747191   0.272293435934189   0.500000000000000
        0.0180020380857845  0.645233580181295   0.500000000000000
        0.129048914298138   0.803653323477456   0.500000000000000
        0.510352849669919   0.378109310535080   0.500000000000000
        0.397866659130036   0.647057991812345   0.500000000000000
        0.554461544560957   0.759411341424498   0.500000000000000
        -0.291760564557441  0.683850293038517   0
        0.0514815980747357  0.957207969090818   0
    ]
    atom_type = [2 1 2 2 1 1 1 1 2 2 1 2 1 1]
    atom_charge = 1.34925 .* atom_type .- 1.9275 # N: -0.57825. Si: 0.7710. Unit: e
    atom_mass = [14.0067 28.085501] # N: 14.0067. Si: 28.085501. Unit: g/mol
    atom_name = split("N Si")
    num_atoms = length(atom_type)
    num_atom_types = 2
    Si3N4(cell_vec, atom_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types)
end

mutable struct Si3N4_Ort <: Str
    cell_vec
    atom_vec
    atom_type
    atom_charge
    atom_mass
    atom_name
    num_atoms
    num_atom_types
end

function Si3N4_Ort()
    cell_vec = [
    7.595   0           0
    0       13.154964   0
    0       0           2.902
    ]
    atom_vec = [
    0                       0.597031356376194   0
    0.000658327847268042    0.521394509327430   0.500000000000000
    0.0308755760368664      0.401100527527099   0.500000000000000
    0.0967083607636604      0.119153727824721   0
    0.126793943383805       0                   0
    0.127583936800527       0.923220010332222   0.500000000000000
    0.146872942725477       0.363320112468571   0
    0.190125082290981       0.676164982283494   0
    0.296115865701119       0.710448618483487   0.500000000000000
    0.309479921000658       0.175330240356416   0
    0.318104015799868       0.844922570673702   0.500000000000000
    0.331599736668861       0.309804192546631   0
    0.437327188940092       0.344011811814916   0.500000000000000
    0.482422646477946       0.657540834015205   0.500000000000000
    0.500000000000000       0.0970328006978962  0
    0.500658327847268       0.0213959536491320  0.500000000000000
    0.530875576036866       0.901099083205397   0.500000000000000
    0.596708360763660       0.619152283503018   0
    0.626793943383805       0.499998555678298   0
    0.627583936800527       0.423221454653924   0.500000000000000
    0.646872942725477       0.863318668146869   0
    0.690125082290981       0.176166426605196   0
    0.796115865701119       0.210450062805189   0.500000000000000
    0.809479921000659       0.675328796034714   0
    0.818104015799869       0.344924014995404   0.500000000000000
    0.831599736668861       0.809802748224929   0
    0.937327188940092       0.844010367493214   0.500000000000000
    0.982422646477946       0.157542278336908   0.500000000000000
    ]
    atom_type = [2 1 2 2 1 2 1 1 2 1 1 2 1 1 2 1 2 2 1 2 1 1 2 1 1 2 1 1]
    atom_charge = 1.34925 .* atom_type .- 1.9275 # N: -0.57825. Si: 0.7710. Unit: e
    atom_mass = [14.0067 28.085501] # N: 14.0067. Si: 28.085501. Unit: g/mol
    atom_name = split("N Si")
    num_atoms = length(atom_type)
    num_atom_types = 2
    Si3N4_Ort(cell_vec, atom_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types)
end

# Type of Data

abstract type Data end
abstract type Unit end

mutable struct Atom <: Unit
    id::Int64
    mol::Int64
    atom_type::Int64
    charge::Float64
    coord::AbstractArray
end

mutable struct Bond <: Unit
    id::Int64
    bond_type::Int64
    atom_01::Int64
    atom_02::Int64
end

function Bond(bond::Bond)
    id = bond.id
    bond_type = bond.bond_type
    atom_01 = bond.atom_01
    atom_02 = bond.atom_02
    Bond(id, bond_type, atom_01, atom_02)
end

mutable struct Angle <: Unit
    id::Int64
    angle_type::Int64
    atom_01::Int64
    atom_02::Int64
    atom_03::Int64
end

function Angle(angle::Angle)
    id = angle.id
    angle_type = angle.angle_type
    atom_01 = angle.atom_01
    atom_02 = angle.atom_02
    atom_03 = angle.atom_03
    Angle(id, angle_type, atom_01, atom_02, atom_03)
end

mutable struct Data_Cell
    cell_mat::Matrix{Int64}
    cell_vec
    num_cells::Int64
end

mutable struct Data_Basic <:Data
    num_atoms::Int64
    num_bonds::Int64
    num_angles::Int64
    num_dihedrals::Int64
    num_impropers::Int64
    num_atom_types::Int64
    num_bond_types::Int64
    num_angle_types::Int64
    num_dihedral_types::Int64
    num_improper_types::Int64
    box_size
    box_tilt
end

mutable struct Data_Atom <: Data
    data_basic::Data_Basic
    data_str::Str
    vec_atom::Vector{Atom}
end

mutable struct Data_Bond <: Data
    data_basic::Data_Basic
    data_str::Str
    vec_atom::Vector{Atom}
    vec_bond::Vector{Bond}
end

mutable struct Data_Angle <: Data
    data_basic::Data_Basic
    data_str::Str
    vec_atom::Vector{Atom}
    vec_bond::Vector{Bond}
    vec_angle::Vector{Angle}
end


## Definations of Functions

# Useful Functions

function max(vec::Vector{Atom})
    max = vec[1].coord
    for i = 2 : length(vec)
        for j = 1 : 3
            if max[j] < vec[i].coord[j]
                max[j] = vec[i].coord[j]
            end
        end
    end
    result = Vector{Float64}(undef, 3)
    result[:] = max[:]
end

function min(vec::Vector{Atom})
    min = vec[1].coord
    for i = 2 : length(vec)
        for j = 1 : 3
            if min[j] > vec[i].coord[j]
                min[j] = vec[i].coord[j]
            end
        end
    end
    result = Vector{Float64}(undef, 3)
    result[:] = min[:]
end

function diag(vec::AbstractArray)
    len = length(vec)
    result = zeros(len, len)
    for i = 1 : len
        result[i, i] = vec[i]
    end
    result
end

function conv(vec::Array{Int64, 1})
    len = length(vec)
    result = convert.(Int64, zeros(1, len))
    for i in 1 : len
        result[i] = vec[i]
    end
    result
end

function conv(vec::Array{Int64, 2})
    len = length(vec)
    result = convert.(Int64, zeros(len, 1))
    for i in 1 : len
        result[i] = vec[i]
    end
    result
end

function dist(atom::Atom, org::Array)
    r = 0
    for dim = 1 : 3
        r += (org[dim] - atom.coord[dim])^2
    end
    sqrt(r)
end

function dist(atom::Atom, org::Array, dim)
    r = 0
    dim_range = [1 2 3]
    dim_range = dim_range[1:end .!= dim]
    for i in dim_range
        r += (org[i] - atom.coord[i])^2
    end
    sqrt(r)
end

function dist(pos::Array, org::Array)
    r = 0
    for dim = 1 : 3
    r += (org[dim] - pos[dim])^2
    end
    sqrt(r)
end

function dist(pos::Array, org::Array, dim)
    r = 0
    dim_range = [1 2 3]
    dim_range = dim_range[1:end .!= dim]
    for i in dim_range
        r += (org[i] - pos[i])^2
    end
    sqrt(r)
end

function delete(mat; id=1, dim=1)
    if dim == 1
        result = mat[1:end .!= id, :]
    elseif dim == 2
        result = mat[:, 1:end .!= id]
    else
        error("Error, dim should be 1 or 2, representing row or column respectively!")
    end
    result
end

# genr_cell

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
    if typeof(cell_vec) == Array{Int64,1}
        cell_vec = conv(cell_vec)
    end
    Data_Cell(cell_mat, cell_vec, num_cells)
end

# genr_atom

function genr_atom(data_cell::Data_Cell, str::Str)
    # Variables Setting
    vec_atom = Vector{Atom}(undef, data_cell.num_cells*str.num_atoms)

    # Generate Atom Data
    for cell in 1 : data_cell.num_cells
        for atom in 1 : str.num_atoms
            id_now = (cell-1) * str.num_atoms + atom
            coord = (data_cell.cell_mat[cell,:]+str.atom_vec[atom,:])' * str.cell_vec
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
    box_size[:, 1] = min(vec_atom)
    box_size[:, 2] = min(vec_atom)' + data_cell.cell_vec*str_vec
    box_tilt = Vector{Float64}(undef, 3)
    box_tilt[1] = data_cell.cell_vec[2] * str.cell_vec[2,1]
    box_tilt[2] = data_cell.cell_vec[3] * str.cell_vec[3,1]
    box_tilt[3] = data_cell.cell_vec[2] * str.cell_vec[3,2]

    data_basic = Data_Basic(num_atoms, num_bonds, num_angles, num_dihedrals, num_impropers, num_atom_types, num_bond_types, num_angle_types, num_dihedral_types, num_improper_types, box_size, box_tilt)

    # Output
    Data_Atom(data_basic, str, vec_atom)
end

# genr_bond

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
            bond_now.atom_01 += atom_tilt
            bond_now.atom_02 += atom_tilt
            vec_bond[id_now] = bond_now
        end
    end

    # Output
    Data_Bond(data_basic, data_str, data.vec_atom, vec_bond)
end

# genr_angle

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
            angle_now.atom_01 += atom_tilt
            angle_now.atom_02 += atom_tilt
            angle_now.atom_03 += atom_tilt
            vec_angle[id_now] = angle_now
        end
    end

    # Output
    Data_Angle(data_basic, data_str, data.vec_atom, data.vec_bond, vec_angle)
end

# move

function move(data::Data, move_vec)
    # Motion of box
    if typeof(move_vec) == Array{Int64,2}
        move_vec_box = conv(move_vec)
        move_vec_atom = move_vec
    else
        move_vec_atom = conv(move_vec)
        move_vec_box = move_vec
    end
    data.data_basic.box_size = broadcast(+, data.data_basic.box_size, move_vec_box)

    # Motion of atoms
    for atom = 1 : data.data_basic.num_atoms
        data.vec_atom[atom].coord += move_vec_atom
    end
    data
end

# select and delete

function select(data_cell::Data_Cell; mode::String, para, dim=3)
    list_mode = split("cylinder sphere")
    if !in(mode, list_mode)
        error(join(["Error, mode: ",mode," is not supported."]))
    end

    center = data_cell.cell_vec ./ 2

    cell_list = 0

    if mode == "cylinder"
        for cell in 1:data_cell.num_cells
            if dist(data_cell.cell_mat[cell, :], center, dim) <= para
                cell_list = vcat(cell_list, cell)
            end
        end
    end

    if mode == "sphere"
        for cell in 1:data_cell.num_cells
            if dist(data_cell.cell_mat[cell, :], center) <= para
                cell_list = vcat(cell_list, cell)
            end
        end
    end

    cell_list[2:end]
end

function select(data_atom::Data; mode::String, para, dim=3)
    list_mode = split("cylinder sphere")
    if !in(mode, list_mode)
        error(join(["Error, mode: ",mode," is not supported."]))
    end

    box_size = data_atom.data_basic.box_size
    box_tilt = data_atom.data_basic.box_tilt
    center = (max(data_atom.vec_atom)-min(data_atom.vec_atom)) ./ 2
    #center[1] += box_tilt[1] + box_tilt[2]
    #center[2] += box_tilt[3]

    atom_list = 0

    if mode == "cylinder"
        for atom in 1:data_atom.data_basic.num_atoms
            if dist(data_atom.vec_atom[atom], center, dim) <= para
                atom_list = vcat(atom_list, atom)
            end
        end
    end

    if mode == "sphere"
        for atom in 1:data_atom.data_basic.num_atoms
            if dist(data_atom.vec_atom[atom], center) <= para
                atom_list = vcat(atom_list, atom)
            end
        end
    end
    atom_list[2:end]
end

function sort(vec::Union{Vector{Atom}, Vector{Bond}, Vector{Angle}})
    for unit in 1 : length(vec)
        vec[unit].id = unit
    end
    vec
end

function delete(vec::Union{Vector{Atom}, Vector{Bond}, Vector{Angle}}, id::Array)
    len = length(vec)
    judge = trues(len)
    for i in id
        judge .&= 1 : len .!= i
    end
    vec[judge]
end

function delete(mat; id=1, dim=1)
    if dim == 1
        result = mat[1:end .!= id, :]
    elseif dim == 2
        result = mat[:, 1:end .!= id]
    else
        error("Error, dim should be 1 or 2, representing row or column respectively!")
    end
    result
end

function delete(data_cell::Data_Cell, list_cell)
    # Reading Input
    cell_mat = data_cell.cell_mat
    len = size(cell_mat)[1]
    # Deleting Cells
    judge = trues(len)
    for id in list_cell
        judge .&= 1:len .!= id
    end
    data_cell.cell_mat = cell_mat[judge, :]
    data_cell.num_cells = size(data_cell.cell_mat)[1]

    # Output
    data_cell
end

function delete(data::Data, list_atom::Array)
    fields = fieldnames(typeof(data))
    name_fields = [string(fields[n]) for n = 1:length(fields)]
    num_fields = length(fields)
    list_fields = findall(x->occursin("vec", x), name_fields)

    for field in  list_fields
        # Find all elements that need to be deleted
        list_id = find(getfield(data, fields[field]), list_atom)
        # Changing # of specifc field
        num_ids = length(list_id) - 1 # First elements is "Array"
        para_now = name_fields[field][findall(x->in('_', x), name_fields[field])[1]+1 : end]
        para_now = Meta.parse(join(["num_",para_now,"s"]))
        para_result = getfield(data.data_basic, para_now) - num_ids
        setfield!(data.data_basic, para_now, para_result)
        # Changing Vector of each field
        setfield!(data, fields[field], delete(getfield(data, fields[field]), list_id))
    end

    data
end

function find(vec_atom::Vector{Atom}, list_atom)
    len = length(vec_atom)
    judge = falses(len)

    id = [vec_atom[n].id for n = 1:len]
    for atom in list_atom
        judge .|= id .== atom
    end

    list = Vector{Int64}
    for vec in vec_atom[judge]
        list = vcat(list, vec.id)
    end
    list
end

function find(vec_bond::Vector{Bond}, list_atom)
    len = length(vec_bond)
    judge_01 = falses(len)
    judge_02 = falses(len)

    atom_01 = [vec_bond[n].atom_01 for n = 1:len]
    atom_02 = [vec_bond[n].atom_02 for n = 1:len]

    for atom in list_atom
        judge_01 .|= atom_01 .== atom
        judge_02 .|= atom_02 .== atom
    end

    judge = judge_01 .| judge_02

    list = Vector{Int64}
    for vec in vec_bond[judge]
        list = vcat(list, vec.id)
    end
    list
end

function find(vec_angle::Vector{Angle}, list_atom)
    len = length(vec_angle)
    judge_01 = falses(len)
    judge_02 = falses(len)
    judge_03 = falses(len)

    atom_01 = [vec_angle[n].atom_01 for n = 1:len]
    atom_02 = [vec_angle[n].atom_02 for n = 1:len]
    atom_03 = [vec_angle[n].atom_02 for n = 1:len]

    for atom in list_atom
        judge_01 .|= atom_01 .== atom
        judge_02 .|= atom_02 .== atom
        judge_03 .|= atom_03 .== atom
    end

    judge = judge_01 .| judge_02 .| judge_03

    list = Vector{Int64}
    for vec in vec_angle[judge]
        list = vcat(list, vec.id)
    end
    list
end

# write_data and write_info

function write_data(data::Data, name_file::AbstractString)
    for info in fieldnames(typeof(data))
        write_info(getfield(data, info), name_file)
    end
end

function write_info(info::Data_Basic, name_file::AbstractString)
    # List Setting
    list_01 = "atoms bonds angles dihedrals impropers"
    list_02 = ["atom types", "bond types", "angle types", "dihedral types", "improper types"]
    list = vcat(split(list_01), list_02)
    num_paras = length(list)

    # Writting Para info
    io = open(name_file, "w")
    write(io, join(["Lammps .data file creat at ", Dates.now(), """ by Julia Package "lmp_str"\n\n"""]))

    fields = fieldnames(typeof(info))

    for para in 1 : num_paras
        write(io, join([string(getfield(info, fields[para])), "\t\t", list[para], "\n"]))
    end

    # Writing Box info

    para = num_paras + 1
    box_range = getfield(info, fields[para])
    write(io, join([string(box_range[1,1]), "\t\t", string(box_range[1,2]), "\t\t xlo xhi\n"]))
    write(io, join([string(box_range[2,1]), "\t\t", string(box_range[2,2]), "\t\t ylo yhi\n"]))
    write(io, join([string(box_range[3,1]), "\t\t", string(box_range[3,2]), "\t\t zlo zhi\n"]))

    para += 1
    box_tilt = getfield(info, fields[para])
    for tilt in box_tilt
        write(io, join([string(tilt), "\t"]))
    end
    write(io, "xy xz yz\n")

    close(io)
end

function write_info(info::Str, name_file::AbstractString)
    # Reading Input
    fields = fieldnames(typeof(info))
    io = open(name_file, "a")

    # Judgement and Output
    if in(:atom_mass, fields)
        write(io, "\n\nMasses\n\n")
        for atom = 1 : info.num_atom_types
            write(io, join([string(info.atom_mass[atom]), "\t\t# ", info.atom_name[atom], "\n"]))
        end
    end
    write(io, "\n\n")

    close(io)
end

function write_info(info::Vector{Atom}, name_file::AbstractString)
    # Reading Input
    io = open(name_file, "a")
    fields = fieldnames(typeof(info[1]))
    num_fields = length(fields)

    # Writing Output
    write(io, "Atoms # full\n\n")

    for atom in info
        for para in 1 : num_fields-1
            write(io, join([string(getfield(atom, fields[para])), " "]))
        end
        para = num_fields
        coord = getfield(atom,fields[para])
        for dim = 1 : 3
            write(io, join([string(coord[dim])," "]))
        end
        write(io, "\n")
    end
    write(io, "\n\n")

    close(io)
end

function write_info(info::Vector{Bond}, name_file::AbstractString)
    # Reading Input
    io = open(name_file, "a")
    fields = fieldnames(typeof(info[1]))
    num_fields = length(fields)

    # Writing Output
    write(io, "Bonds\n\n")

    for bond in info
        for para in fields
            write(io, join([string(getfield(bond, para)), " "]))
        end
        write(io, "\n")
    end
    write(io, "\n\n")

    close(io)
end

function write_info(info::Vector{Angle}, name_file::AbstractString)
    # Reading Input
    io = open(name_file, "a")
    fields = fieldnames(typeof(info[1]))
    num_fields = length(fields)

    # Writing Output
    write(io, "Angles\n\n")

    for angle in info
        for para in fields
            write(io, join([string(getfield(angle, para)), " "]))
        end
        write(io, "\n")
    end
    write(io, "\n\n")

    close(io)
end

end # module
