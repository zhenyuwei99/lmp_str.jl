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

mutable struct Atom
    id::Int64
    mol::Int64
    atom_type::Int64
    charge::Float64
    coord::AbstractArray
end

mutable struct Bond
    id::Int64
    bond_type::Int64
    atom_01::Int64
    atom_02::Int64
end

mutable struct Angle
    id::Int64
    angle_type::Int64
    atom_01::Int64
    atom_02::Int64
    atom_03::Int64
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

function max(vec::Vector{Atom})
    max = vec[1].coord
    for i = 2 : length(vec)
        for j = 1 : 3
            if max[j] < vec[i].coord[j]
                max[j] = vec[i].coord[j]
            end
        end
    end
    max
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
    min
end

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
                cell_mat[id, :] = [x, y, z]
            end
        end
    end
    print(typeof(cell_vec))
    Data_Cell(cell_mat, cell_vec, num_cells)
end

function genr_atom(data_cell::Data_Cell, str::Str; pbc="x y z", bond_arg=1, r_cut=1.9)
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
    box_size[:, 2] = min(vec_atom) + data_cell.cell_vec*str_vec
    box_tilt = Vector{Float64}(undef, 3)
    box_tilt[1] = data_cell.cell_vec[2] * str.cell_vec[2,1]
    box_tilt[2] = data_cell.cell_vec[3] * str.cell_vec[3,1]
    box_tilt[3] = data_cell.cell_vec[2] * str.cell_vec[3,2]

    data_basic = Data_Basic(num_atoms, num_bonds, num_angles, num_dihedrals, num_impropers, num_atom_types, num_bond_types, num_angle_types, num_dihedral_types, num_improper_types, box_size, box_tilt)

    # Output
    Data_Atom(data_basic, str, vec_atom)
end

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



end # module
