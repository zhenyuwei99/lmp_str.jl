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

"""
    Tip3p()

Do this will generate a Famliy_Wat type which contains all information needed
to build a Tip3p moedel.
"""
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
    bond_mode[1] = Bond(0, 1, [1, 2])
    bond_mode[2] = Bond(0, 1, [2, 3])
    num_bonds = length(bond_mode)
    num_bond_types = 1

    # Setting Angle Modes
    angle_mode = Vector{Angle}(undef, 1)
    angle_mode[1] = Angle(1, 1, [1, 2, 3])
    num_angles = length(angle_mode)
    num_angle_types = 1

    Family_Wat(atom_vec, cell_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types)
end

"""
    SPC()

Do this will generate a Famliy_Wat type which contains all information needed
to build a SPC moedel.
"""
function SPC()
    # Parameters of water
    density = 1 / cm2an^3;      # Unit: g/A^3
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
        ( (atom_mass[1]*2+atom_mass[2]) * gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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

    Family_Wat(atom_vec, cell_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types)
end

"""
    SPCE()

Do this will generate a Famliy_Wat type which contains all information needed
to build a SPCE moedel.
"""
function SPCE()
    # Parameters of water
    density = 1 / cm2an^3;      # Unit: g/A^3
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
        ( (atom_mass[1]*2+atom_mass[2]) * gm2g) / (density*ratio[1]*bond_len*sin(angle/2) * ratio[2]*bond_len*cos(angle/2))
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

    Family_Wat(atom_vec, cell_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types, bond_mode, num_bonds, num_bond_types, angle_mode, num_angles, num_angle_types)
end

"""
    mutable struct Family_Si  <: Str

This contains information for model constructing of Si and its compounds.

# Supported list:
- Si
- SiN4
- SiN4_Ort
- SiO2
Notice:Case sensetive
"""
mutable struct Family_Si <: Str
    atom_vec::Matrix{Float64}
    cell_vec::Matrix{Float64}
    atom_type
    atom_charge
    atom_mass
    atom_name
    num_atoms
    num_atom_types
end

"""
    Si()

Do this will generate a Famliy_Si type which contains all information needed
to build a Si moedel.
"""
function Si()
    atom_vec = [
        0       0       0
        0       0.5     0.5
        0.5     0       0.5
        0.5     0.5     0
        0.25    0.25    0.25
        0.25    0.75    0.75
        0.75    0.25    0.75
        0.75    0.75    0.25
    ]
    cell_vec = [
        5.4310      0           0
        0           5.4310      0
        0           0           5.4310
    ]
    atom_type = [1 1 1 1 1 1 1]
    atom_charge = 0 .* atom_type
    atom_mass = [28.085501] # Si: 28.085501. Unit: g/mol
    atom_name = split("Si")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Si(atom_vec, cell_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types)
end

"""
    Si3N4()

Do this will generate a Famliy_Si type which contains all information needed
to build a Si3N4 moedel.
"""
function Si3N4()
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
    cell_vec = [
        7.595  0         0
        3.7975 6.577463  0
        0      0         2.902
    ]
    atom_type = [2 1 2 2 1 1 1 1 2 2 1 2 1 1]
    atom_charge = 1.34925 .* atom_type .- 1.9275 # N: -0.57825. Si: 0.7710. Unit: e
    atom_mass = [14.0067 28.085501] # N: 14.0067. Si: 28.085501. Unit: g/mol
    atom_name = split("N Si")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Si(atom_vec, cell_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types)
end

"""
    Si3N4_Ort()

Do this will generate a Famliy_Si type which contains all information needed
to build a Si3N4 moedel with orthogonal unit cell.
"""
function Si3N4_Ort()
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
    cell_vec = [
    7.595   0           0
    0       13.154964   0
    0       0           2.902
    ]
    atom_type = [2 1 2 2 1 2 1 1 2 1 1 2 1 1 2 1 2 2 1 2 1 1 2 1 1 2 1 1]
    atom_charge = 1.34925 .* atom_type .- 1.9275 # N: -0.57825. Si: 0.7710. Unit: e
    atom_mass = [14.0067 28.085501] # N: 14.0067. Si: 28.085501. Unit: g/mol
    atom_name = split("N Si")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Si(atom_vec, cell_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types)
end

"""
    SiO2()

Do this will generate a Famliy_Si type which contains all information needed
to build a SiO2 moedel.
"""
function SiO2()
    atom_vec = [
        0.196866211329851   0.196866211329851   0
        0.136601044596223   0                   0.178468624064479
        0.596625150662917   0.596625150662917   0.500000000000000
        0.656890317396545   0.793491361992768   0.678468624064479
        0.0966251506629169  0.696866211329852   0.250000000000000
        0.293491361992768   0.636601044596224   0.428468624064479
        0.696866211329852   0.0966251506629169  0.750000000000000
        0.500000000000000   0.156890317396545   0.928468624064479
        0.156890317396545   0.500000000000000   0.0715313759355211
        0.636601044596224   0.293491361992768   0.571531375935521
        0                   0.136601044596223   0.821531375935521
        0.793491361992768   0.656890317396545   0.321531375935521
    ]
    cell_vec = [
        4.9780  0           0
        0       4.9780      0
        0       0           6.9480
    ]
    atom_type = [2 1 2 1 2 1 2 1 1 1 1 1]
    atom_charge = 1.5 .* atom_type .- 2 # N: -0.57825. Si: 0.7710. Unit: e
    atom_mass = [15.999400 28.085501] # O: 15.999400. Si: 28.085501. Unit: g/mol
    atom_name = split("O Si")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Si(atom_vec, cell_vec, atom_type, atom_charge, atom_mass, atom_name, num_atoms, num_atom_types)
end

# Type of Data
abstract type Data end
abstract type Unit end

"""
    mutable struct Atom

This contains information needed to describe an atom in Lammps Data File.
"""
mutable struct Atom <: Unit
    atom::Int64
    mol::Int64
    type::Int64
    charge::Float64
    coord::Vector{Float64}
end

"""
    Atom(atom::Atom)

Do this will generate a new variable of `Atom` type same as variable `atom`.
"""
function Atom(atom::Atom)
    id = atom.atom
    mol = atom.mol
    type = atom.type
    charge = atom.charge
    coord = Vector{Float64}(undef, 3)
    coord[:] = atom.coord[:]
    Atom(id, mol, type, charge, coord)
end

"""
    mutable struct Bond <: Unit

This contains information needed to describe a bond in Lammps Data File.
"""
mutable struct Bond <: Unit
    id::Int64
    type::Int64
    atom::Vector{Int64}
end

"""
    Bond(bond::Bond)

Do this will generate a new variable of `Bond` type same as variable `bond`.
"""
function Bond(bond::Bond)
    id = bond.id
    type = bond.type
    atom = Vector{Int64}(undef, length(bond.atom))
    atom[:] = bond.atom[:]
    Bond(id, type, atom)
end

"""
    mutable struct Angle <: Unit

This contains information needed to describe an angle in Lammps Data File.
"""
mutable struct Angle <: Unit
    id::Int64
    type::Int64
    atom::Vector{Int64}
end

"""
    Angle(angle::Angle)

Do this will generate a new variable of `Angle` type same as variable `angle`.
"""
function Angle(angle::Angle)
    id = angle.id
    type = angle.type
    atom = Vector{Int64}(undef, length(angle.atom))
    atom[:] = angle.atom[:]
    Angle(id, type, atom)
end

"""
    mutable struct Data_Cell

This contains information of cell matrix.

# Arguments
- `cell_mat`: postion of each cell point.
- `cell_vec`: # of cells in 3 direction.
- `num_cells`: # of cells in total.
"""
mutable struct Data_Cell
    cell_mat::Matrix{Int64}
    cell_vec::Vector{Int64}
    num_cells::Int64
end

"""
    mutable struct Data_Basic <:Data

This contains basic information of Lammps Model.

# Arguments
- `num_atoms`: # of atoms in total.
- ... Same formate for bonds, angles, dihedrals, and impropers.
- `num_atom_types`: # of types of atoms.
- ... Same fromate for types of bonds, angles, dihedrals, and impropers.
- `box_size`: 3 x 2 Matrix of lowest and highest postion.
- `box_tilt`: Represent xy xz yz respectively.
"""
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
    box_size::Matrix{Float64}
    box_tilt::Vector{Float64}
end

"""
    mutable struct Data_Atom <:Data

This contains all information needed for  Lammps Model with only atoms (no bond,
angles, impropers, and dihedrals).
"""
mutable struct Data_Atom <: Data
    data_basic::Data_Basic
    data_str::Str
    vec_atom::Vector{Atom}
end

"""
    mutable struct Data_Bond <:Data

This contains all information needed for Lammps Model with atoms and bonds
(no angle impropers, and dihedrals).
"""
mutable struct Data_Bond <: Data
    data_basic::Data_Basic
    data_str::Str
    vec_atom::Vector{Atom}
    vec_bond::Vector{Bond}
end

"""
    mutable struct Data_Atom <:Data

This contains all information needed for Lammps Model with atoms, bonds,
and angles (no impropers and dihedrals).
"""
mutable struct Data_Angle <: Data
    data_basic::Data_Basic
    data_str::Str
    vec_atom::Vector{Atom}
    vec_bond::Vector{Bond}
    vec_angle::Vector{Angle}
end


## Definations of Functions

# Useful Functions
"""
    max(vec::Vector{Atom}, para)

Do this will return the maximum data of `para` in `vec_atom`.

# Example
```julia-repl
max_data = lmp_str.max(vec_atom, "coord")
```
"""
function max(vec::Vector{Atom}, para)
    fields = fieldnames(typeof(vec[1]))
    para = Meta.parse(para)
    if !in(para, fields)
        error(join(["Data_Type Atom doesn't have field: ", string(para)]))
    end
    data_para = getfield(vec[1], para)
    len = length(data_para)
    if len != 1
        max = Vector{typeof(data_para).parameters[1]}(undef, len)
        max[:] = data_para[:]
        for i = 2 : length(vec)
            data_now = getfield(vec[i], para)
            for j = 1 : len
                if max[j] <= data_now[j]
                    max[j] = data_now[j]
                end
            end
        end
    else
        max = getfield(vec[1], para)
        for i = 2 : length(vec)
            data_now = getfield(vec[i], para)
            if max <= data_now
                max = data_now
            end
        end
    end
    max
end

"""
    min(vec::Vector{Atom}, para)

Do this will return the minimum data of `para` in `vec_atom`.

# Example
```julia-repl
in_data = lmp_str.min(vec_atom, "coord")
```
"""
function min(vec::Vector{Atom}, para)
    fields = fieldnames(typeof(vec[1]))
    para = Meta.parse(para)
    if !in(para, fields)
        error(join(["Data_Type Atom doesn't have field: ", string(para)]))
    end
    data_para = getfield(vec[1], para)
    len = length(data_para)
    if len != 1
        min = Vector{typeof(data_para).parameters[1]}(undef, len)
        min[:] = data_para[:]
        for i = 2 : length(vec)
            data_now = getfield(vec[i], para)
            for j = 1 : len
                if min[j] >= data_now[j]
                    min[j] = data_now[j]
                end
            end
        end
    else
        min = getfield(vec[1], para)
        for i = 2 : length(vec)
            data_now = getfield(vec[i], para)
            if min >= data_now
                min = data_now
            end
        end
    end
    min
end

function diag(vec::AbstractArray)
    len = length(vec)
    result = zeros(len, len)
    for i = 1 : len
        result[i, i] = vec[i]
    end
    result
end

function conv(vec::Vector, dim)
    len = length(vec)
    para = typeof(vec).parameters
    if dim == 1
        result = convert.(para[1], zeros(1, len))
    elseif dim == 2
        result = convert.(para[1], zeros(len, 1))
    else
        error("dim should be 1 or 2, representing Row or Column vector respectively")
    end
    for i in 1 : len
        result[i] = vec[i]
    end
    result
end

function conv(mat::Matrix, redu=0)
    para = typeof(mat).parameters
    if redu == 1
        (row, col) = size(mat)
        len = row>col ? row : col
        result = convert.(para[1], zeros(len))
    elseif redu == 0
        result = convert.(para[1], zeros(size(mat)[2], size(mat)[1]))
    else
        error("redu should be 0 or 1, representing wether Martix will be reduced to Vector.")
    end
    result[:] = mat[:]
end

"""
    dist(atom::Atom, org::Array)

Do this will return the distance between coord of `atom` and coord of `org`

# Example
```julia-repl
r = lmp_str.dist(atom, [0 0 0])
```
"""
function dist(atom::Atom, org::Array)
    r = 0
    for dim = 1 : 3
        r += (org[dim] - atom.coord[dim])^2
    end
    sqrt(r)
end

"""
    dist(atom::Atom, org::Array, dim)

Do this will return the 2-d distance between coord of `atom` and coord of `org`
in plane vertical to axies `dim`.
`dim = 1, 2, 3` represents axies x, y, z respectively

# Example
```julia-repl
r = lmp_str.dist(atom, [0 0 0], 1)
```
"""
function dist(atom::Atom, org::Array, dim)
    r = 0
    dim_range = [1 2 3]
    dim_range = dim_range[1:end .!= dim]
    for i in dim_range
        r += (org[i] - atom.coord[i])^2
    end
    sqrt(r)
end

"""
    dist(pos::Array, org::Array)

Do this will return the distance between coord of `pos` and coord of `org`

# Example
```julia-repl
r = lmp_str.dist([1 2 3], [0 0 0], 1)
```
"""
function dist(pos::Array, org::Array)
    r = 0
    for dim = 1 : 3
    r += (org[dim] - pos[dim])^2
    end
    sqrt(r)
end

"""
    dist(pos::Array, org::Array, dim)

Do this will return the 2-d distance between coord of `pos` and coord of `org`
in plane vertical to axies `dim`.
`dim = 1, 2, 3` represents axies x, y, z respectively

# Example
```julia-repl
r = lmp_str.dist([1 2 3], [0 0 0], 1)
```
"""
function dist(pos::Array, org::Array, dim)
    r = 0
    dim_range = [1 2 3]
    dim_range = dim_range[1:end .!= dim]
    for i in dim_range
        r += (org[i] - pos[i])^2
    end
    sqrt(r)
end


# genr_cell
"""
    genr_cell(cell_vec)

Do this will return a variable in `Data_Cell` type

# Example
```julia-repl
data_cell = lmp_str.genr_cell([1 2 3])
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

# genr_atom
"""
    genr_atom(data_cell::Data_Cell, str::Str)

Do this will return a variable in `Data_Atom` type

# Example
```julia-repl
data_cell = lmp_str.genr_cell([1 1 1])
str = lmp_str.Si3N4()
data_atom = lmp_str.genr_atom(data_cell, str)
```
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
    Data_Atom(data_basic, str, vec_atom)
end

# genr_bond
"""
    genr_bond(data_cell::Data_Cell, data::Data)

Do this will add information of bonds to origin model, returning a variable in
`Date_Bond` type.

# Notice
This function can only be used when `bond_mode` is contained in `data.data_str`

# Example
```julia-repl
data_cell = lmp_str.genr_cell([1 1 1])
str = lmp_str.Si3N4()
data_atom = lmp_str.genr_atom(data_cell, str)
data_bond = lmp_str.genr_bond(data_cell, data_atom)
```
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
    Data_Bond(data_basic, data_str, data.vec_atom, vec_bond)
end

# genr_angle
"""
    genr_angle(data_cell::Data_Cell, data::Data)

Do this will add information of angles to origin model, returning a variable in
`Data_Angle` type.

# Notice
This function can only be used when `angle_mode` is contained in `data.data_str`
# Example
```julia-repl
data_cell = lmp_str.genr_cell([1 1 1])
str = lmp_str.Si3N4()
data_atom = lmp_str.genr_atom(data_cell, str)
data_bond = lmp_str.genr_bond(data_cell, data_atom)
data_angle = lmp_str.genr_angle(data_cell, data_bond)
```
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
    Data_Angle(data_basic, data_str, data.vec_atom, data.vec_bond, vec_angle)
end

# move
"""
    move(data::Data, move_vec::Array)

Do this will move the model uniformly respect to `move_vec` in unit of Angstrom

# Example
```julia-repl
data_cell = lmp_str.genr_cell([1 1 1])
str = lmp_str.Si3N4()
data_atom = lmp_str.genr_atom(data_cell, str)
data_move = lmp_str.move(data_atom, [10 10 10])
```
"""
function move(data::Data, move_vec::Array)
    # Motion of box
    if typeof(move_vec) == Array{Int64,2}
        move_vec = conv(move_vec, 1)
    end
    data.data_basic.box_size = broadcast(+, data.data_basic.box_size, move_vec)

    # Motion of atoms
    for atom = 1 : data.data_basic.num_atoms
        data.vec_atom[atom].coord += move_vec
    end
    data
end

# addions
"""
    addions(data::Data, ion_type::String, conc::Float64)

Do this to addions into water model

# Supported List
- Type          Mass(g/mol)
- `Na^+`        22.99
- `K^+`         39.09
- `Al^{3+}`     26.98
- `Cl^-`        35.453
Notcice: `ion_type` is case sensetive and space is needed for recogonition

# Example
```julia-repl
cell = genr_cell([30,20,3])
str = Tip3p()
data = genr_atom(cell, str)
data = addions(data, "Na Cl", 0.5)
```
"""
function addions(data::Data, ion_type::String, conc::Float64)
    # Supported List
    list_ion = split("K Na Al Cl")
    list_charge = [1, 1, 3, -1]
    list_mass = [39.09, 22.99, 26.98, 35.453]

    # Reading Input
    ion = split(ion_type)
    len = length(ion)
    ion_id =  [findall(x->x==ion[n], list_ion)[1] for n = 1:len]
    ion_charge = [list_charge[ion_id[n]] for n = 1:len]

    # Genrate ion_mode
    num_unit_ions = lcm(ion_charge[1], ion_charge[2])
    num_ions = convert.(Int64, [abs(num_unit_ions/ion_charge[1]), abs(num_unit_ions/ion_charge[2])])
    num_unit_ions = sum(num_ions)

    unit_type = vcat([1 for i = 1:num_ions[1]], [2 for i = 1:num_ions[2]])
    unit_charge = vcat([ion_charge[1] for i = 1:num_ions[1]], [ion_charge[2] for i = 1:num_ions[2]])
    ion_mode = [Atom(i, i, unit_type[i], unit_charge[i], [0, 0, 0]) for i = 1:num_unit_ions]

    # Calculate # of solute molecules
    ratio_sol2wat = conc / (density_wat * kg2g / m2dm^3 / sum(data.data_str.atom_mass))

    num_wat = data.data_basic.num_atoms
    num_sol = convert(Int64, round(num_wat * ratio_sol2wat))
    num_atoms = num_sol * num_unit_ions

    # Generate new Atoms
    vec_atom = data.vec_atom
    vec_atom_new = Vector{Atom}(undef, num_atoms)
    box_vec = data.data_basic.box_size[:,2] - data.data_basic.box_size[:,1]
    for sol = 1:num_sol
        for ion = 1:num_unit_ions
            id_now = (sol-1) * num_unit_ions + ion
            atom_now = Atom(ion_mode[ion])
            atom_now.atom = data.data_basic.num_atoms + id_now
            atom_now.mol = max(vec_atom, "mol") + id_now
            atom_now.type += max(vec_atom, "type")
            atom_now.coord .+= rand(3) .* box_vec
            vec_atom_new[id_now] = atom_now
        end
    end

    # Generate new Data
    data.data_basic.num_atoms +=  num_atoms
    data.data_basic.num_atom_types += 2
    data.vec_atom = vcat(vec_atom, vec_atom_new)

    # Output
    data
end

# select and delete

"""
    select(data_cell::Data_Cell; mode::String, para)

Do this will return a list of id of cells in specific region controled by
variables `mode` and `para`

# Arguments
- `data_cell`:
- `mode`: "cylinder" and "sphere" are supported currently.
- `para`: parameters needed to describe a region correspond to each modes

# Parameters List

- `"cylinder"`: `[radius dim]`, `dim` represents which axies cylinder will be prallel to
- `"sphere"`: `[radius]`

# Example
```julia-repl
data_cell = lmp_str.genr_cell([10 10 2])
data_select = lmp_str.select(data_cell, mode="cylinder", para=[3, 3])
```
"""
function select(data_cell::Data_Cell; mode::String, para)
    list_mode = split("cylinder sphere")
    if !in(mode, list_mode)
        error(join(["Error, mode: ",mode," is not supported."]))
    end

    center = data_cell.cell_vec ./ 2

    cell_list = 0

    if mode == "cylinder"
        for cell in 1:data_cell.num_cells
            if dist(data_cell.cell_mat[cell, :], center, para[2]) <= para[1]
                cell_list = vcat(cell_list, cell)
            end
        end
    end

    if mode == "sphere"
        for cell in 1:data_cell.num_cells
            if dist(data_cell.cell_mat[cell, :], center) <= para[1]
                cell_list = vcat(cell_list, cell)
            end
        end
    end

    cell_list[2:end]
end

"""
    select(data_atom::Data; mode::String, para)

Do this will return a list of id of atoms in specific region controled by
variables `mode` and `para`

# Arguments
- `data_cell`:
- `mode`: "cylinder" and "sphere" are supported currently.
- `para`: parameters needed to describe a region correspond to each modes

# Parameters List

- `"cylinder"`: `[radius dim]`, `dim` represents which axies cylinder will be prallel to
- `"sphere"`: `[radius]`

# Example
```julia-repl
data_cell = lmp_str.genr_cell([10 10 2])
str = lmp_str.Si3N4()
data_atom = lmp_str.genr_atom(data_cell, str)
data_select = lmp_str.select(data_atom, mode="cylinder", para=[10, 3])
```
"""
function select(data_atom::Data; mode::String, para)
    list_mode = split("cylinder sphere")
    if !in(mode, list_mode)
        error(join(["Error, mode: ",mode," is not supported."]))
    end

    box_size = data_atom.data_basic.box_size
    box_tilt = data_atom.data_basic.box_tilt
    center = (max(data_atom.vec_atom, "coord")-min(data_atom.vec_atom, "coord")) ./ 2
    #center[1] += box_tilt[1] + box_tilt[2]
    #center[2] += box_tilt[3]

    atom_list = 0

    if mode == "cylinder"
        for atom in 1:data_atom.data_basic.num_atoms
            if dist(data_atom.vec_atom[atom], center, para[2]) <= para[1]
                atom_list = vcat(atom_list, atom)
            end
        end
    end

    if mode == "sphere"
        for atom in 1:data_atom.data_basic.num_atoms
            if dist(data_atom.vec_atom[atom], center) <= para[1]
                atom_list = vcat(atom_list, atom)
            end
        end
    end
    atom_list[2:end]
end

function sort(vec::Vector{T}) where T<:Unit
    for unit in 1 : length(vec)
        vec[unit].id = unit
    end
    vec
end

"""
    delete(data_cell::Data_Cell, list_cell::Array)

Do this will delete all cells in `list_cell` from `data_cell`

# Example
```julia-repl
data_cell = lmp_str.genr_cell([10 10 2])
data_select = lmp_str.select(data_cell, mode="cylinder", para=[3, 3])
data_new = lmp_str.delete(data_cell, data_select)
```
"""
function delete(data_cell::Data_Cell, list_cell::Array)
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

"""
    delete(data::Data, list_atom::Array)

Do this will delete all atoms in `list_atom` themselves and other infomation
related to them  from `data`

# Example
```julia-repl
data_cell = lmp_str.genr_cell([10 10 2])
str = lmp_str.Si3N4()
data_atom = lmp_str.genr_atom(data_cell, str)
data_select = lmp_str.select(data_atom, mode="cylinder", para=[3, 3])
data_new = lmp_str.delete(data_atom, data_select)
```
"""
function delete(data::Data, list_atom::Array)
    fields = fieldnames(typeof(data))
    name_fields = [string(fields[n]) for n = 1:length(fields)]
    num_fields = length(fields)
    list_fields = findall(x->occursin("vec", x), name_fields)

    for field in  list_fields
        # Find all elements that need to be deleted
        list_id = find(getfield(data, fields[field]), list_atom)

        # Changing # of specifc field
        num_ids = length(list_id)
        para_now = name_fields[field][findall(x->in('_', x), name_fields[field])[1]+1 : end]
        para_now = Meta.parse(join(["num_", para_now, "s"]))
        para_result = getfield(data.data_basic, para_now) - num_ids
        setfield!(data.data_basic, para_now, para_result)

        # Changing Vector of each field
        setfield!(data, fields[field], delete(getfield(data, fields[field]), list_id)) # delete(vec::Vector{T} , id::Array) where T<:Unit
    end

    data
end

function delete(vec::Vector{T} , id::Array) where T <: Unit
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

function find(vec_unit::Vector{T}, list_atom) where T <: Unit
    len = length(vec_unit)
    atom = [vec_unit[n].atom for n = 1:len]
    list = 0
    id = 0
    for i in atom
        id += 1
        for j in i
            if sum(list_atom .== j) == 1
                list = vcat(list, id)
                break
            end
        end
    end
    list[2:end]
end

# write_data and write_info
"""
    write_data(data::Data, name_file::AbstractString)

Do this will write all information in `data` to file `name_file`

# Example
```julia-repl
data_cell = lmp_str.genr_cell([10 10 3])
str = lmp_str.Si3N4()
data_atom = lmp_str.genr_atom(data_cell, str)
write(data_atom, "test.data")
```
"""
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

function write_info(info::Vector{T}, name_file::AbstractString) where T <: Union{Bond, Angle}
    # Reading Input
    io = open(name_file, "a")
    fields = fieldnames(typeof(info[1]))
    num_fields = length(fields)

    # Writing Output
    para = string(typeof(info[1]))
    write(io, join([para[findall(x->in('.', x), para)[1]+1 : end], "s\n\n"])) # Get "xxx" from "lmp_str.xxx"

    for bond in info
        for para in 1 : num_fields - 1
            write(io, join([string(getfield(bond, fields[para])), " "]))
        end
        atom = getfield(bond, fields[end])
        for para in atom
            write(io, join([string(para), " "]))
        end
        write(io, "\n")
    end
    write(io, "\n\n")

    close(io)
end

end # module
