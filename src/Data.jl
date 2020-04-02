"""
    mutable struct Atom

This contains information needed to describe an atom in Lammps Data File.
"""
mutable struct Atom <: Unit
    atom::Int64
    mol::Int64
    typ::Int64
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
    typ = atom.typ
    charge = atom.charge
    coord = Vector{Float64}(undef, 3)
    coord[:] = atom.coord[:]
    Atom(id, mol, typ, charge, coord)
end

"""
    mutable struct Bond <: Unit

This contains information needed to describe a bond in Lammps Data File.
"""
mutable struct Bond <: Unit
    id::Int64
    typ::Int64
    atom::Vector{Int64}
end

"""
    Bond(bond::Bond)

Do this will generate a new variable of `Bond` type same as variable `bond`.
"""
function Bond(bond::Bond)
    id = bond.id
    typ = bond.typ
    atom = Vector{Int64}(undef, length(bond.atom))
    atom[:] = bond.atom[:]
    Bond(id, typ, atom)
end

"""
    mutable struct Angle <: Unit

This contains information needed to describe an angle in Lammps Data File.
"""
mutable struct Angle <: Unit
    id::Int64
    typ::Int64
    atom::Vector{Int64}
end

"""
    Angle(angle::Angle)

Do this will generate a new variable of `Angle` type same as variable `angle`.
"""
function Angle(angle::Angle)
    id = angle.id
    typ = angle.typ
    atom = Vector{Int64}(undef, length(angle.atom))
    atom[:] = angle.atom[:]
    Angle(id, typ, atom)
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

function Data_Basic(data::Data_Basic)
    box_size = [data.box_size[i, j] for i = 1 : 3, j = 1 :2]
    box_tilt = [data.box_tilt[i] for i = 1 : 3]
    Data_Basic(data.num_atoms, data.num_bonds, data.num_angles, data.num_dihedrals, data.num_impropers, data.num_atom_types, data.num_bond_types, data.num_angle_types, data.num_dihedral_types, data.num_improper_types, box_size, box_tilt)
end

"""
    mutable struct Data_Unit <:Data

This contains all information needed for an unit Lammps Model (only one structure).
"""
mutable struct Data_Unit <: Data
    data_basic::Data_Basic
    data_str::Str
    vec_atom::Union{Vector{Atom}, Int64}
    vec_bond::Union{Vector{Bond}, Int64}
    vec_angle::Union{Vector{Angle}, Int64}
end

function Data_Unit(data::Data_Unit)
    vec_atom = [Atom(data.vec_atom[n]) for n = 1 : length(data.vec_atom)]
    if typeof(data.vec_bond) != Int64
        vec_bond = [Bond(data.vec_bond[n]) for n = 1 : length(data.vec_bond)]
    else
        vec_bond = 0
    end
    if typeof(data.vec_angle) != Int64
        vec_angle = [Angle(data.vec_angle[n]) for n = 1 : length(data.vec_angle)]
    else
        vec_angle = 0
    end
    Data_Unit(Data_Basic(data.data_basic), data.data_str, vec_atom[:], vec_bond, vec_angle)
end

mutable struct Data_Sum <: Data
    data_basic::Data_Basic
    vec_str::Vector{Str}
    vec_atom::Union{Vector{Atom}, Int64}
    vec_bond::Union{Vector{Bond}, Int64}
    vec_angle::Union{Vector{Angle}, Int64}
end

function Data_Sum(data::Data_Unit)
    Data_Sum(Data_Basic(data.data_basic), [data.data_str], data.vec_atom[:], typeof(data.vec_bond)==Int64 ? data.vec_bond : data.vec_bond[:], typeof(data.vec_angle)==Int64 ? data.vec_angle : data.vec_angle[:])
end

function Data_Sum(data::Data_Sum)
    Data_Sum(Data_Basic(data.data_basic), [data.vec_str[n] for n = 1:length(data.vec_str)], data.vec_atom[:], typeof(data.vec_bond)==Int64 ? data.vec_bond : data.vec_bond[:], typeof(data.vec_angle)==Int64 ? data.vec_angle : data.vec_angle[:])
end