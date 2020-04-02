"""
    mutable struct Family_Basic  <: Str

This contains information for model constructing of the famliy of basic structures.

# Supported list:
 - SC (Simple Cubic)
 - BCC (Body Centered Cubic)
 - FCC (Face Centered Cubic)
 - DC (Diamond Cubic)
Notice:Case sensetive
"""
mutable struct Family_Basic <: Str
    atom_vec::Matrix{Float64}
    cell_vec::Matrix{Float64}
    atom_type
    atom_charge
    num_atoms
    num_atom_types
end

function transform(atom_vec, atom_type, atom_charge)
    num_atoms = size(atom_vec)[1]
    base = convert.(Int64, ones(num_atoms))

    if length(atom_type) == 1
        atom_type = base .* atom_type
        atom_charge = base .* atom_charge
        num_atom_types = 1
    else
        num_atom_types = length(unique(atom_type))
        if length(atom_charge) == 1
            atom_charge = base .* atom_charge
        end
    end

    atom_type, atom_charge, num_atoms, num_atom_types
end

"""
    SC(;atom_type=1, atom_charge=0, lattice_const=1)

Do this will return a Family_Basic type which contains all information needed to build a SC structure.

# Example
```julia-replr
str = SC(atom_type=1, atom_charge=0, lattice_const=3)
```
"""
function SC(;atom_type=1, atom_charge=0, lattice_const=1)
    atom_vec = [
        0   0   0
    ]
    cell_vec = diag([
        1   1   1
    ] .* lattice_const)
    atom_type, atom_charge, num_atoms, num_atom_types = transform(atom_vec, atom_type, atom_charge)
    Family_Basic(atom_vec, cell_vec, atom_type, atom_charge, num_atoms, num_atom_types)
end

"""
    BCC(;atom_type=1, atom_charge=0, lattice_const=1)

Do this will return a Family_Basic type which contains all information needed to build a BCC structure.

# Arguments
- `atom_type`: type of each atoms, could be a interger or a vector represent type of all atoms or each atoms repectively.
- `atom_charge`: charge of each atoms, could be a float or a vector, same as `atom_type` Unit: e.
- `lattice_const`: lattice constance of structure

# Example
```julia-repl
str = BCC(atom_type=1, atom_charge=-1, lattice_const=3) # All atom in cell will have same type and charge
```
Or
```julia-replr
str = BCC(atom_type=[1 2], atom_charge=[0 2], lattice_const=3)
```
"""
function BCC(;atom_type=1, atom_charge=0, lattice_const=1)
    atom_vec = [
        0   0   0
        0.5 0.5 0.5
    ]
    cell_vec = diag([
        1   1   1
    ] .* lattice_const)
    atom_type, atom_charge, num_atoms, num_atom_types = transform(atom_vec, atom_type, atom_charge)
    Family_Basic(atom_vec, cell_vec, atom_type, atom_charge, num_atoms, num_atom_types)
end

"""
    FCC(;atom_type=1, atom_charge=0, lattice_const=1)

Do this will return a Family_Basic type which contains all information needed to build a FCC structure.

# Arguments
- `atom_type`: type of each atoms, could be a interger or a vector represent type of all atoms or each atoms repectively.
- `atom_charge`: charge of each atoms, could be a float or a vector, same as `atom_type` Unit: e.
- `lattice_const`: lattice constance of structure

# Example
```julia-repl
str = FCC(atom_type=1, atom_charge=-1, lattice_const=3) # All atom in cell will have same type and charge
```
Or
```julia-replr
str = FCC(atom_type=[1 2 3 4], atom_charge=[0 0 0 0], lattice_const=3)
```
"""
function FCC(;atom_type=1, atom_charge=0, lattice_const=1)
    atom_vec = [
        0   0   0
        0   0.5 0.5
        0.5 0   0.5
        0.5 0.5 0
    ]
    cell_vec = diag([
        1   1   1
    ] .* lattice_const)
    atom_type, atom_charge, num_atoms, num_atom_types = transform(atom_vec, atom_type, atom_charge)
    Family_Basic(atom_vec, cell_vec, atom_type, atom_charge, num_atoms, num_atom_types)
end

"""
    DC(;atom_type=1, atom_charge=0, lattice_const=1)

Do this will return a Family_Basic type which contains all information needed to build a DC structure.

# Arguments
- `atom_type`: type of each atoms, could be a interger or a vector represent type of all atoms or each atoms repectively.
- `atom_charge`: charge of each atoms, could be a float or a vector, same as `atom_type` Unit: e.
- `lattice_const`: lattice constance of structure

# Example
```julia-repl
str = DC(atom_type=1, atom_charge=-1, lattice_const=3) # All atom in cell will have same type and charge
```
Or
```julia-replr
str = DC(atom_type=[1 2 3 4 5 6 7 8], atom_charge=[0 0 0 0 0 0 0 0], lattice_const=3)
```
"""
function DC(;atom_type=1, atom_charge=0, lattice_const=1)
    atom_vec = [
        0    0    0
        0    0.5  0.5
        0.5  0    0.5
        0.5  0.5  0
        0.25 0.25 0.25
        0.25 0.75 0.75
        0.75 0.25 0.75
        0.75 0.75 0.25
    ]
    cell_vec = diag([
        1   1   1
    ] .* lattice_const)
    atom_type, atom_charge, num_atoms, num_atom_types = transform(atom_vec, atom_type, atom_charge)
    Family_Basic(atom_vec, cell_vec, atom_type, atom_charge, num_atoms, num_atom_types)
end
