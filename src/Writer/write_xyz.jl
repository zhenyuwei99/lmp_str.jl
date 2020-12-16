
"""
write_xyz(data::Data, name_file::AbstractString)
Do this will write all information in `data` to file `name_file` in `.xyz` formation

# Example
```julia-repl
data_cell = genr_cell([10 10 3])
str = TiO2()
data_atom = genr_atom(data_cell, str)
write_xyz(data_atom, "test.xyz")
```
"""
function write_xyz(data::Data, name_file::AbstractString)
io = open(name_file, "w")
num_atoms = data.data_basic.num_atoms
name_atoms = data.data_str.atom_name
@printf(io, "%.d\n", num_atoms)
@printf(io, ".XYZ file of %s and %s created by lmp_str \n", name_atoms[1], name_atoms[2])
for atom = 1:num_atoms
    @printf(io, "%s", name_atoms[data.vec_atom[atom].typ])
    for dim = 1:3
        @printf(io, "\t%4.4f", data.vec_atom[atom].coord[dim])
    end
    @printf(io, "\n")
end
close(io)
end