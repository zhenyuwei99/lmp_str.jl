# write_data and write_info
"""
    write_data(data::Data, name_file::AbstractString)

Do this will write all information in `data` to file `name_file` in `.data` formation

# Example
```julia-repl
data_cell = genr_cell([10 10 3])
str = Si3N4()
data_atom = genr_atom(data_cell, str)
write_data(data_atom, "test.data")
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
    if !(length(findall(x->x==0, box_tilt)) == length(box_tilt))
        for tilt in box_tilt
            write(io, join([string(tilt), "\t"]))
        end
        write(io, "xy xz yz\n")
    end

    close(io)
end

function write_info(info::Str, name_file::AbstractString)
    # Reading Input
    fields = fieldnames(typeof(info))
    io = open(name_file, "a")

    # Judgement and Output
    if in(:para_mass, fields)
        write(io, "\n\nMasses\n\n")
        for atom = 1 : info.num_atom_types
            write(io, join([string(atom), "\t", string(info.para_mass[atom]), "\t\t# ", info.atom_name[atom], "\n"]))
        end
    end

    if in(:para_pair, fields)
        write(io, "\n\nPair Coeffs\n\n")
        for atom = 1 : info.num_atom_types
            para_now = info.para_pair[atom, :]
            write(io, join([string(atom)]))
            for para in para_now
                write(io, join(["\t\t", string(para)]))
            end
            write(io, "\n")
        end
    end

    if in(:para_bond, fields)
        if info.num_bond_types == 0
            break
        end
        write(io, "\n\nBond Coeffs\n\n")
        for bond = 1 : info.num_bond_types
            para_now = info.para_bond[bond, :]
            write(io, join([string(bond)]))
            for para in para_now
                write(io, join(["\t\t", string(para)]))
            end
            write(io, "\n")
        end
    end

    if in(:para_angle, fields)
        if info.num_angle_types == 0
            break
        end
        write(io, "\n\nAngle Coeffs\n\n")
        for angle = 1 : info.num_angle_types
            para_now = info.para_angle[angle, :]
            write(io, join([string(angle)]))
            for para in para_now
                write(io, join(["\t\t", string(para)]))
            end
            write(io, "\n")
        end
    end

    if in(:para_dihedral, fields)
        if info.num_dihedral_types == 0
            break
        end
        write(io, "\n\nDihedral Coeffs\n\n")
        for dihedral = 1 : info.num_dihedral_types
            para_now = info.para_dihedral[dihedral, :]
            write(io, join([string(dihedral)])) 
            write(io, join(["\t\t", string(para_now[1])]))
            write(io, join(["\t\t", string(Int(para_now[2]))]))
            write(io, join(["\t\t", string(Int(para_now[3]))]))
            write(io, join(["\t\t", string(para_now[4])]))
            write(io, "\n")
        end
    end

    if in(:para_improper, fields)
        if info.num_improper_types == 0
            break
        end
        write(io, "\n\nImproper Coeffs\n\n")
        for improper = 1 : info.num_improper_types
            para_now = info.para_improper[improper, :]
            write(io, join([string(improper)]))
            for para in para_now
                write(io, join(["\t\t", string(para)]))
            end
            write(io, "\n")
        end
    end


    write(io, "\n\n")

    close(io)
end

function write_info(info::Vector{Str}, name_file::AbstractString)
    # Reading Input
    fields = fieldnames(typeof(info[1]))
    io = open(name_file, "a")

    # Judgement and Output
    if in(:para_mass, fields)
        write(io, "\n\nMasses\n\n")
        id = 1
        for str = 1 : length(info)
            for atom = 1 : info[str].num_atom_types
                write(io, join([string(id), "\t", string(info[str].para_mass[atom]), "\t\t# ", info[str].atom_name[atom], "\n"]))
                id += 1
            end
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

function write_info(info::Vector{T}, name_file::AbstractString) where T <: Union{Bond, Angle, Dihedral, Improper}
    # Reading Input
    io = open(name_file, "a")
    fields = fieldnames(typeof(info[1]))
    num_fields = length(fields)
    
    # Writing Output
    para = string(typeof(info[1]))
    write(io, join([para, "s\n\n"]))
    #para[findall(x->in('.', x), para)+1 : end]
    for item in info
        for para in 1 : num_fields - 1
            write(io, join([string(getfield(item, fields[para])), " "]))
        end
        atom = getfield(item, fields[end])
        for para in atom
            write(io, join([string(para), " "]))
        end
        write(io, "\n")
    end
    write(io, "\n\n")
    
    close(io)
end

function write_info(info::Int64, name_file::AbstractString)
    # do nothing for vec_bond or vec_angle = 0
    # It is very important although no code inside
end

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