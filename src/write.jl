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
            write(io, join([string(atom), "\t", string(info.atom_mass[atom]), "\t\t# ", info.atom_name[atom], "\n"]))
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
    if in(:atom_mass, fields)
        write(io, "\n\nMasses\n\n")
        id = 1
        for str = 1 : length(info)
            for atom = 1 : info[str].num_atom_types
                write(io, join([string(id), "\t", string(info[str].atom_mass[atom]), "\t\t# ", info[str].atom_name[atom], "\n"]))
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

function write_info(info::Vector{T}, name_file::AbstractString) where T <: Union{Bond, Angle}
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

