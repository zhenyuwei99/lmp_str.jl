function write_pdb(data::Data, name_file::AbstractString)
    io = open(name_file, "w")
    write_pdb_CRYST1(io, data)
    write_pdb_ATOM(io, data)
    write_pdb_END(io)
    close(io)
end

function write_pdb_CRYST1(io, data::Data)
    box_size = data.data_basic.box_size[:, 2] - data.data_basic.box_size[:, 1]
    box_tilt = data.data_basic.box_tilt
    x = sqrt(box_size[1]^2 + box_tilt[1]^2 + box_tilt[2]^2)
    y = sqrt(box_size[2]^2 + box_tilt[3]^2)
    z = box_size[3]
    
    @printf(io, "CRYST1")
    @printf(io, "%9.3f", x)
    @printf(io, "%9.3f", y)
    @printf(io, "%9.3f", z)
    @printf(io, "%7.2f", 90)
    @printf(io, "%7.2f", 90)
    @printf(io, "%7.2f", 90)
    @printf(io, "\n")
end

function write_pdb_ATOM(io, data::Data_Unit)
    atom_vec = data.vec_atom
    atom_name = data.data_str.atom_name
    res_name = data.data_str.res_name
    if isa(data.data_str, Structure_Wat)
        type_name = "WT1"
        chain_id = "W"
    elseif isa(data.data_str, Structure_Ion)
        type_name = "ION"
        chain_id = "I"
    else
        type_name = "U0"
        chain_id = "U"
    end
    for (id_atom, atom) in enumerate(atom_vec)
        @printf(io, "%-6s", "ATOM")                            #  1 -  6        Record name   "ATOM  "
        @printf(io, "%5d", atom.atom)                   #  7 - 11        Integer       serial       Atom  serial number.
        @printf(io, "%2s%-3s", " ", uppercase(atom_name[atom.typ]))
        #@printf(io, "%5s", atom_name[atom.typ])        # 13 - 16        Atom          name         Atom name.
        #@printf(io, " ")                                # 17             Character     altLoc       Alternate location indicator.
        @printf(io, "%1s%-4s", " ", res_name)                    # 18 - 20        Residue name  resName      Residue name.
        @printf(io, "%1s", chain_id)                   # 22             Character     chainID      Chain identifier.
        @printf(io, "%4d", atom.mol)                    # 23 - 26        Integer       resSeq       Residue sequence number.
        @printf(io, "%4s", " ")                        # 27             AChar         iCode        Code for insertion of residues.
        @printf(io, "%8.3f", atom.coord[1])             # 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        @printf(io, "%8.3f", atom.coord[2])             # 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        @printf(io, "%8.3f", atom.coord[3])             # 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        @printf(io, "%6.2f", 1.0)                       # 55 - 60        Real(6.2)     occupancy    Occupancy.
        @printf(io, "%6.2f", 0.0)                       # 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
        @printf(io, "%6s", " ")
        @printf(io, "%-3s%3s", type_name, uppercase(atom_name[atom.typ][1:end-1]))# 73 - 78        LString(2)    element      Element symbol, right-justified.
        # 79 - 80        LString(2)    charge       Charge  on the atom.
        @printf(io, "\n")
    end
end

function write_pdb_ATOM(io, data::Data_Sum)
    atom_vec = data.vec_atom
    atom_name = []
    res_name = []
    type_name = []
    chain_name = []
    chain_id = []
    for (id, str) in enumerate(data.vec_str)
        atom_name = vcat(atom_name, str.atom_name)
        res_name = vcat(res_name, str.res_name) # Residue name
        append!(chain_id, ones(Int64, str.num_atom_types)*id)
        if isa(str, Structure_Wat)
            type_name = vcat(type_name, "WT1")
            chain_name = vcat(chain_name, "W")
        elseif isa(str, Structure_Ion)
            type_name = vcat(type_name, "ION")
            chain_name = vcat(chain_name, "I")
        else
            type_name = vcat(type_name, "U0")
            chain_name = vcat(chain_name, "U")
        end
    end
    for (id_atom, atom) in enumerate(atom_vec)
        @printf(io, "%-6s", "ATOM")                                         #  1 -  6        Record name   "ATOM  "
        @printf(io, "%5d", atom.atom)                                       #  7 - 11        Integer       serial       Atom  serial number.
        @printf(io, "%2s%-3s", " ", uppercase(atom_name[atom.typ]))         # 13 - 16        Atom          name         Atom name.
        @printf(io, "%1s%-4s", " ",res_name[chain_id[atom.typ]])            # 17 - 20        Residue name  resName      Residue name.
        @printf(io, "%1s", chain_name[chain_id[atom.typ]])                  # 22             Character     chainID      Chain identifier.
        @printf(io, "%4d", atom.mol)                                        # 23 - 26        Integer       resSeq       Residue sequence number.
        @printf(io, "%4s", " ")                                             # 27             AChar         iCode        Code for insertion of residues.
        @printf(io, "%8.3f", atom.coord[1])                                 # 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        @printf(io, "%8.3f", atom.coord[2])                                 # 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        @printf(io, "%8.3f", atom.coord[3])                                 # 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        @printf(io, "%6.2f", 0.5)                                           # 55 - 60        Real(6.2)     occupancy    Occupancy.
        @printf(io, "%6.2f", 300)                                           # 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
        @printf(io, "%6s", " ")
        @printf(io, "%-3s%3s", type_name[chain_id[atom.typ]], uppercase(atom_name[atom.typ][1:end-1]))# 73 - 78        LString(2)    element      Element symbol, right-justified.
        #@printf(io, "%9s%6s", type_name[chain_id[atom.typ]], atom_name[atom.typ])# 73 - 78        LString(2)    element      Element symbol, right-justified.
        # 79 - 80        LString(2)    charge       Charge  on the atom.
        @printf(io, "\n")
    end
end

function write_pdb_END(io)
    @printf(io, "END")
end