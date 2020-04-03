# cat_data
"""
    cat_data(vec_data::Data...)

Do this will concatinate all model in `vec_data`, returning a variable in `Data_Sum` type
"""
function cat_data(vec_data::Data...)
    data = Data_Sum(vec_data[1])
    num_data = length(vec_data)
    if num_data == 1
        return
    else
        for id = 2 : num_data

            flag_str = 0
            # Str Info
            str_now = [typeof(data.vec_str[n]) for n = 1:length(data.vec_str)]
            if !in(typeof(vec_data[id].data_str), str_now)
                flag_str = 1  # Differnet Structure
                typ_tilt = max(data.vec_atom, "typ")
                vec_data[id].data_str.vec_type_id .+= typ_tilt
                data.vec_str = vcat(data.vec_str, vec_data[id].data_str)
                add!(vec_data[id].vec_atom, typ_tilt, "typ")
            else
                flag_str = 0 # Same Structure
                typ_tilt = data.vec_str[findall(x->x==typeof(vec_data[id].data_str), str_now)[1]].vec_type_id[1]
                typ_tilt = typ_tilt - vec_data[id].vec_atom[1].typ
                add!(vec_data[id].vec_atom, typ_tilt, "typ")
            end

            # Basic Info
            data.data_basic.num_atoms += vec_data[id].data_basic.num_atoms
            data.data_basic.num_bonds += vec_data[id].data_basic.num_bonds
            data.data_basic.num_angles += vec_data[id].data_basic.num_angles
            if flag_str == 1
                data.data_basic.num_atom_types += vec_data[id].data_basic.num_atom_types
                data.data_basic.num_bond_types += vec_data[id].data_basic.num_bond_types
                data.data_basic.num_angle_types += vec_data[id].data_basic.num_angle_types
            end
            data.data_basic.box_size[:, 2] = max(data.data_basic.box_size[:, 2], vec_data[id].data_basic.box_size[:, 2])
            data.data_basic.box_size[:, 1] = min(data.data_basic.box_size[:, 1], vec_data[id].data_basic.box_size[:, 1])
            data.data_basic.box_tilt = conv(max(data.data_basic.box_tilt, vec_data[id].data_basic.box_tilt), 1)

            atom_tilt = max(data.vec_atom, "atom")
            # Atom Info
            add!(vec_data[id].vec_atom, atom_tilt, "atom")
            data.vec_atom = vcat(data.vec_atom, vec_data[id].vec_atom)

            # Bond Info
            if typeof(data.vec_bond) == Int64
                if typeof(vec_data[id].vec_bond) != Int64
                    add!(vec_data[id].vec_bond, atom_tilt, "atom")
                    data.vec_bond = vec_data[id].vec_bond
                end
            else
                if typeof(vec_data[id].vec_bond) != Int64
                    add!(vec_data[id].vec_bond, atom_tilt, "atom")
                    data.vec_bond = vcat(data.vec_bond, vec_data[id].vec_bond)
                end
            end

            # Angle Info
            if typeof(data.vec_angle) == Int64
                if typeof(vec_data[id].vec_angle) != Int64
                    add!(vec_data[id].vec_angle, atom_tilt, "atom")
                    data.vec_angle = vec_data[id].vec_angle
                end
            else
                if typeof(vec_data[id].vec_angle) != Int64
                    add!(vec_data[id].vec_angle, atom_tilt, "atom")
                    data.vec_angle = vcat(data.vec_angle, vec_data[id].vec_angle)
                end
            end
        end
    end
    data
end

# sort_data
"""
    sort_data!(data::Data_Unit, list_atom::Array)

Do this will rearrange the atom id to consecutive one after calling `delete!`
"""
function sort_data!(data::Data_Unit, list_atom::Array)
    fields = fieldnames(typeof(data))
    name_fields = [string(fields[n]) for n = 1:length(fields)]
    num_fields = length(fields)
    list_fields = findall(x->occursin("vec", x), name_fields)
    list_atom .-= 1    # Find atom next to the deleted one
    num_atoms = length(list_atom)

    for field in list_fields
        vec_now = getfield(data, fields[field])
        if typeof(vec_now) != Int64
            len = length(vec_now)
            dims = length(vec_now[1].atom)
            vec_id = [vec_now[n].atom[i] for n = 1 : len, i = 1 : dims]
            atom_now = 1
            for id = 1 : len-1
                for dim = 1 : dims
                    diff = vec_id[id+1, dim] - vec_id[id, dim]
                    if diff >= 1
                        vec_id[1:id, dim] .+= diff - 1
                    end
                end
            end
            vec_id .-= num_atoms
            change(vec_now, vec_id, "atom")
        end
    end
    data
end


