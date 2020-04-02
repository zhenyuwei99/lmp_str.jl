# move
"""
    move(data::Data, move_vec::Array)

Do this will move the model uniformly respect to `move_vec` in unit of Angstrom

# Example
```julia-repl
data_cell = genr_cell([1 1 1])
str = Si3N4()
data_atom = genr_atom(data_cell, str)
data_move = move(data_atom, [10 10 10])
```
"""
function move(data::Data, move_vec)
    # Motion of box
    data_new = Data_Unit(data)
    if typeof(move_vec) == Array{Int64,2}
        move_vec = conv(move_vec, 1)
    end
    data_new.data_basic.box_size = broadcast(+, data_new.data_basic.box_size, move_vec)
    # Motion of atoms
    for atom = 1 : data_new.data_basic.num_atoms
        data_new.vec_atom[atom].coord += move_vec
    end
    data_new
end

# select and delete

function find(vec_unit::Union{Vector{T}, Int64}, list_atom) where T <: Unit
    len = length(vec_unit)
    if typeof(vec_unit) == Int64
        return 0
    else
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
end


"""
    select(data_cell::Data_Cell; mode::String, para)

Do this will return a list of id of cells in specific region controled by variables `mode` and `para`

# Arguments
- `data_cell`:
- `mode`: "cylinder" and "sphere" are supported currently.
- `para`: parameters needed to describe a region correspond to each modes

# Parameters List

- `"cylinder"`: `[radius dim]`, `dim` represents which axies cylinder will be prallel to
- `"sphere"`: `[radius]`

# Example
```julia-repl
data_cell = genr_cell([10 10 2])
data_select = select(data_cell, mode="cylinder", para=[3, 3])
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

Do this will return a list of id of atoms in specific region controled by variables `mode` and `para`

# Arguments
- `data_cell`:
- `mode`: "cylinder" and "sphere" are supported currently.
- `para`: parameters needed to describe a region correspond to each modes

# Parameters List

- `"cylinder"`: `[radius dim]`, `dim` represents which axies cylinder will be prallel to
- `"sphere"`: `[radius]`

# Example
```julia-repl
data_cell = genr_cell([10 10 2])
str = Si3N4()
data_atom = genr_atom(data_cell, str)
data_select = select(data_atom, mode="cylinder", para=[10, 3])
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

"""
    delete(data_cell::Data_Cell, list_cell::Array)

Do this will delete all cells in `list_cell` from `data_cell`

# Example
```julia-repl
data_cell = genr_cell([10 10 2])
data_select = select(data_cell, mode="cylinder", para=[3, 3])
data_new = delete(data_cell, data_select)
```
"""
function delete!(data_cell::Data_Cell, list_cell)
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

Do this will delete all atoms in `list_atom` themselves and other infomation related to them  from `data`

# Example
```julia-repl
data_cell = lmp_str.genr_cell([10 10 2])
str = Si3N4()
data_atom = genr_atom(data_cell, str)
data_select = select(data_atom, mode="cylinder", para=[3, 3])
data_new = delete(data_atom, data_select)
```
"""
function delete!(data::Data, list_atom::Array)
    fields = fieldnames(typeof(data))
    name_fields = [string(fields[n]) for n = 1:length(fields)]
    num_fields = length(fields)
    list_fields = findall(x->occursin("vec", x), name_fields)

    for field in  list_fields
        # Find all elements that need to be deleted
        list_id = find(getfield(data, fields[field]), list_atom)

        if list_id == 0
            continue
        else
            # Changing # of specifc field
            num_ids = length(list_id)
            para_now = name_fields[field][findall(x->in('_', x), name_fields[field])[1]+1 : end]
            para_now = Meta.parse(join(["num_", para_now, "s"]))
            para_result = getfield(data.data_basic, para_now) - num_ids
            setfield!(data.data_basic, para_now, para_result)

            # Changing Vector of each field
            setfield!(data, fields[field], delete(getfield(data, fields[field]), list_id)) # delete(vec::Vector{T} , id::Array) where T<:Unit
        end
    end
    sort_data(data, list_atom)
end

function delete!(vec::Union{Vector{T}, Int64}, id::Array) where T <: Unit
    len = length(vec)
    judge = trues(len)
    for i in id
        judge .&= 1 : len .!= i
    end
    vec[judge]
end

function delete!(mat; id=1, dim=1)
    if dim == 1
        result = mat[1:end .!= id, :]
    elseif dim == 2
        result = mat[:, 1:end .!= id]
    else
        error("Error, dim should be 1 or 2, representing row or column respectively!")
    end
    result
end