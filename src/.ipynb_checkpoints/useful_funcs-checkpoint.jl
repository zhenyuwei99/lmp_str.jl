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

function max(array_01::AbstractArray, array_02::AbstractArray)
    if typeof(array_01).parameters[2] == 1
        num_cols = 1
        num_rows = length(array_01)
    else
        num_rows, num_cols = size(array_01)
    end
    [array_01[row, col]>array_02[row, col] ? array_01[row, col] : array_02[row, col] for row = 1:num_rows, col = 1:num_cols]
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

function min(array_01::AbstractArray, array_02::AbstractArray)
    if typeof(array_01).parameters[2] == 1
        num_cols = 1
        num_rows = length(array_01)
    else
        num_rows, num_cols = size(array_01)
    end
    [array_01[row, col]<array_02[row, col] ? array_01[row, col] : array_02[row, col] for row = 1:num_rows, col = 1:num_cols]
end

"""
    add(vec_unit::Vector{T}, tilt, para::AbstractString) where T <: Unit

Do this will add `tilt` to variable `para` for all types in `vec_unit`

# Example
```julia-repl
str = lmp_str.Si3N4()
cell = lmp_str.genr_cell([10, 10, 10])
data = genr(cell, str)
add(data.vec_bond, 3, "typ")
```
"""
function add(vec_unit::Vector{T}, tilt, para::AbstractString) where T <: Unit
    fields = fieldnames(typeof(vec_unit[1]))
    para = Meta.parse(para)
    if !in(para, fields)
        error(join(["Data_Type Atom doesn't have field: ", string(para)]))
    end
    num_units = length(vec_unit)

    for unit = 1 : num_units
        data = getfield(vec_unit[unit], para)
        len = length(data)
        if len == 1
            data += tilt
        else
            data .+= tilt
        end
        setfield!(vec_unit[unit], para, data)
    end
end

"""
    change(vec_unit::Vector{T}, tilt, para::AbstractString) where T <: Unit

Do this will change variable `para` of all types in `vec_unit` to `tilt`

# Example
```julia-repl
str = lmp_str.Si3N4()
cell = lmp_str.genr_cell([10, 10, 10])
data = genr(cell, str)
change(data.vec_bond, 3, "typ")
```
"""
function change(vec_unit::Vector{T}, goal, para::AbstractString) where T <: Unit
    fields = fieldnames(typeof(vec_unit[1]))
    para = Meta.parse(para)
    if !in(para, fields)
        error(join(["Data_Type Atom doesn't have field: ", string(para)]))
    end
    num_units = length(vec_unit)

    for unit = 1 : num_units
        data = getfield(vec_unit[unit], para)
        data = goal
        setfield!(vec_unit[unit], para, data)
    end
end

function change(vec_unit::Vector{T}, goal::AbstractArray, para::AbstractString) where T <: Unit
    fields = fieldnames(typeof(vec_unit[1]))
    para = Meta.parse(para)
    dims = size(goal)[2]
    if !in(para, fields)
        error(join(["Data_Type Atom doesn't have field: ", string(para)]))
    end
    num_units = length(vec_unit)

    for unit = 1 : num_units
        data = getfield(vec_unit[unit], para)
        if dims == 1
            data = goal[unit]
        else
            data = goal[unit, :]
        end
        setfield!(vec_unit[unit], para, data)
    end
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

Do this will return the 2-d distance between coord of `atom` and coord of `org` in plane vertical to axies `dim`.
`dim = 1, 2, 3` represents axies x, y, z respectively.

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

Do this will return the 2-d distance between coord of `pos` and coord of `org` in plane vertical to axies `dim`.
`dim = 1, 2, 3` represents axies x, y, z respectively.

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
