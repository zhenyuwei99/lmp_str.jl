# Useful Functions

"""
    max(vec::Vector{Atom}, para)

Do this will return the maximum data of `para` in `vec_atom`.

# Example
```julia-repl
max_data = max(vec_atom, "coord")
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
in_data = min(vec_atom, "coord")
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
    get_data(vec_unit::Vector{T}, para::AbstractString) where T <: Data_Unit
Do this will return a vector of variable `para` for all elements in `vec_unit`.
# Example
```julia-repl
str = Si3N4()
cell = genr_cell([10, 10, 10])
data = genr(cell, str)
get_data(data.vec_angle, "atom")
```
"""
function get_data(vec_unit::Vector{T}, para::AbstractString) where T <: Unit
    fields = fieldnames(typeof(vec_unit[1]))
    para = Meta.parse(para)
    if !in(para, fields)
        error(join(["Data_Type Atom doesn't have field: ", string(para)]))
    end
    num_units = length(vec_unit)
    data = getfield(vec_unit[1], para)
    res = zeros(typeof(data[1]), num_units, length(data))

    for unit = 1 : num_units
        res[unit, :] .= getfield(vec_unit[unit], para) # try a=zeros(3, 1) a[1, :] = 1 
    end
    return res
end

"""
    add(vec_unit::Vector{T}, tilt, para::AbstractString) where T <: Data_Unit

Do this will add `tilt` to variable `para` for all elements in `vec_unit`

# Example
```julia-repl
str = Si3N4()
cell = genr_cell([10, 10, 10])
data = genr(cell, str)
add!(data.vec_bond, 3, "typ")
```
"""
function add!(vec_unit::Vector{T}, tilt, para::AbstractString) where T <: Unit
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
    change!(vec_unit::Vector{T}, tilt, para::AbstractString) where T <: Unit

Do this will change variable `para` of all elements in `vec_unit` to `tilt`

# Example
```julia-repl
str = Si3N4()
cell = genr_cell([10, 10, 10])
data = genr(cell, str)
change(data.vec_bond, 3, "typ")
```
"""
function change!(vec_unit::Vector{T}, goal, para::AbstractString) where T <: Unit
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

function change!(vec_unit::Vector{T}, goal::AbstractArray, para::AbstractString) where T <: Unit
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
r = dist(atom, [0 0 0])
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
r = dist(atom, [0 0 0], 1)
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
r = dist([1 2 3], [0 0 0], 1)
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
r = dist([1 2 3], [0 0 0], 1)
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

"""
    norm_vec(vec::Vector)
Do this will return a normalized vector along the direction of `vec`
"""
function norm_vec(vec::Vector)
    return vec ./ norm(vec)
end

"""
    function rot_mat(θx, θy, θz)
This will return a rotation matrix, coord times which will have a rotation along x, y, and z axies with angle θx, θy, θz respectively.

More details in [Wiki](https://zh.wikipedia.org/wiki/%E6%97%8B%E8%BD%AC%E7%9F%A9%E9%98%B5)
"""
function rot_mat(θx, θy, θz)
    x_mat = [
        1  0        0
        0  cos(θx)  -sin(θx)
        0  sin(θx)  cos(θx)
    ]
    
    y_mat = [
        cos(θy)  0  sin(θy)
        0        1  0
        -sin(θy) 0  cos(θy)
    ]
    
    z_mat = [
        cos(θz)  -sin(θz) 0
        sin(θz)  cos(θz)  0
        0        0        1
    ]
    
    return x_mat * y_mat * z_mat
end

"""
    function central_point_atom(data::Data)
This will return a vector of the coordinate of the center points, based on the atom infomation.
**Recommand** to use when model is small.
"""
function central_point_atom(data::Data)
	coord_vec = hcat([data.vec_atom[n].coord for n = 1:data.data_basic.num_atoms]...)
	return mean(coord_vec, dims=2)
end

"""
    function central_point_atom(data::Data)
This will return a vector of the coordinate of the center points, based on the box infomation.
**Recommand** to use when model is large.
"""
function central_point_box(data::Data)
	box_size = data.data_basic.box_size
	return (box_size[:, 2] .- box_size[:, 1])/2
end

"""
    genr_box_diag(data::Data)
This will return the matrix of box size information.
"""
function genr_box_diag(data::Data)
	box_info = data.data_basic.box_size
	box_info = box_info[:, 2] - box_info[:, 1]
	box_diag = diag(box_info)
	try
		data.data_basic.box_tilt
	catch
		return box_diag
	end
	box_diag[2, 1] = data.data_basic.box_tilt[1]
	box_diag[3, 1] = data.data_basic.box_tilt[2]
	box_diag[3, 2] = data.data_basic.box_tilt[3]
	return box_diag
end

"""
    genr_box_inv(data::Data)
This will return the inverse matrix of box size information.
"""
function genr_box_inv(data::Data)
	box_diag = genr_box_diag(data)
	return inv(box_diag)
end

function copy_array(goal)
    res = zeros(typeof(goal).parameters[1], size(goal))
    res[:] = goal[:]
    return res
end

"""
    function grid(x, y)
This will generate two array from two vector for contour.
"""
function grid(x, y)
    len_x = length(x)
    len_y = length(y)

    x_mat = zeros(len_y, len_x)
    y_mat = zeros(len_y, len_x)
    z_mat = zeros(len_y, len_x)

    for i = 1:len_y
        for j = 1:len_x
            x_mat[i, j] = x[j]
            y_mat[i, j] = y[i]
        end
    end

    return x_mat, y_mat
end

function in_str(data_str::Str, vec_str::Vector{T}) where T <: Str 
    for str in vec_str
        if str.atom_name == data_str.atom_name
            return true
        end
    end
    return false
end