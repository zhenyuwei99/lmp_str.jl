"""
    function judgeequal(potential_unit::T, topo_name) where T <: Potential_Unit
This judges wether the topology (e.g bond) is same as `potential_unit`

# Example
- `judgeequal(Potential_Charmm36_Unit(["CTL3", "CL"], [200.0, 1.522]), ["CTL3", "CL"]) == true`
- `judgeequal(Potential_Charmm36_Unit(["CTL3", "CL"], [200.0, 1.522]), ["CL", "CTL3"]) == true`
- `judgeequal(Potential_Charmm36_Unit(["CTL3", "CL"], [200.0, 1.522]), ["CTL3", "CTL3"]) == false`
"""
function judgeequal(potential_unit::T, topo_name) where T <: Potential_Unit
    if length(potential_unit.flag) != length(topo_name)
        return false
    end
    
    flag_foward = true
    flag_backward = true
    for (id, atom) in enumerate(potential_unit.flag)
        if atom == "X"
            continue
        elseif atom == topo_name[id]
            continue
        else
            flag_foward = false
        end
    end

    for id = 1:length(potential_unit.flag)
        atom_now = potential_unit.flag[id]
        if atom_now == "X"
            continue
        elseif atom_now == topo_name[end-id+1]
            continue
        else
            flag_backward = false
        end
    end

    return flag_foward | flag_backward
end

function strvec2intvec(raw_str_vec)
    res = [0]
    for i in raw_str_vec
        append!(res, parse(Int, i))
    end
    return res[2:end]
end