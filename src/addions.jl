# addions
"""
    addions(data::Data, ion_type::String, conc::Float64)

Do this will return a variable in "Data_Unit" type containing information of ions to make a solution with `conc`.

# Supported List
- Type          Mass(g/mol)
- ``Na⁺``      22.99
- ``K⁺``         39.09
- ``Al³⁺``     26.98
- ``Cl⁻``        35.453
Notcice: `ion_type` is case sensetive and space is needed for recogonition

# Example
```julia-repl
cell = genr_cell([30,20,3])
str = Tip3p()
data = genr_atom(cell, str)
data = addions(data, "Na Cl", 0.5)
```
"""
function addions(data::Data, ion_type::String, conc::Float64)
    # Supported List
    list_ion = split("K Na Fe Al Cl")
    list_charge = [1, 1, 2, 3, -1]
    list_mass = [39.09, 22.99, 55.845, 26.98, 35.453]

    # Reading Input
    ion = split(ion_type)
    len = length(ion)
    ion_id =  [findall(x->x==ion[n], list_ion)[1] for n = 1:len]
    ion_charge = [list_charge[ion_id[n]] for n = 1:len]
    ion_mass = [list_mass[ion_id[n]] for n = 1:len]

    # Genrate ion_mode
    num_unit_ions = lcm(ion_charge[1], ion_charge[2])
    num_ions = convert.(Int64, [abs(num_unit_ions/ion_charge[1]), abs(num_unit_ions/ion_charge[2])])
    num_unit_ions = sum(num_ions)

    unit_type = vcat([1 for i = 1:num_ions[1]], [2 for i = 1:num_ions[2]])
    unit_charge = vcat([ion_charge[1] for i = 1:num_ions[1]], [ion_charge[2] for i = 1:num_ions[2]])

    center = (max(data.vec_atom, "coord")-min(data.vec_atom, "coord")) ./ 2
    ion_mode = [Atom(i, i, unit_type[i], unit_charge[i], center) for i = 1:num_unit_ions]

    # Calculate # of solute molecules
    ratio_sol2wat = conc / (Const_density_wat * Const_kg2g / Const_m2dm^3 / sum(data.data_str.atom_mass))

    num_wat = data.data_basic.num_atoms
    num_sol = convert(Int64, round(num_wat * ratio_sol2wat))
    num_atoms = num_sol * num_unit_ions

    # Generate new Atoms
    vec_atom = data.vec_atom
    vec_atom_new = Vector{Atom}(undef, num_atoms)
    box_vec = data.data_basic.box_size[:,2] - data.data_basic.box_size[:,1]
    #max_mol = max(vec_atom, "mol")
    #max_type = max(vec_atom, "typ")
    for sol = 1:num_sol
        for ion = 1:num_unit_ions
            id_now = (sol-1) * num_unit_ions + ion
            atom_now = Atom(ion_mode[ion])
            atom_now.atom = id_now
            atom_now.mol = id_now
            atom_now.coord .+= (rand(3).-0.5) .* box_vec
            vec_atom_new[id_now] = atom_now
        end
    end

    # Generate new Data
    data_str = Family_Ion(ion, ion_charge, ion_mass, len, [n for n = 1:len])
    data_basic = Data_Basic(data.data_basic)
    data_basic.num_atoms = num_atoms
    data_basic.num_atom_types = len
    data_basic.num_bonds = 0
    data_basic.num_bond_types = 0
    data_basic.num_angles = 0
    data_basic.num_angle_types = 0

    # Output
    data_ion = Data_Unit(data_basic, data_str, vec_atom_new, 0, 0)
end
