"""
    mutable struct Family_Si  <: Str

This contains information for model constructing of Si and its compounds.

# Supported list:
- Si
- SiN4
- SiN4_Ort
- SiO2
Notice:Case sensetive
"""
mutable struct Family_Si <: Str
    atom_vec::Matrix{Float64}
    cell_vec::Matrix{Float64}
    atom_type
    atom_name
    atom_charge
    atom_mass
    num_atoms
    num_atom_types
    vec_type_id
end

"""
    Si()

Do this will generate a Famliy_Si type which contains all information needed to build a Si moedel.
"""
function Si()
    atom_vec = [
        0       0       0
        0       0.5     0.5
        0.5     0       0.5
        0.5     0.5     0
        0.25    0.25    0.25
        0.25    0.75    0.75
        0.75    0.25    0.75
        0.75    0.75    0.25
    ]
    cell_vec = [
        5.4310      0           0
        0           5.4310      0
        0           0           5.4310
    ]
    atom_type = [1 1 1 1 1 1 1]
    atom_charge = 0 .* atom_type
    atom_mass = [28.085501] # Si: 28.085501. Unit: g/mol
    atom_name = split("Si")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Si(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, [1, 2])
end

"""
    Si3N4()

Do this will generate a Famliy_Si type which contains all information needed to build a Si3N4 moedel.
"""
function Si3N4()
    atom_vec = [
        -0.0977580565637542 0.195516113127508   0
        0.0588368288671674  0.307869462739661   0
        -0.0535176961032620 0.576818144016926   0
        0.327654573699065   0.151274131074550   0
        0.403949967083608   0                   0.500000000000000
        0.438645801493259   0.309541840068124   0
        0.749567567747191   0.272293435934189   0.500000000000000
        0.0180020380857845  0.645233580181295   0.500000000000000
        0.129048914298138   0.803653323477456   0.500000000000000
        0.510352849669919   0.378109310535080   0.500000000000000
        0.397866659130036   0.647057991812345   0.500000000000000
        0.554461544560957   0.759411341424498   0.500000000000000
        -0.291760564557441  0.683850293038517   0
        0.0514815980747357  0.957207969090818   0
    ]
    cell_vec = [
        7.595  0         0
        3.7975 6.577463  0
        0      0         2.902
    ]
    atom_type = [2 1 2 2 1 1 1 1 2 2 1 2 1 1]
    atom_charge = 0 .* atom_type
    atom_mass = [14.0067 28.0855] # N: 14.0067. Si: 28.085501. Unit: g/mol
    atom_name = split("N Si")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Si(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, [1, 2])
end

"""
    Si3N4_Ort()

Do this will generate a Famliy_Si type which contains all information needed to build a Si3N4 moedel with orthogonal unit cell.
"""
function Si3N4_Ort()
    atom_vec = [
    0                       0.597031356376194   0
    0.000658327847268042    0.521394509327430   0.500000000000000
    0.0308755760368664      0.401100527527099   0.500000000000000
    0.0967083607636604      0.119153727824721   0
    0.126793943383805       0                   0
    0.127583936800527       0.923220010332222   0.500000000000000
    0.146872942725477       0.363320112468571   0
    0.190125082290981       0.676164982283494   0
    0.296115865701119       0.710448618483487   0.500000000000000
    0.309479921000658       0.175330240356416   0
    0.318104015799868       0.844922570673702   0.500000000000000
    0.331599736668861       0.309804192546631   0
    0.437327188940092       0.344011811814916   0.500000000000000
    0.482422646477946       0.657540834015205   0.500000000000000
    0.500000000000000       0.0970328006978962  0
    0.500658327847268       0.0213959536491320  0.500000000000000
    0.530875576036866       0.901099083205397   0.500000000000000
    0.596708360763660       0.619152283503018   0
    0.626793943383805       0.499998555678298   0
    0.627583936800527       0.423221454653924   0.500000000000000
    0.646872942725477       0.863318668146869   0
    0.690125082290981       0.176166426605196   0
    0.796115865701119       0.210450062805189   0.500000000000000
    0.809479921000659       0.675328796034714   0
    0.818104015799869       0.344924014995404   0.500000000000000
    0.831599736668861       0.809802748224929   0
    0.937327188940092       0.844010367493214   0.500000000000000
    0.982422646477946       0.157542278336908   0.500000000000000
    ]
    cell_vec = [
    7.595   0           0
    0       13.154964   0
    0       0           2.902
    ]
    atom_type = [2 1 2 2 1 2 1 1 2 1 1 2 1 1 2 1 2 2 1 2 1 1 2 1 1 2 1 1]
    atom_charge = 0 .* atom_type
    atom_mass = [14.0067 28.0855] # N: 14.0067. Si: 28.085501. Unit: g/mol
    atom_name = split("N Si")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Si(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, [1, 2])
end

"""
    SiO2()

Do this will generate a Famliy_Si type which contains all information needed to build a SiO2 moedel.
"""
function SiO2()
    atom_vec = [
        0.196866211329851   0.196866211329851   0
        0.136601044596223   0                   0.178468624064479
        0.596625150662917   0.596625150662917   0.500000000000000
        0.656890317396545   0.793491361992768   0.678468624064479
        0.0966251506629169  0.696866211329852   0.250000000000000
        0.293491361992768   0.636601044596224   0.428468624064479
        0.696866211329852   0.0966251506629169  0.750000000000000
        0.500000000000000   0.156890317396545   0.928468624064479
        0.156890317396545   0.500000000000000   0.0715313759355211
        0.636601044596224   0.293491361992768   0.571531375935521
        0                   0.136601044596223   0.821531375935521
        0.793491361992768   0.656890317396545   0.321531375935521
    ]
    cell_vec = [
        4.9780  0           0
        0       4.9780      0
        0       0           6.9480
    ]
    atom_type = [2 1 2 1 2 1 2 1 1 1 1 1]
    atom_charge = 0 .* atom_type
    atom_mass = [15.999400 28.0855] # O: 15.999400. Si: 28.085501. Unit: g/mol
    atom_name = split("O Si")
    num_atoms = length(atom_type)
    num_atom_types = length(atom_name)
    Family_Si(atom_vec, cell_vec, atom_type, atom_name, atom_charge, atom_mass, num_atoms, num_atom_types, [1, 2])
end