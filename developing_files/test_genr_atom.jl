using lmp_str

cell = lmp_str.genr_cell([20 20 20])
str = lmp_str.Si3N4()

data = lmp_str.genr_atom(cell, str)

lmp_str.write_data(data, "test.data")
