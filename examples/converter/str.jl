using lmp_str

file_out = "test.data"
file_pdb = "test.pdb"
file_psf = "test.psf"

data = converter_vmd(file_pdb, file_psf, potential_charmm36)
write_data(data, file_out)
