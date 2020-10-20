# 4.2 converter.jl

## example abstract
This example shows the simple conversion from `.pdb` and `.psf` file created by VMD into `.data` file, which can be used in LAMMPS. The parameters for pair potential, bond, angle, dihedral, improper are all generated automatically.

## Code

```julia-repl
file_pdb = "test.pdb"
file_psf = "test.psf"
file_out = "test.data"

data = converter_vmd(file_pdb, file_psf, potential_charmm36)
write_data(data, file_out)
```

## Result
Result is shown below:
![result of example in converter folder](result.png)