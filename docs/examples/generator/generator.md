# 4.1 generator.jl

## example abstract
This example shows the basic pipline used in generator part of `lmp_str` to create lattice.

## Code
```julia-repl
file_out = "test.data"

data_cell = genr_cell([5, 5, 5])
str = structure_si3n4()
data = genr(data_cell, str)
write_data(data, file_out)
```
Or try run `str.jl` file in examples/generator folder. 

## Results
Results is shown below:
![result of example in generator folder](result.png)