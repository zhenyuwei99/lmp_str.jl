# 2.1 Install
## Install with Julia Pkg
The package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:
```julia-repl
julia> ]
(@v1.x) pkg> add https://github.com/zhenyuwei99/lmp_str.jl.git
```
This method is convenient, while the path of `lmp_str` may take some time to found. Usually it will be `$HOME/.julia/packages/lmp_str/{version code}`. Meanwhile, files in this path are read-only.

More about Pkg: [Julia Pkg Manual](https://docs.julialang.org/en/v1/stdlib/Pkg/), [中文文档](https://cn.julialang.org/JuliaZH.jl/latest/stdlib/Pkg/#%E6%B7%BB%E5%8A%A0%E6%9C%AC%E5%9C%B0%E5%8C%85)

## Install from source code
- Download the source code with git clone, not Zips!
```sh
cd {target_path}
git clone https://github.com/zhenyuwei99/lmp_str.git lmp_str
```
- Add package in julia
```julia-repl
julia> ]
(@v1.x) pkg> add "{target_path}/lmp_str"
```
This method require one more steps but more flexible.