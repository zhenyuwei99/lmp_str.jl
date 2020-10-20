# lmp_str.jl
***molecule builder for LAMMPS in Julia***

| **Documentation** |  **Licence** |
|:---------------------------:|:----------------------:|
|[![Documentation Status](https://readthedocs.org/projects/lmp-strjl/badge/?version=latest)](https://lmp-strjl.readthedocs.io/en/latest/?badge=latest)| [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.png?v=103)](https://opensource.org/licenses/mit-license.php)

## Requirements
The current version of `lmp_str.jl` requires **Julia 1.0** or higher.

## Quick start
### (a) Install with Julia Pkg
The package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:
```julia-repl
julia> ]
(@v1.x) pkg> add https://github.com/zhenyuwei99/lmp_str.jl.git
```
This method is convenient, while the path of `lmp_str` may take some time to found. Usually it will be `$HOME/.julia/packages/lmp_str/{version code}`. Meanwhile, files in this path are read-only.

More about Pkg: [Julia Pkg Manual](https://docs.julialang.org/en/v1/stdlib/Pkg/), [中文文档](https://cn.julialang.org/JuliaZH.jl/latest/stdlib/Pkg/#%E6%B7%BB%E5%8A%A0%E6%9C%AC%E5%9C%B0%E5%8C%85)

### (b) Install from source code
- Download the source code with git clone, not Zips!
  ```
  cd {target_path}
  git clone https://github.com/zhenyuwei99/lmp_str.git lmp_str
  ```
- Add package in julia

  ```julia-repl
  julia> ]
  (@v1.x) pkg> add "{target_path}/lmp_str"
  ```
This method require one more steps but more flexible.

