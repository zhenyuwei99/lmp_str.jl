# v-0.0.5 [2020-10-09]

## Added
- READMD.md
- CHANGELOG.md
- converter series
  - `converter_vmd()`
- `Stucture_VMD`: Along with converter_vmd(), used to store potential information for `write_data()`
- Potential series
  - `potential_charmm36`: variables contain all information in Potential_Charmm36.jl
- `Angle`
- `Dihedral`
- `Improper`
- examples folder
  - generator
  - converter

## Changed
- `write_data()`: Now support writting poetential parameters for pair, bond, angle, dihedral, improper.
- Family_xxx -> Structure_xxx
  - `Family_C` -> `Structure_C`
- Xxx() -> stucture_xxx(): case sensetive -> lower case
  - `Tip3p()` -> `structure_tip3p()`