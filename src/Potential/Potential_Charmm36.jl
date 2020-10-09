mutable struct Potential_Charmm36_Unit <: Potential_Unit
    flag::Vector{String}
    para::Vector{Float64}
end

"""
    function Potential_Charmm36_Unit(para_info, complete_mode="normal", para_complete=0)
This will transform a string list like "HAL2   0.028     1.3400" to a `Potential_Charmm36_Unit` variable which is the formation used in `Potential_Charmm36`
- `para_info`: string list.
- `complete_mode`: 
   - `normal`: normally transformation without modification
   - `zero`: automatically fullfill zeros make `para` have the same length specified by `para_complete`
   - `copy`: automatically fullfill `para` with exist parameters in `para_info`
"""
function Potential_Charmm36_Unit(para_info, complete_mode="normal", para_complete=0)
    res = split(para_info)
    num_info = length(res)
    col_para = 0
    for info = 1:num_info
        try
            parse(Float64, String(res[info]))
            col_para = info
            break
        catch
            continue
        end
    end
    num_paras = num_info - col_para + 1  # Number of parameters
    flag = [String(res[i]) for i=1:col_para-1] # Store the atom info
    para = [parse(Float64, String(res[i])) for i=col_para:num_info]

    if complete_mode == "zero"
        if para_complete == 0
            error("`para_complete` should be specified while `zero` mode is choosen")
        elseif num_paras <= para_complete
            para = cat(para, zeros(para_complete - num_paras), dims=1)
        end
    elseif complete_mode == "copy"
        if para_complete == 0
            error("`para_complete` should be specified while `copy` mode is choosen")
        elseif num_paras < para_complete
            para = cat(para, para, dims=1)
        end
    elseif complete_mode == "one"
        if para_complete == 0
            para = cat(para, ones(1), dims=1)
        elseif num_paras < para_complete
            para = cat(para, ones(para_complete - num_paras), dims=1)
        end
    end
    
    return Potential_Charmm36_Unit(flag, para)

end

mutable struct Potential_Charmm36 <: Potential
    para_mass::Array{Potential_Charmm36_Unit}
    para_pair::Array{Potential_Charmm36_Unit}
    para_bond::Array{Potential_Charmm36_Unit}
    para_angle::Array{Potential_Charmm36_Unit}
    para_dihedral::Array{Potential_Charmm36_Unit}
    para_improper::Array{Potential_Charmm36_Unit}
end


"""
    function potential_charmm36_split_para(para_string, atom_name_vec, col_atom_end, num_cols=0)

This will convert the Charmm36 force filed data into an array which assign the atom ID according to the order of atom_name_vec.
- `para_string`: String contains the parameters of Charmm36 (like :
    "CT2A CC    200.000     1.5220 
    CT1  CS    190.000     1.5380")
- `atom_name_vec`: String vectors contains all atom names that used in the potential file.
- `col_atom_end`: The last columns of atom name in `para_string`, e.g: For bonds parameters which account two atoms interaction, atom_range=2
- `num_cols`: The # of columns of parameters, specified only while number of parameters varies between differnet rows.
"""
function Potential_Charmm36()
    para_mass = [Potential_Charmm36_Unit(info) for info in mass_para_info]
    para_pair = [Potential_Charmm36_Unit(info, "copy", 4) for info in pair_para_info]
    para_bond = [Potential_Charmm36_Unit(info) for info in bond_para_info]
    para_angle = [Potential_Charmm36_Unit(info, "zero", 4) for info in angle_para_info]
    para_dihedral = [Potential_Charmm36_Unit(info, "one") for info in dihedral_para_info]
    para_improper = [Potential_Charmm36_Unit(info) for info in improper_para_info]
    return Potential_Charmm36(para_mass, para_pair, para_bond, para_angle, para_dihedral, para_improper)
end

## Potential Info

mass_para_info = split("HL         1.00800
    HCL        1.00800
    HOL        1.00800
    HAL1       1.00800
    HAL2       1.00800
    HAL3       1.00800
    HEL1       1.00800
    HEL2       1.00800
    HBL        1.00800
    CL        12.01100
    CTL1      12.01100
    CTL2      12.01100
    CTL3      12.01100
    CTL5      12.01100
    CEL1      12.01100
    CEL2      12.01100
    CCL       12.01100
    NTL       14.00700
    NH3L      14.00700
    OBL       15.99940
    OCL       15.99940
    OSL       15.99940
    O2L       15.99940
    OHL       15.99940
    OSLP      15.99940
    PL        30.97400
    SL        32.06000
    CRL1      12.01100
    CRL2      12.01100
    H          1.00800
    HC         1.00800
    HA         1.00800
    HP         1.00800
    HB1        1.00800
    HB2        1.00800
    HR1        1.00800
    HR2        1.00800
    HR3        1.00800
    HS         1.00800
    HE1        1.00800
    HE2        1.00800
    HA1        1.00800
    HA2        1.00800
    HA3        1.00800
    C         12.01100
    CA        12.01100
    CT        12.01100
    CT1       12.01100
    CT2       12.01100
    CT2A      12.01100
    CT3       12.01100
    CPH1      12.01100
    CPH2      12.01100
    CPT       12.01100
    CY        12.01100
    CP1       12.01100
    CP2       12.01100
    CP3       12.01100
    CC        12.01100
    CD        12.01100
    CS        12.01100
    CE1       12.01100
    CE2       12.01100
    CAI       12.01100
    N         14.00700
    NR1       14.00700
    NR2       14.00700
    NR3       14.00700
    NH1       14.00700
    NH2       14.00700
    NH3       14.00700
    NC2       14.00700
    NY        14.00700
    NP        14.00700
    O         15.99900
    OB        15.99900
    OC        15.99900
    OH1       15.99900
    OS        15.99940
    S         32.06000
    SM        32.06000
    SS        32.06000
    HL         1.00800
    HCL        1.00800
    HOL        1.00800
    HAL1       1.00800
    HAL2       1.00800
    HAL3       1.00800
    HEL1       1.00800
    HEL2       1.00800
    HBL        1.00800
    CL        12.01100
    CTL1      12.01100
    CTL2      12.01100
    CTL3      12.01100
    CTL5      12.01100
    CEL1      12.01100
    CEL2      12.01100
    CCL       12.01100
    NTL       14.00700
    NH3L      14.00700
    OBL       15.99940
    OCL       15.99940
    OSL       15.99940
    O2L       15.99940
    OHL       15.99940
    OSLP      15.99940
    PL        30.97400
    SL        32.06000
    CRL1      12.01100
    CRL2      12.01100
    HN1        1.0080
    HN2        1.0080
    HN3        1.0080
    HN4        1.0080
    HN5        1.0080
    HN6        1.0080
    HN7        1.0080
    HN8        1.0080
    HN9        1.0080
    CN1       12.0110
    CN1T      12.0110
    CN2       12.0110
    CN3       12.0110
    CN3T      12.0110
    CN4       12.0110
    CN5       12.0110
    CN5G      12.0110
    CN7       12.0110
    CN7B      12.0110
    CN8       12.0110
    CN8B      12.0110
    CN9       12.0110
    NN1       14.0070
    NN2       14.0070
    NN2B      14.0070
    NN2U      14.0070
    NN2G      14.0070
    NN3       14.0070
    NN3A      14.0070
    NN3G      14.0070
    NN4       14.0070
    NN6       14.0070
    ON1       15.9994
    ON1C      15.9994
    ON2       15.9994
    ON3       15.9994
    ON4       15.9994
    ON5       15.9994
    ON6       15.9994
    ON6B      15.9994
    P         30.9740
    P2        30.9740
    P3        30.9740
    P4        30.9740
    HT    1.00800
    HX    1.00800
    OT   15.99940
    OX   15.99940
    LIT  	6.94100
    SOD  22.98977
    MG   24.30500
    POT  39.09830
    CAL  40.08000
    RUB  85.46780
    CES 132.90545
    BAR 137.32700
    ZN   65.37000
    CAD 112.41100
    CLA  35.45000
    OTMG 15.99940
    OTCA 15.99940", "\n")

pair_para_info = split("HOL      0.046     0.2245
    HAL1   0.022     1.3200
    HAL2   0.028     1.3400
    HAL3   0.024     1.3400
    HBL    0.022     1.3200
    HCL    0.046     0.2245
    HL     0.046     0.7   
    HEL1   0.031     1.25  
    HEL2   0.026     1.26  
    CL     0.0700    2.00  
    CCL    0.0700    2.00  
    CTL1   0.0200    2.275 0.01 1.9
    CTL2   0.0560    2.010 0.01 1.9
    CTL3   0.0780    2.040 0.01 1.9
    CTL5   0.0800    2.06  0.01 1.9                                          
    CEL1   0.068     2.09
    CEL2   0.064     2.08
    CRL1   0.0360    2.010 0.01 1.9
    CRL2   0.0600    2.020 0.01 1.9
    OBL    0.12      1.70  0.12 1.4
    OCL    0.12      1.70
    O2L    0.12      1.70
    OHL    0.1521    1.77
    OSL    0.1000    1.6500
    OSLP   0.1000    1.6500
    NH3L   0.20      1.85 
    NTL    0.20      1.85 
    SL     0.47      2.1  
    PL     0.585     2.15 
    HN1    0.0460    0.2245 
    HN2    0.0460    0.2245 
    HN3    0.046     1.1000
    HN4    0.0460    0.2245
    HN5    0.0460    0.2245
    HN6    0.0220    1.3200
    HN7    0.0220    1.3200
    HN8    0.0280    1.3400
    HN9    0.0240    1.3400
    NN1    -0.20     1.85
    NN2    -0.20     1.85
    NN2B   -0.20     1.85
    NN2G   -0.20     1.85
    NN2U   -0.20     1.85
    NN3    -0.20     1.85
    NN3A   -0.20     1.85 
    NN3G   -0.20     1.85 
    NN4    -0.20     1.85
    NN6    -0.20     1.85
    ON1    0.1200    1.70  
    ON1C   0.1200    1.70 
    ON2    0.1521    1.77   
    ON3    0.1200    1.70  
    ON4    0.1521    1.77  
    ON5    0.1521    1.77  
    ON6    0.1521    1.77  
    ON6B   0.1521    1.77  
    CN1    0.10      1.9000
    CN1T   0.10      1.9000
    CN2    0.10      1.9000
    CN3    0.09      1.9000
    CN3T   0.09      1.9000
    CN4    0.075     1.9000
    CN5    0.075     1.9000
    CN5G   0.075     1.9000
    CN7    0.02      2.275  0.01 1.90
    CN7B   0.02      2.275  0.01 1.90
    CN8    0.0560    2.010  0.01 1.90
    CN8B   0.0560    2.010  0.01 1.90
    CN9    0.0780    2.040  0.01 1.90
    P      0.585     2.15  
    P2     0.585     2.15  
    P3     0.585     2.15 
    P4     0.585     2.15 
    C      0.110000     2.000000      
    CA     0.070000     1.992400      
    CC     0.070000     2.000000      
    CD     0.070000     2.000000      
    CE1    0.068000     2.090000
    CE2    0.064000     2.080000	
    CP1    0.020000     2.275000   0.010000     1.900000      
    CP2    0.055000     2.175000   0.010000     1.900000      
    CP3    0.055000     2.175000   0.010000     1.900000      
    CPH1   0.050000     1.800000      
    CPH2   0.050000     1.800000      
    CS     0.110000     2.200000      
    CPT    0.099000     1.860000
    CY     0.073000     1.990000
    CAI    0.073000     1.990000      
    CT     0.0200    2.275 0.01 1.9
    CT1    0.0320    2.000 0.01 1.9
    CT2    0.0560    2.010 0.01 1.9
    CT2A   0.0560    2.010 0.01 1.9
    CT3    0.0780    2.040 0.01 1.9
    H      0.046000     0.224500      
    HA     0.022000     1.320000      
    HB1    0.022000     1.320000      
    HB2    0.028000     1.340000      
    HE1    0.031000     1.250000	
    HE2    0.026000     1.260000	      
    HC     0.046000     0.224500      
    HP     0.030000     1.358200  0.030000     1.358200      
    HR1    0.046000     0.900000
    HR2    0.046000     0.700000
    HR3    0.007800     1.468000
    HS     0.100000     0.450000
    HA1    0.045     1.3400
    HA2    0.034     1.3400
    HA3    0.024     1.3400
    N      0.200000     1.850000  0.000100     1.850000      
    NC2    0.200000     1.850000
    NH1    0.200000     1.850000  0.200000     1.550000      
    NH2    0.200000     1.850000
    NH3    0.200000     1.850000
    NP     0.200000     1.850000
    NR1    0.200000     1.850000
    NR2    0.200000     1.850000
    NR3    0.200000     1.850000
    NY     0.200000     1.850000
    O      0.120000     1.700000  0.120000     1.400000      
    OB     0.120000     1.700000  0.120000     1.400000      
    OC     0.120000     1.700000      
    OH1    0.152100     1.770000      
    OS     0.152100     1.770000      
    S      0.450000     2.000000      
    SM     0.380000     1.975000      
    SS     0.470000     2.200000
    HT     0.046     0.2245
    OT     0.1521    1.768
    OTMG   0.152100     1.76820
    OTCA   0.152100     1.76820
    OX     0.120000     1.70000  
    HX     0.046000     0.22450  
    LIT    0.00233       1.2975      
    SOD    0.0469    1.41075      
    MG     0.0150    1.18500       
    POT    0.0870    1.76375       
    CAL    0.120      1.367        
    RUB    0.15      1.90          
    CES    0.1900    2.100         
    BAR    0.150     1.890         
    ZN     0.250000     1.09000      
    CAD    0.120000     1.35700      
    CLA    0.150      2.27   ", "\n")

bond_para_info = split("CTL3  CL    200.0       1.522
    CTL2  CL    200.0       1.522
    CTL1  CL    200.0       1.522
    CTL1  CCL   200.0       1.522
    OBL   CL    750.0       1.220
    OCL   CL    525.0       1.260
    OCL   CCL   525.0       1.260
    OSL   CL    150.0       1.334
    OSLP  CL    150.0       1.334
    OHL   CL    230.0       1.40 
    HOL   OHL   545.0       0.960
    CTL1  HAL1  309.00      1.111
    CTL1  HBL   330.00      1.080
    CTL2  HAL2  309.00      1.111
    CTL3  HAL3  322.00      1.111
    CTL3  OSL   340.0       1.43 
    CTL2  OSL   340.0       1.43 
    CTL1  OSL   340.0       1.43 
    CTL3  OSLP  340.0       1.43 
    CTL2  OSLP  340.0       1.43 
    CTL1  OSLP  340.0       1.43 
    OSL   PL    270.0       1.60 
    OSLP  PL    270.0       1.60 
    O2L   PL    580.0       1.48 
    OHL   PL    237.0       1.59 
    NH3L  HCL   410.0       1.04 
    NH3L  CTL1  200.0       1.48 
    NH3L  CTL2  261.0       1.51 
    NTL   CTL2  215.00      1.51 
    NTL   CTL5  215.00      1.51 
    CTL5  HL    300.00      1.08 
    CTL2  HL    300.00      1.08 
    CTL1  CTL1  222.500     1.500
    CTL1  CTL2  222.500     1.538
    CTL1  CTL3  222.500     1.538
    CTL2  CTL2  222.500     1.530
    CTL2  CTL3  222.500     1.528
    CTL3  CTL3  222.500     1.530
    OHL   CTL1  428.0       1.420
    OHL   CTL2  428.0       1.420
    OHL   CTL3  428.0       1.420
    SL    O2L   540.0       1.448
    SL    OSL   250.0       1.575
    CEL2  CEL2  510.000     1.330
    HEL2  CEL2  365.000     1.100
    CEL1  CTL3  383.000     1.504
    CEL1  CEL2  500.000     1.342
    HEL1  CEL1  360.500     1.100
    CEL1  CTL2  365.000     1.502
    CEL1  CEL1  440.000     1.340
    CN8  NN6    200.000     1.480
    NN6  HN1    403.000     1.040
    ON6  CN8B   260.0       1.420
    CN8  CN8B   222.50      1.528
    CN1  CN3    302.0       1.409
    CN1  CN3T   302.0       1.403
    CN1  CN5G   302.0       1.360
    CN1  NN2    380.0       1.367
    CN1T NN2B   302.0       1.348
    CN1  NN2G   340.0       1.396
    CN1  NN2U   340.0       1.389
    CN1T NN2U   340.0       1.383
    CN1  NN3    350.0       1.335
    CN1T ON1    860.0       1.230
    CN1  ON1    660.0       1.234
    CN1  ON1C   620.0       1.245
    CN2  CN3    320.0       1.406
    CN2  CN5    360.0       1.358
    CN2  NN1    360.0       1.366
    CN2  NN2G   400.0       1.392
    CN2  NN3    450.0       1.343
    CN2  NN3A   400.0       1.342
    CN2  NN3G   320.0       1.326
    CN3  CN3    500.0       1.326
    CN3  CN3T   560.0       1.320
    CN3T CN9    230.0       1.478
    CN3  HN3    350.0       1.09 
    CN3T HN3    350.0       1.09 
    CN3  NN2    302.0       1.343
    CN3  NN2B   320.0       1.343
    CN4  HN3    380.0       1.09 
    CN4  NN2    320.0       1.374
    CN4  NN2B   300.0       1.378
    CN4  NN3A   420.0       1.322
    CN4  NN4    400.0       1.305
    CN5  CN5    310.0       1.361
    CN5  CN5G   320.0       1.350
    CN5  NN2    300.0       1.375
    CN5  NN2B   302.0       1.375
    CN5  NN3A   350.0       1.312
    CN5  NN3G   350.0       1.315
    CN5  NN4    310.0       1.355
    CN5G NN4    310.0       1.365
    CN8  CN8    222.50      1.528
    CN8  CN9    222.50      1.528
    CN8  NN2    400.0       1.460
    CN8  ON5    428.0       1.42 
    CN9  HN9    322.0       1.111
    CN9  ON2    340.0       1.43 
    HN1  NN1    488.0       1.00 
    HN2  NN2    474.0       1.01 
    HN2  NN2B   474.0       1.01 
    HN2  NN2G   471.0       1.01 
    HN2  NN2U   474.0       1.01 
    HN4  ON4    545.0       0.960
    ON2  P      270.0       1.60 
    ON2  P2     270.0       1.60 
    ON3  P      580.0       1.48 
    ON3  P2     580.0       1.48 
    ON4  P      237.0       1.58 
    ON4  P2     237.0       1.58 
    CN7B ON6    260.0       1.420
    CN7B CN8    200.0       1.518
    CN7  ON6    240.0       1.446
    CN7  CN7    222.5       1.529
    CN7  CN8    222.5       1.516
    CN7  CN9    222.5       1.516
    CN7  HN7    309.0       1.111
    CN8  HN8    309.0       1.111
    CN7B HN7    309.0       1.111
    CN7B ON6B   260.0       1.420
    CN7  ON6B   240.0       1.480
    CN7B CN7B   200.0       1.450
    CN7  CN7B   222.5       1.460
    CN7  CN8B   222.5       1.512
    CN8B ON2    320.0       1.44 
    CN8B ON5    428.0       1.42 
    CN7  ON2    310.0       1.433
    CN7B ON2    310.0       1.433
    CN7  ON5    428.0       1.42 
    CN9  NN2    400.0       1.456
    CN8  NN2B   400.0       1.458
    CN9  NN2B   400.0       1.458
    CN7B NN2    220.0       1.456
    CN7B NN2B   220.0       1.458
    CN8B  HN8    309.0       1.11
    ON5   HN5    545.0       0.96
    CN7B ON5    428.0       1.400
    CN8  ON2    340.0       1.44 
    ON2  P3     300.0       1.68 
    ON3  P3     480.0       1.53 
    ON2  P4     300.0       1.68 
    ON3  P4     480.0       1.53 
    ON4  P4     237.0       1.58 
    NH2  CT1   240.000     1.4550
    CA   CAI   305.000     1.3750
    CAI  CAI   305.000     1.3750
    CPT  CA    300.000     1.3600
    CPT  CAI   300.000     1.3600
    CPT  CPT   360.000     1.3850
    CY   CA    350.000     1.3650
    CY   CAI   350.000     1.3650
    CY   CPT   350.000     1.4300
    CY   CT3   375.000     1.4920
    CY   CT2   375.000     1.4920
    HP   CAI   340.000     1.0800
    HP   CY    350.000     1.0800
    NY   CA    270.000     1.3700
    NY   CPT   270.000     1.3700
    NY   H     537.500     0.9760
    CA   CA    305.000     1.3750
    CE1  CE1   440.000     1.3400	
    CE1  CE2   500.000     1.3420		
    CE1  CT2   365.000     1.5020	
    CE1  CT3   383.000     1.5040		
    CE2  CE2   510.000     1.3300		
    CP1  C     250.000     1.4900
    CP1  CC    250.000     1.4900
    CP1  CD    200.000     1.4900
    CP2  CP1   222.500     1.5270
    CP2  CP2   222.500     1.5370
    CP3  CP2   222.500     1.5370
    CPH1 CPH1  410.000     1.3600
    CT1  C     250.000     1.4900
    CT1  CC    200.000     1.5220
    CT1  CD    200.000     1.5220
    CT1  CT1   222.500     1.5000
    CT2  C     250.000     1.4900
    CT2  CA    230.000     1.4900
    CT2  CC    200.000     1.5220
    CT2  CD    200.000     1.5220
    CT2  CPH1  229.630     1.5000
    CT2  CT1   222.500     1.5380
    CT2  CT2   222.500     1.5300
    CT3  C     250.000     1.4900
    CT3  CA    230.000     1.4900
    CT3  CC    200.000     1.5220
    CT3  CD    200.000     1.5220
    CT3  CPH1  229.630     1.5000
    CT3  CS    190.000     1.5310
    CT3  CT1   222.500     1.5380
    CT3  CT2   222.500     1.5280
    CT3  CT3   222.500     1.5300
    H    CD    330.000     1.1100
    HA1  CC    317.130     1.1000
    HA2  CP2   309.000     1.1110
    HA2  CP3   309.000     1.1110
    HA2  CS    300.000     1.1110
    HA3  CS    300.000     1.1110
    HA1  CT1   309.000     1.1110
    HA2  CT2   309.000     1.1110
    HA3  CT3   322.000     1.1110
    HE1  CE1   360.500     1.1000		
    HE2  CE2   365.000     1.1000		
    HB1  CP1   330.000     1.0800
    HB1  CT1   330.000     1.0800
    HB2  CT2   330.000     1.0800
    HP   CA    340.000     1.0800
    HR1  CPH1  375.000     1.0830
    HR1  CPH2  340.000     1.0900
    HR2  CPH2  333.000     1.0700
    HR3  CPH1  365.000     1.0830
    N    C     260.000     1.3000
    N    CP1   320.000     1.4340
    N    CP3   320.000     1.4550
    NC2  C     450.000     1.3650
    NC2  CT2   390.000     1.4900
    NC2  CT3   390.000     1.4900
    NC2  HC    455.000     1.0000
    NH1  C     370.000     1.3450
    NH1  CT1   320.000     1.4300
    NH1  CT2   320.000     1.4300
    NH1  CT3   320.000     1.4300
    NH1  H     440.000     0.9970
    NH1  HC    405.000     0.9800
    NH2  CC    430.000     1.3600
    NH2  CT2   240.000     1.4550
    NH2  CT3   240.000     1.4550
    NH2  H     480.000     1.0000
    NH2  HC    460.000     1.0000
    NH3  CT1   200.000     1.4800
    NH3  CT2   200.000     1.4800
    NH3  CT3   200.000     1.4800
    NH3  HC    403.000     1.0400
    NP   CP1   320.000     1.4850
    NP   CP3   320.000     1.5020
    NP   HC    460.000     1.0060
    NR1  CPH1  400.000     1.3800
    NR1  CPH2  400.000     1.3600
    NR1  H     466.000     1.0000
    NR2  CPH1  400.000     1.3800
    NR2  CPH2  400.000     1.3200
    NR3  CPH1  380.000     1.3700
    NR3  CPH2  380.000     1.3200
    NR3  H     453.000     1.0000
    O    C     620.000     1.2300
    O    CC    650.000     1.2300
    OB   CC    750.000     1.2200
    OB   CD    750.000     1.2200
    OC   CA    525.000     1.2600
    OC   CC    525.000     1.2600
    OC   CT2   450.000     1.3300
    OC   CT3   450.000     1.3300
    OH1  CA    334.300     1.4110
    OH1  CD    230.000     1.4000
    OH1  CT1   428.000     1.4200
    OH1  CT2   428.000     1.4200
    OH1  CT3   428.000     1.4200
    OH1  H     545.000     0.9600
    OS   CD    150.000     1.3340
    OS   CT3   340.000     1.4300
    S    CT2   198.000     1.8180
    S    CT3   240.000     1.8160
    S    HS    275.000     1.3250
    SM   CT2   214.000     1.8160
    SM   CT3   214.000     1.8160
    SM   SM    173.000     2.0290
    SS   CS    205.000     1.8360
    HR1  CD    330.000     1.1100
    O    CD    720.000     1.2050
    CT2A CT1   222.500     1.5380
    CT2  CT2A  222.500     1.5300
    CT2A HA2   309.000     1.1110
    CT2A CPH1  229.630     1.5000
    CT2A CC    200.000     1.5220
    CT1  CS    190.000     1.5380
    HT    HT      0.0       1.5139
    HT    OT    450.0       0.9572
    OX    HX    545.0       0.9700
    MG      OTMG    500.00  1.9400    
    OTMG    HT     450.000     0.9572
    CAL     OTCA    500.00  1.9400    
    OTCA    HT     450.000     0.9572", "\n")

angle_para_info = split("OBL  CL   CTL3    70.0     125.0     20.0   2.442 
    OBL  CL   CTL2    70.0     125.0     20.0   2.442 
    OBL  CL   CTL1    70.0     125.0     20.0   2.442 
    OSL  CL   OBL     90.0     125.9    160.0   2.2576
    CL   OSL  CTL1    40.0     109.6     30.0   2.2651
    CL   OSL  CTL2    40.0     109.6     30.0   2.2651
    CL   OSL  CTL3    40.0     109.6     30.0   2.2651
    HAL2 CTL2 CL      33.00    109.50   30.00   2.163 
    HAL3 CTL3 CL      33.00    109.50   30.00   2.163 
    CTL2 CTL2 CL      52.0     108.00  
    CTL2 CTL1 CCL     52.0     108.00  
    CTL3 CTL2 CL      52.0     108.00  
    OSL  CL   CTL3    55.0     109.0    20.00   2.3260
    OSL  CL   CTL2    55.0     109.0    20.00   2.3260
    OSL  CL   CTL1    55.0     109.0    20.00   2.3260
    OHL  CL   OBL     50.0     123.0    210.0   2.2620
    OCL  CL   CTL2    40.0     118.0     50.0   2.3880
    OCL  CL   CTL3    40.0     118.0     50.0   2.3880
    OCL  CL   OCL    100.0     124.0     70.0   2.2250
    OCL  CCL  OCL    100.0     124.0     70.0   2.2250
    OCL  CCL  CTL1    40.0     118.0     50.0   2.3880
    OHL  CL   CTL3    55.0     110.50  
    OHL  CL   CTL2    55.0     110.50  
    HOL  OHL  CL     55.0      115.0   
    OSL  CTL1 CTL1   75.700    110.10  
    OSL  CTL1 CTL2   75.700    110.10  
    OSL  CTL1 CTL3   75.700    110.10  
    OSL  CTL2 CTL1   75.700    110.10  
    OSL  CTL2 CTL2   75.700    110.10  
    OSL  CTL2 CTL3   75.700    110.10  
    OSLP CTL1 CTL1   75.700    110.10  
    OSLP CTL1 CTL2   75.700    110.10  
    OSLP CTL1 CTL3   75.700    110.10  
    OSLP CTL2 CTL1   75.700    110.10  
    OSLP CTL2 CTL2   75.700    110.10  
    OSLP CTL2 CTL3   75.700    110.10  
    HAL2 CTL2 HAL2   35.500    109.00    5.40   1.80200
    HAL3 CTL3 HAL3   35.500    108.40    5.40   1.80200
    HAL1 CTL1 OSL    60.0      109.5   
    HAL2 CTL2 OSL    60.0      109.5   
    HAL3 CTL3 OSL    60.0      109.5   
    HAL1 CTL1 OSLP   60.0      109.5   
    HAL2 CTL2 OSLP   60.0      109.5   
    HAL3 CTL3 OSLP   60.0      109.5   
    CTL1 OSL  PL     20.0      120.0    35.0    2.33  
    CTL2 OSL  PL     20.0      120.0    35.0    2.33  
    CTL3 OSL  PL     20.0      120.0    35.0    2.33  
    CTL1 OSLP PL     20.0      120.0    35.0    2.33  
    CTL2 OSLP PL     20.0      120.0    35.0    2.33  
    CTL3 OSLP PL     20.0      120.0    35.0    2.33  
    HOL  OHL  PL     30.0      115.0    40.0    2.30  
    OSL  PL   OSL    80.0      104.3   
    OSL  PL   O2L    98.9      111.6   
    OSL  PL   OHL    48.1      108.0   
    OSLP PL   OSLP   80.0      104.3   
    OSLP PL   O2L    98.9      111.6   
    OSLP PL   OHL    48.1      108.0   
    O2L  PL   O2L   120.0      120.0   
    O2L  PL   OHL    98.9      108.23  
    NTL  CTL2 HL     40.0      109.5    27.     2.13  
    NTL  CTL5 HL     40.0      109.5    27.     2.13  
    HL   CTL2 HL     24.0      109.50   28.     1.767 
    HL   CTL5 HL     24.0      109.50   28.     1.767 
    CTL2 NTL  CTL2   60.0      109.5    26.     2.466 
    CTL5 NTL  CTL2   60.0      109.5    26.     2.466 
    CTL5 NTL  CTL5   60.0      109.5    26.     2.466 
    HL   CTL2 CTL2   33.430    110.10   22.53   2.179 
    HL   CTL2 CTL3   33.430    110.10   22.53   2.179 
    HAL1 CTL1 CTL1   34.500    110.10   22.53   2.179 
    HAL1 CTL1 CTL2   34.500    110.10   22.53   2.179 
    HAL1 CTL1 CTL3   34.500    110.10   22.53   2.179 
    HAL2 CTL2 CTL1   26.500    110.10   22.53   2.179 
    HAL2 CTL2 CTL2   26.500    110.10   22.53   2.179 
    HAL2 CTL2 CTL3   34.600    110.10   22.53   2.179 
    HAL3 CTL3 CTL1   33.430    110.10   22.53   2.179 
    HAL3 CTL3 CTL2   34.600    110.10   22.53   2.179 
    HAL3 CTL3 CTL3   37.500    110.10   22.53   2.179 
    HBL  CTL1 CCL    50.000    109.50  
    HBL  CTL1 CTL2   35.000    111.00  
    NTL  CTL2 CTL2   67.7      115.00  
    NTL  CTL2 CTL3   67.7      115.00  
    HCL  NH3L CTL2   33.0      109.50    4.00   2.056 
    HCL  NH3L CTL1   30.0      109.50   20.00   2.074 
    HCL  NH3L HCL    41.0      109.50  
    NH3L CTL2 CTL2   67.7      110.00  
    NH3L CTL2 HAL2   45.0      107.50   35.00   2.0836
    CTL1 CTL1 CTL1   53.350    111.00    8.00   2.561 
    NH3L CTL1 CCL    43.7      110.00  
    NH3L CTL1 CTL2   67.7      110.00  
    NH3L CTL1 HBL    51.5      107.50  
    CTL1 CTL1 CTL2   58.350    113.50   11.16   2.561 
    CTL1 CTL1 CTL3   53.350    108.50    8.00   2.561 
    CTL1 CTL2 CTL1   58.350    113.50   11.16   2.561 
    CTL1 CTL2 CTL2   58.350    113.50   11.16   2.561 
    CTL1 CTL2 CTL3   58.350    113.50   11.16   2.561 
    CTL2 CTL1 CTL2   58.350    113.50   11.16   2.561 
    CTL2 CTL1 CTL3   58.350    113.50   11.16   2.561 
    CTL2 CTL2 CTL2   58.350    113.60   11.16   2.561 
    CTL2 CTL2 CTL3   58.000    115.00    8.00   2.561 
    CTL3 CTL1 CTL3   58.350    113.50   11.16   2.561 
    HOL  OHL  CTL1   57.500    106.00  
    HOL  OHL  CTL2   57.500    106.00  
    HOL  OHL  CTL3   57.500    106.00  
    OHL  CTL1 CTL1   75.700    110.10  
    OHL  CTL1 CTL2   75.700    110.10  
    OHL  CTL2 CTL1   75.700    110.10  
    OHL  CTL2 CTL2   75.700    110.10  
    OHL  CTL2 CTL3   75.700    110.10  
    OHL  CTL1 HAL1   45.900    108.89  
    OHL  CTL2 HAL2   45.900    108.89  
    OHL  CTL3 HAL3   45.900    108.89  
    O2L  SL   O2L   130.0      109.47  35.0    2.45
    O2L  SL   OSL    85.0       98.0           
    CTL2 OSL  SL     15.0      109.0   27.00   1.90
    CTL3 OSL  SL     15.0      109.0   27.00   1.90
    CEL1 CEL1 CTL2   48.00     123.50  
    CEL1 CEL1 CTL3   48.00     123.50  
    CEL2 CEL1 CTL2   48.00     126.00  
    CEL2 CEL1 CTL3   47.00     125.20  
    HEL1 CEL1 CEL1   52.00     119.50  
    HEL1 CEL1 CEL2   42.00     118.00  
    HEL1 CEL1 CTL2   40.00     116.00  
    HEL1 CEL1 CTL3   22.00     117.00  
    HEL2 CEL2 CEL1   45.00     120.50  
    HEL2 CEL2 CEL2   55.50     120.50  
    HEL2 CEL2 HEL2   19.00     119.00  
    CEL1 CTL2 CTL2   32.00     112.20  
    CEL1 CTL2 CTL3   32.00     112.20  
    HAL2 CTL2 CEL1   45.00     111.50  
    HAL3 CTL3 CEL1   42.00     111.50  
    CEL1 CTL2 CEL1   30.0      114.0   
    CN7  CN8  CN8      58.35    113.60   11.16   2.561
    CN8  CN7  CN8      58.35    113.60   11.16   2.561
    CN8  CN8  CN8      58.35    113.60   11.16   2.561
    HN1  NN6  CN8      30.00    109.50   20.00   2.074
    NN6  CN8  HN8      45.00    107.50   35.00   2.101
    CN7  CN8  ON2     115.00    109.70 
    NN6  CN8  CN8      67.70    110.00 
    HN1  NN6  HN1      44.00    109.50 
    ON2  CN8  CN8     115.0     109.7  
    CN7  ON6  CN8B    110.0     109.0
    ON6  CN8B CN8      90.0     106.0
    CN8B CN8  CN7      80.0     106.0
    ON6  CN8B HN8      45.2     107.24
    HN8  CN8B CN8      34.53    110.10  22.53   2.179
    HN8  CN8  CN8B     34.53    110.10  22.53   2.179
    CN2  NN3A CN4     90.0     117.8 
    NN3A CN4  NN3A    60.0     133.0 
    CN4  NN3A CN5     90.0     110.1 
    CN5  CN5  NN3A    60.0     127.4 
    CN2  CN5  CN5     60.0     121.0 
    CN5  CN2  NN3A    60.0     110.7 
    CN5  CN5  NN2    100.0     105.7 
    CN5  CN5  NN4    100.0     110.0 
    CN4  NN4  CN5    120.0     104.6 
    NN2  CN4  NN4    100.0     113.4 
    CN4  NN2  CN5    100.0     106.3 
    NN2  CN5  NN3A   100.0     126.9 
    CN2  CN5  NN4    100.0     129.0 
    HN3  CN4  NN3A    38.0     113.5 
    NN3A CN2  NN1     50.0     130.7 
    CN5  CN2  NN1     50.0     118.6 
    CN2  NN1  HN1     40.0     121.5 
    HN1  NN1  HN1     31.0     117.0 
    NN4  CN4  HN3     39.0     124.8 
    NN2  CN4  HN3     39.0     121.8 
    CN5  NN2  HN2     30.0     129.4 
    CN4  NN2  HN2     30.0     125.0 
    CN1  NN2G CN2     70.0     131.1 
    NN2G CN2  NN3G    70.0     122.2 
    CN2  NN3G CN5     90.0     109.4 
    CN5G CN5  NN3G    70.0     129.9 
    CN1  CN5G CN5     70.0     119.6 
    CN5G CN1  NN2G    70.0     107.8 
    CN5G CN5  NN2B   100.0     104.6 
    CN5  CN5G NN4    100.0     111.4 
    CN4  NN4  CN5G   120.0     103.8 
    NN2B CN4  NN4    100.0     113.0 
    CN4  NN2B CN5    100.0     107.2 
    NN2B CN5  NN3G   140.0     125.5 
    CN1  CN5G NN4    125.0     129.0 
    CN1  NN2G HN2     45.0     113.3 
    CN2  NN2G HN2     45.0     115.6 
    NN1  CN2  NN2G    95.0     115.4 
    NN1  CN2  NN3G    95.0     122.4 
    NN2G CN1  ON1     50.0     127.5 
    CN5G CN1  ON1     50.0     124.7 
    HN3  CN4  NN2B    40.0     122.2 
    CN4  NN2B HN2     30.0     124.6 
    CN5  NN2B HN2     30.0     129.3 
    CN1  NN2  CN3     50.0     124.1 
    NN2  CN1  NN3     50.0     116.8 
    CN1  NN3  CN2     85.0     119.1 
    CN3  CN2  NN3     85.0     119.3 
    CN2  CN3  CN3     85.0     117.8 
    CN3  CN3  NN2     85.0     122.9 
    CN1  NN2  HN2     37.0     121.2 
    CN3  NN2  HN2     37.0     114.7 
    NN2  CN1  ON1C   130.0     119.4 
    NN3  CN1  ON1C   130.0     123.8 
    NN3  CN2  NN1     81.0     122.3 
    CN3  CN2  NN1     81.0     118.4 
    CN2  CN3  HN3     38.0     120.1 
    CN3  CN3  HN3     38.0     122.1 
    HN3  CN3  NN2     44.0     115.0 
    CN1T NN2B CN3     70.0     122.0 
    NN2B CN1T NN2U    50.0     114.0 
    CN1T NN2U CN1     50.0     130.2 
    NN2U CN1  CN3     70.0     112.6 
    CN1  CN3  CN3    100.0     117.6 
    CN3  CN3  NN2B   100.0     123.6 
    CN1T NN2B HN2     40.5     122.0 
    CN3  NN2B HN2     32.0     116.0 
    NN2B CN1T ON1    100.0     121.6 
    NN2U CN1T ON1    100.0     124.4 
    CN1T NN2U HN2     40.5     114.4 
    CN1  NN2U HN2     40.5     115.4 
    NN2U CN1  ON1    100.0     121.9 
    CN3  CN1  ON1    100.0     125.5 
    CN1  CN3  HN3     30.0     120.3 
    HN3  CN3  NN2B    30.0     114.3 
    CN3T CN1 NN2U     70.0     113.5 
    CN1  CN3T CN3    120.0     116.7 
    CN3T CN3  NN2B   120.0     123.6 
    CN3T CN1  ON1    100.0     124.6 
    CN1  CN3T CN9     38.0     118.7 
    CN3  CN3T CN9     38.0     124.6 
    CN3T CN3  HN3     30.0     122.1 
    CN1T NN2B CN9     70.0     116.0 
    CN3  NN2B CN9     70.0     122.0 
    CN1  NN2  CN9     70.0     115.4 
    CN3  NN2  CN9     70.0     120.5 
    CN5  NN2  CN9     70.0     125.9 
    CN4  NN2  CN9     70.0     127.8 
    CN5  NN2B CN9     70.0     125.9 
    CN4  NN2B CN9     70.0     126.9 
    CN5  NN2B CN8     70.0     125.9 
    CN4  NN2B CN8     70.0     126.9 
    NN2B CN8  CN9     70.0     113.7 
    CN1T NN2B CN7B    45.0     118.4 
    CN3  NN2B CN7B    45.0     119.6 
    CN1  NN2  CN7B    45.0     120.0 
    CN3  NN2  CN7B    45.0     115.9 
    CN5  NN2  CN7B    45.0     126.1 
    CN4  NN2  CN7B    45.0     127.6 
    CN5  NN2B CN7B    45.0     126.5 
    CN4  NN2B CN7B    45.0     126.3 
    ON6  CN7B NN2    110.0     108.0 
    ON6B CN7B NN2    110.0     112.0 
    CN8  CN7B NN2    110.0     113.7 
    CN7B CN7B NN2    110.0     111.0 
    ON6  CN7B NN2B   110.0     108.0 
    ON6B CN7B NN2B   110.0     112.0 
    CN8  CN7B NN2B   110.0     113.7 
    CN7B CN7B NN2B   110.0     111.0 
    HN7  CN7B NN2     43.0     111.0 
    HN7  CN7B NN2B    43.0     111.0 
    CN9  CN8  HN8     34.6     110.10  22.53 2.179
    CN9  CN7  HN7     34.6     110.10  22.53 2.179
    HN8  CN8  NN2     33.43    110.1 
    HN8  CN8  ON5     45.9     108.89
    CN3  CN9  HN9     33.43    110.10  22.53 2.179
    CN3T CN9  HN9     33.43    110.10  22.53 2.179
    CN8  CN9  HN9     34.60    110.10  22.53 2.179
    HN9  CN9  CN7     33.43    110.1   22.53 2.179
    HN9  CN9  NN2     33.43    110.1 
    HN9  CN9  NN2B    33.43    110.1 
    HN8  CN8  NN2B    33.43    110.1 
    HN9  CN9  ON2     60.0     109.5 
    CN9  ON2  P       20.0     120.0  35.   2.33
    HN4  ON4  P       30.0     115.0  40.0  2.35
    HN4  ON4  P2      30.0     115.0  40.0  2.35
    HN5  ON5  CN8     57.5     106.0 
    HN5  ON5  CN9     57.5     106.0 
    ON2  P    ON2     80.0     104.3 
    ON2  P2   ON2     80.0     104.3 
    ON2  P    ON4     48.1     108.0 
    ON2  P2   ON4     48.1     108.0 
    ON3  P    ON4     98.9     108.23
    ON3  P2   ON4     98.9     108.23
    ON4  P    ON4     98.9     104.0 
    ON4  P2   ON4     98.9     104.0 
    CN7  CN8  ON5     75.7     110.10
    HN9  CN9  HN9     35.500   108.40  5.40 1.802
    CN7  ON6  CN7B    110.0    108.0  
    ON6  CN7B CN8      90.0    102.0  
    CN7B CN8  CN7      80.00   100.0  
    CN8  CN7  CN7      60.00   102.0   8.0   2.561 
    CN9  CN7  CN7      60.00   102.0   8.0   2.561 
    CN7  CN7  ON6     100.0    104.0  
    HN7  CN7  ON6      45.2    107.24 
    HN7  CN7B ON6      45.2    107.24 
    HN7  CN7  CN7      40.0    108.00 
    CN7B CN8  HN8      33.4    110.10  22.53   2.179
    CN8  CN7B HN7      33.4    110.10  22.53   2.179
    HN7  CN7  CN8      34.5    110.1   22.53   2.179
    HN8  CN8  CN7      34.53   110.10  22.53   2.179
    HN8  CN8  CN8      34.53   110.10  22.53   2.179
    HN8  CN8  HN8      35.5    109.00   5.40   1.802
    HN7  CN7  HN7      35.5    109.00   5.40   1.802
    CN7  ON6B CN7B    110.0    115.0  
    CN7  CN7  ON6B    100.0    110.0  
    ON6B CN7B CN7B     90.0    106.0  
    CN7B CN7B CN7     110.0     96.0  
    CN7B CN7  CN7      60.0    100.0    8.00   2.561 
    HN7  CN7  ON6B     45.2    107.24 
    HN7  CN7B ON6B     45.2    107.24 
    CN7B CN7B HN7      33.4    110.10  22.53   2.179
    HN7  CN7B HN7      35.5    109.00   5.40   1.802
    ON6  CN7  CN8B     90.0    108.2 
    ON6  CN7  CN9      90.0    108.2 
    CN7  CN7  CN8B     45.0    110.0 
    CN8  CN7  CN8B    58.35    113.60 11.16 2.561
    CN7  CN8B ON2      70.0    108.4 
    CN7  CN7  ON2     115.0    109.7 
    CN7B CN7B ON2     115.0    109.7 
    CN8  CN7  ON2     115.0    109.7 
    CN8B ON2  P        20.0    120.0   35.00   2.33
    CN8B ON2  P2       20.0    120.0   35.00   2.33
    CN7  ON2  P        20.0    120.0   35.00   2.33
    CN7  ON2  P2       20.0    120.0   35.00   2.33
    CN7B ON2  P        20.0    120.0   35.00   2.33
    CN7B ON2  P2       20.0    120.0   35.00   2.33
    HN7  CN7  CN8B     34.5    110.1   22.53   2.179
    HN8  CN8B ON2      60.0    109.5 
    HN5  ON5  CN8B     57.5    106.0 
    HN8  CN8B HN8      35.5    109.0    5.40   1.802 
    HN8  CN8B CN7      34.53   110.1   22.53   2.179 
    HN7  CN7  ON2      60.0    109.5 
    HN7  CN7B ON2      60.0    109.5 
    CN7  CN8B ON5      75.7    110.10
    CN8B CN7  ON5      90.0    108.2 
    HN8  CN8B ON5      45.9    108.89
    ON5  CN7  CN8      75.7    110.0 
    ON5  CN7  CN7      75.7    110.1 
    HN7  CN7  ON5      60.0    109.5 
    HN5  ON5  CN7      57.5    109.0 
    ON6B CN7  CN8B     90.0    108.2 
    ON6B CN7  CN9      90.0    108.2 
    ON2  CN7  CN7B     90.0    110.0 
    ON5  CN7  CN7B     90.0    110.0 
    ON5  CN7B CN7B     80.0    108.4 
    ON5  CN7B CN7      90.0    108.0 
    HN7  CN7B ON5      60.0    109.5 
    HN5  ON5  CN7B     57.5    109.0 
    HN7  CN7B CN7      34.53   110.10  22.53   2.179
    HN7  CN7  CN7B     34.53   110.10  22.53   2.179
    CN8  ON2  P        20.0    120.0  35.   2.33
    CN8  ON2  P2       20.0    120.0  35.   2.33
    ON2  P    ON3      98.9    111.6 
    ON2  P2   ON3      98.9    111.6 
    ON3  P    ON3     120.0    120.0 
    ON3  P2   ON3     120.0    120.0 
    HN8  CN8  ON2      60.0    109.5 
    ON5  P    ON3      98.9    111.6 
    ON6  CN7B CN7     120.0     106.25  
    CN7B CN7  CN8      58.35    113.6   11.16 2.561
    P    ON2  P       15.0     140.0 -40.0  2.80
    CN9  ON2  P3      20.0     120.0  35.   2.33
    P3   ON2  P       15.0     140.0 -40.0  2.80
    P4   ON2  P       45.0     140.0  40.0  3.15
    P3   ON2  P3      15.0     140.0 -40.0  2.80
    HN4  ON4  P4      30.0     120.0  40.0  2.35
    ON2  P4   ON4     48.1     100.0 
    ON3  P4   ON4     98.9     108.23
    ON2  P3   ON2     80.0     104.3 
    ON2  P3   ON3     88.9     111.6 
    ON2  P4   ON3     88.9     105.0 
    ON3  P3   ON3    104.0     120.0 
    ON3  P4   ON3    104.0     120.0 
    H    NH2  CT1   50.000    111.00         
    H    NH2  CT2   50.000    111.00         
    NH2  CT1  CT1   67.700    110.00         
    NH2  CT1  CT2   67.700    110.00         
    NH2  CT1  CT2A  67.700    110.00         
    NH2  CT1  CT3   67.700    110.00         
    CT1  CD   OH1   55.000    110.50         
    CT3  CT1  CD    52.000    108.00         
    NH2  CT1  HB1   38.000    109.50   50.00   2.1400
    NH2  CT1  C     50.000    107.00         
    NH2  CT1  CD    50.000    107.00         
    NH2  CT2  C     50.000    107.00         
    HB2  CT1  HB2   36.000    115.00         
    HB2  CT1  CD    50.000    109.50         
    NH1  CT1  HB2   48.000    108.00         
    CAI  CAI  CA    40.000    120.00   35.00   2.41620
    CAI  CA   CA    40.000    120.00   35.00   2.41620
    CPT  CA   CA    50.000    113.20
    CPT  CPT  CA    50.000    110.00
    CPT  CAI  CA    50.000    113.20
    CPT  CPT  CAI   50.000    110.00
    CPT  CY   CA    85.000    106.40   25.00   2.26100
    CPT  NY   CA    85.000    112.00
    CT2  CY   CA    30.000    127.00
    CT2  CY   CPT   30.000    126.70
    CT3  CY   CA    30.000    127.00
    CT3  CY   CPT   30.000    126.70
    CY   CPT  CA   130.000    133.50
    CY   CPT  CAI  130.000    133.50
    CY   CPT  CPT   85.000    108.00
    CY   CT2  CT1   58.350    114.00
    CY   CT2  CT3   58.350    114.00
    H    NY   CA    28.000    126.00
    H    NY   CAI   28.000    126.00
    H    NY   CPT   28.000    126.00
    HA2  CT2  CY    55.000    109.50
    HA3  CT3  CY    55.000    109.50
    HP   CA   CAI   30.000    120.00   22.00   2.15250
    HP   CAI  CA    30.000    120.00   22.00   2.15250
    HP   CA   CPT   30.000    122.00   22.00   2.14600
    HP   CAI  CPT   30.000    122.00   22.00   2.14600
    HP   CA   CY    32.000    125.00   25.00   2.17300
    HP   CY   CA    32.000    126.40   25.00   2.18600
    HP   CY   CPT   32.000    126.40   25.00   2.25500
    NY   CA   CY    85.000    110.50   25.00   2.24000
    NY   CA   HP    32.000    125.00   25.00   2.17700
    NY   CPT  CA   130.000    129.50
    NY   CPT  CAI  130.000    129.50
    NY   CPT  CPT   95.000    107.40
    CA   CA   CA    40.000    120.00   35.00   2.41620
    CE1  CE1  CT2    48.00    123.50  
    CE1  CE1  CT3    48.00    123.50  	
    CE1  CT2  CT3    32.00    112.20  	
    CE2  CE1  CT2    48.00    126.00  	
    CE2  CE1  CT3    47.00    125.20  
    CP1  N    C      60.000   117.0000
    CP2  CP1  C      52.000   112.3000
    CP2  CP1  CC     52.000   112.3000
    CP2  CP1  CD     50.000   112.3000
    CP2  CP2  CP1    70.000   108.5000
    CP3  CP2  CP2    70.000   108.5000
    CP3  N    C      60.000   117.0000
    CP3  N    CP1   100.000   114.2000
    CP3  NP   CP1   100.000   111.0000
    CPH2 NR1  CPH1  130.000   107.5000
    CPH2 NR2  CPH1  130.000   104.0000
    CPH2 NR3  CPH1  145.000   108.0000
    CT1  CT1  C      52.000   108.0000
    CT1  CT1  CC     52.000   108.0000
    CT1  CT1  CD     52.000   108.0000
    CT1  CT1  CT1   53.350    111.00    8.00   2.56100
    CT1  CT2  CA     51.800   107.5000
    CT1  CT2  CC     52.000   108.0000
    CT1  CT2  CD     52.000   108.0000
    CT1  CT2  CPH1   58.350   113.0000
    CT1  CT2  CT1   58.350    113.50   11.16   2.56100
    CT1  NH1  C      50.000   120.0000
    CT2  CA   CA     45.800   122.3000
    CT2  CPH1 CPH1   45.800   130.0000
    CT2  CT1  C      52.000   108.0000
    CT2  CT1  CC     52.000   108.0000
    CT2A CT1  CC     52.000   108.0000
    CT2  CT1  CD     52.000   108.0000
    CT2  CT1  CT1   53.350    111.00    8.00   2.56100
    CT2  CT2  C      52.000   108.0000
    CT2  CT2  CC     52.000   108.0000
    CT3  CT2  CC     52.000   108.0000
    CT2  CT2  CD     52.000   108.0000
    CT2A CT2  CD     52.000   108.0000
    CT2  CT2  CT1   58.350    113.50   11.16   2.56100
    CT2  CT2  CT2   58.350    113.60   11.16   2.56100
    CT2  CT3  CT1   58.350    113.50   11.16   2.56100
    CT2  NC2  C      62.300   120.0000
    CT2  NH1  C      50.000   120.0000
    CT2  OS   CD    40.000    109.60   30.00   2.26510
    CT3  CA   CA     45.800   122.3000
    CT3  CPH1 CPH1   45.800   130.0000
    CT3  CT1  C      52.000   108.0000
    CT3  CT1  CC     52.000   108.0000
    CT3  CT1  CT1   53.350    108.50    8.00   2.56100
    CT3  CT1  CT2   53.350    114.00    8.00   2.56100
    CT3  CT1  CT3   53.350    114.00    8.00   2.56100
    CT3  CT2  CA     51.800   107.5000
    CT3  CT2  CPH1   58.350   113.0000
    CT3  CT2  CT1   58.350    113.50   11.16   2.56100
    CT3  CT2  CT2   58.000    115.00    8.00   2.56100
    CT3  CT2  CT3   53.350    114.00    8.00   2.56100
    CT3  NC2  C      62.300   120.0000
    CT3  NH1  C      50.000   120.0000
    CT3  OS   CD    40.000    109.60   30.00   2.26510
    CT3  S    CT2    34.000    95.0000
    H    NH1  C      34.000   123.0000
    H    NH1  CT1    35.000   117.0000
    H    NH1  CT2    35.000   117.0000
    H    NH1  CT3    35.000   117.0000
    H    NH2  CC     50.000   120.0000
    H    NH2  H      23.000   120.0000
    H    NR1  CPH1  30.000    125.50   20.00   2.15000
    H    NR1  CPH2  30.000    127.00   20.00   2.14000
    H    NR3  CPH1  25.000    126.00   15.00   2.13000
    H    NR3  CPH2  25.000    126.00   15.00   2.09000
    H    OH1  CA     65.000   108.0000
    H    OH1  CD     55.000   115.0000
    H    OH1  CT1    57.500   106.0000
    H    OH1  CT2    57.500   106.0000
    H    OH1  CT3    57.500   106.0000
    HA2  CP2  CP1   33.430    110.10   22.53   2.17900
    HA2  CP2  CP2   26.500    110.10   22.53   2.17900
    HA2  CP2  CP3   26.500    110.10   22.53   2.17900
    HA2  CP2  HA2   35.500    109.00    5.40   1.80200
    HA2  CP3  CP2   26.500    110.10   22.53   2.17900
    HA2  CP3  HA2   35.500    109.00    5.40   1.80200
    HA2  CS   CT3   34.600    110.10   22.53   2.17900
    HA2  CS   HA2   35.500    108.40   14.00   1.77500
    HA3  CS   HA3   35.500    108.40   14.00   1.77500
    HA1  CT1  C     33.000    109.50   30.00   2.16300
    HA1  CT1  CD    33.000    109.50   30.00   2.16300
    HA1  CT1  CT1   34.500    110.10   22.53   2.17900
    HA1  CT1  CT2   34.500    110.10   22.53   2.17900
    HA1  CT1  CT3   34.500    110.10   22.53   2.17900
    HA1  CT1  HA1   35.500    109.00    5.40   1.80200
    HA2  CT2  C     33.000    109.50   30.00   2.16300
    HA2  CT2  CA     49.300   107.5000
    HA2  CT2  CC    33.000    109.50   30.00   2.16300
    HA2  CT2  CD    33.000    109.50   30.00   2.16300
    HA2  CT2  CE1    45.00    111.50  	
    HA2  CT2  CPH1   33.430   109.5000
    HA2  CT2  CT1   26.500    110.10   22.53   2.17900
    HA2  CT2  CT2   26.500    110.10   22.53   2.17900
    HA2  CT2  CT3   34.600    110.10   22.53   2.17900
    HA2  CT2  HA2   35.500    109.00    5.40   1.80200
    HA3  CT3  C     33.000    109.50   30.00   2.16300
    HA3  CT3  CA     49.300   107.5000
    HA3  CT3  CC    33.000    109.50   30.00   2.16300
    HA3  CT3  CD    33.000    109.50   30.00   2.16300
    HA3  CT3  CE1    42.00    111.50  	
    HA3  CT3  CPH1   33.430   109.5000
    HA3  CT3  CS    34.600    110.10   22.53   2.17900
    HA3  CT3  CT1   33.430    110.10   22.53   2.17900
    HA3  CT3  CT2   34.600    110.10   22.53   2.17900
    HA3  CT3  CT3   37.500    110.10   22.53   2.17900
    HA3  CT3  HA3   35.500    108.40    5.40   1.80200
    HE1  CE1  CE1    52.00    119.50  	
    HE1  CE1  CE2    42.00    118.00  	
    HE1  CE1  CT2    40.00    116.00  	
    HE1  CE1  CT3    22.00    117.00  	
    HE2  CE2  CE1    45.00    120.50  	
    HE2  CE2  CE2    55.50    120.50  	
    HE2  CE2  HE2    19.00    119.00  	
    HB1  CP1  C      50.000   112.0000
    HB1  CP1  CC     50.000   112.0000
    HB1  CP1  CD     50.000   112.0000
    HB1  CP1  CP2    35.000   118.0000
    HB1  CT1  C      50.000   109.5000
    HB1  CT1  CC     50.000   109.5000
    HB1  CT1  CD     50.000   109.5000
    HB1  CT1  CT1    35.000   111.0000
    HB1  CT1  CT2    35.000   111.0000
    HB1  CT1  CT3    35.000   111.0000
    HB2  CT2  C      50.000   109.5000
    HB2  CT2  CC     50.000   109.5000
    HB2  CT2  CD     50.000   109.5000
    HB2  CT2  HB2    36.000   115.0000
    HC   NC2  C      49.000   120.0000
    HC   NC2  CT2    40.400   120.0000
    HC   NC2  CT3    40.400   120.0000
    HC   NC2  HC     25.000   120.0000
    HC   NH2  CT2    50.000   111.0000
    HC   NH2  CT3    50.000   111.0000
    HC   NH2  HC     39.000   106.5000
    HC   NH3  CT1   30.000    109.50   20.00   2.07400
    HC   NH3  CT2   30.000    109.50   20.00   2.07400
    HC   NH3  CT3   30.000    109.50   20.00   2.07400
    HC   NH3  HC     44.000   109.5000
    HC   NP   CP1   33.000    109.50    4.00   2.05600
    HC   NP   CP3   33.000    109.50    4.00   2.05600
    HC   NP   HC     51.000   107.5000
    HP   CA   CA    30.000    120.00   22.00   2.15250
    HR1  CPH1 CPH1  22.000    130.00   15.00   2.21500
    HR3  CPH1 CPH1  25.000    130.00   20.00   2.20000
    HS   S    CT2    38.800    95.0000
    HS   S    CT3    43.000    95.0000
    N    C    CP1    20.000   112.5000
    N    C    CT1    20.000   112.5000
    N    C    CT2    20.000   112.5000
    N    C    CT3    20.000   112.5000
    N    CP1  C      50.000   108.2000
    N    CP1  CC     50.000   108.2000
    N    CP1  CD     50.000   108.2000
    N    CP1  CP2    70.000   110.8000
    N    CP1  HB1    48.000   112.0000
    N    CP3  CP2    70.000   110.5000
    N    CP3  HA2    48.000   108.0000
    NC2  C    NC2   40.000    120.00   70.00   2.31
    NC2  CT2  CT2    67.700   107.5000
    NC2  CT2  HA2    56.500   107.5000
    NC2  CT3  HA3    56.5000   107.5000
    NH1  C    CP1    80.000   116.5000
    NH1  C    CT1    80.000   116.5000
    NH1  C    CT2    80.000   116.5000
    NH1  C    CT3    80.000   116.5000
    NH1  CT1  C      50.000   107.0000
    NH1  CT1  CC     50.000   107.0000
    NH1  CT1  CD     50.000   107.0000
    NH1  CT1  CT1    70.000   113.5000
    NH1  CT1  CT2    70.000   113.5000
    NH1  CT1  CT3    70.000   113.5000
    NH1  CT1  HB1    48.000   108.0000
    NH1  CT2  C      50.000   107.0000
    NH1  CT2  CC     50.000   107.0000
    NH1  CT2  CD     50.000   107.0000
    NH1  CT2  CT2    70.000   113.5000
    NH1  CT2  HA2    51.500   109.5000
    NH1  CT2  HB2    48.000   108.0000
    NH1  CT3  HA3    51.500   109.5000
    NH2  CC   CP1    80.000   112.5000
    NH2  CC   CT1   50.000    116.50   50.00   2.45000
    NH2  CC   CT2   50.000    116.50   50.00   2.45000
    NH2  CC   CT3   50.000    116.50   50.00   2.45000
    NH2  CC   HA1   44.000    111.00   50.00   1.98000
    NH2  CT2  HB2   38.000    109.50   50.00   2.14000
    NH2  CT2  CD    52.000   108.0000
    NH2  CT2  CT2    67.700   110.0000
    NH2  CT2  HA2   38.000    109.50   50.00   2.14000
    NH2  CT3  HA3   38.000    109.50   50.00   2.14000
    NH3  CT1  C      43.700   110.0000
    NH3  CT1  CC     43.700   110.0000
    NH3  CT1  CT1    67.700   110.0000
    NH3  CT1  CT2    67.700   110.0000
    NH3  CT1  CT3    67.700   110.0000
    NH3  CT1  HB1    51.500   107.5000
    NH3  CT2  C      43.700   110.0000
    NH3  CT2  CC     43.700   110.0000
    NH3  CT2  CD     43.700   110.0000
    NH3  CT2  CT2    67.700   110.0000
    NH3  CT2  CT3    67.700   110.0000
    NH3  CT2  HA2   45.000    107.50   35.00   2.10100
    NH3  CT2  HB2    51.500   107.5000
    NH3  CT3  HA3   45.000    107.50   35.00   2.10100
    NP   CP1  C      50.000   106.0000
    NP   CP1  CC     50.000   106.0000
    NP   CP1  CD     50.000   106.0000
    NP   CP1  CP2    70.000   108.5000
    NP   CP1  HB1    51.500   107.5000
    NP   CP3  CP2    70.000   108.5000
    NP   CP3  HA2    51.500   109.1500
    NR1  CPH1 CPH1  130.000   106.0000
    NR1  CPH1 CT2    45.800   124.0000
    NR1  CPH1 CT3    45.800   124.0000
    NR1  CPH1 HR3   25.000    124.00   20.00   2.14000
    NR1  CPH2 HR1   25.000    122.50   20.00   2.14000
    NR2  CPH1 CPH1  130.000   110.0000
    NR2  CPH1 CT2    45.800   120.0000
    NR2  CPH1 HR3   25.000    120.00   20.00   2.14000
    NR2  CPH2 HR1   25.000    125.00   20.00   2.12000
    NR2  CPH2 NR1   130.000   112.5000
    NR3  CPH1 CPH1  145.000   108.0000
    NR3  CPH1 CT2    45.800   122.0000
    NR3  CPH1 HR1   22.000    122.00   15.00   2.18000
    NR3  CPH2 HR2   32.000    126.00   25.00   2.14000
    NR3  CPH2 NR3   145.000   108.0000
    O    C    CP1    80.000   118.0000
    O    C    CT1    80.000   121.0000
    O    C    CT2    80.000   121.0000
    O    C    CT3    80.000   121.0000
    O    C    H      50.000   121.7000
    O    C    N      80.000   122.5000
    O    C    NH1    80.000   122.5000
    O    CC   CP1    80.000   118.0000
    O    CC   CT1   15.000    121.00   50.00   2.44000
    O    CC   CT2   15.000    121.00   50.00   2.44000
    O    CC   CT3   15.000    121.00   50.00   2.44000
    O    CC   HA1    44.000   122.0000
    O    CC   NH2   75.000    122.50   50.00   2.37000
    OB   CD   CP1   70.000    125.00   20.00   2.44200
    OB   CD   CT1   70.000    125.00   20.00   2.44200
    OB   CD   CT2   70.000    125.00   20.00   2.44200
    OB   CD   CT3   70.000    125.00   20.00   2.44200
    OC   CA   CA     40.000   120.0000
    OC   CC   CP1   40.000    118.00   50.00   2.38800
    OC   CC   CT1   40.000    118.00   50.00   2.38800
    OC   CC   CT2   40.000    118.00   50.00   2.38800
    OC   CC   CT3   40.000    118.00   50.00   2.38800
    OC   CC   OC   100.000    124.00   70.00   2.22500
    OC   CT2  CT3    65.000   122.0000
    OC   CT2  HA2    65.000   118.3000
    OC   CT3  HA3    65.000   118.3000
    OH1  CA   CA     45.200   120.0000
    OH1  CD   CT2    55.000   110.5000
    OH1  CD   CT3    55.000   110.5000
    OH1  CD   OB    50.000    123.00  210.00   2.26200
    OH1  CT1  CT1    75.700   110.1000
    OH1  CT1  CT3    75.700   110.1000
    OH1  CT1  HA1    45.900   108.8900
    OH1  CT2  CT1    75.700   110.1000
    OH1  CT2  CT2    75.700   110.1000
    OH1  CT2  CT3    75.700   110.1000
    OH1  CT2  HA2    45.900   108.8900
    OH1  CT3  HA3    45.900   108.8900
    OS   CD   CP1   55.000    109.00   20.00   2.32600
    OS   CD   CT1   55.000    109.00   20.00   2.32600
    OS   CD   CT2   55.000    109.00   20.00   2.32600
    OS   CD   CT3   55.000    109.00   20.00   2.32600
    OS   CD   OB    90.000    125.90  160.00   2.25760
    OS   CT2  HA2    60.000   109.5000
    OS   CT3  HA3    60.000   109.5000
    S    CT2  CT1    58.000   112.5000
    S    CT2  CT2    58.000   114.5000
    S    CT2  CT3    58.000   114.5000
    S    CT2  HA2    46.100   111.3000
    S    CT3  HA3    46.100   111.3000
    SM   CT2  CT1    58.000   112.5000
    SM   CT2  CT3    58.000   112.5000
    SM   CT2  HA2    38.000   111.0000
    SM   CT3  HA3    38.000   111.0000
    SM   SM   CT2    72.500   103.3000
    SM   SM   CT3    72.500   103.3000
    SS   CS   CT3    55.000   118.0000
    SS   CS   HA2    40.000   112.3000
    SS   CS   HA3    40.000   112.3000
    O    CD   HR1    75.000   121.0000
    NH1  CT1  CT2A   70.000   113.5000
    HB1  CT1  CT2A   35.000   111.0000
    CT2A CT1  C      52.000   108.0000
    CT1  CT2A HA2    26.500   110.1000  22.53   2.17900
    CT1  CT2A CT2    58.350   113.5000  11.16   2.56100
    HA2  CT2A HA2    35.500   109.0000   5.40   1.80200
    HA2  CT2A CT2    26.500   110.1000  22.53   2.17900
    CT2A CT2  HA2    26.500   110.1000  22.53   2.17900
    CT2A CT2  CC     52.000   108.0000
    CT1  CT2A CPH1   58.350   113.0000
    HA2  CT2A CPH1   33.430   109.5000
    CT2A CPH1 CPH1   45.800   130.0000
    CT2A CPH1 NR3    45.800   122.0000
    CT1  CT2A CC     52.000   108.0000
    HA2  CT2A CC     33.000   109.5000  30.00   2.16300
    OC   CC   CT2A   40.000   118.0000  50.00   2.38800
    NH3  CT1  CT2A   67.700   110.0000
    CT2A CT1  CD     52.000   108.0000
    NH2  CT1  CS     67.700   110.0000
    CS   CT1  C      52.000   108.0000
    CS   CT1  CC     52.000   108.0000
    CS   CT1  CD     52.000   108.0000
    HB1  CT1  CS     35.000   111.0000
    NH1  CT1  CS     70.000   113.5000
    NH3  CT1  CS     67.700   110.0000
    SS   CS   CT1    55.000   118.0000
    HA2  CS   CT1    34.600   110.10    22.53   2.17900
    OC   CT2  CT1    65.000   122.0000
    HT   OT   HT     55.0      104.52  
    HT   OTMG   HT     55.000   104.5200
    HT   OTMG   MG      0.000   104.5200
    OTMG MG    OTMG     0.000   104.5200
    HT   OTCA   HT     55.000   104.5200
    HT   OTCA   CAL     0.000   104.5200
    OTCA CAL   OTCA     0.000   104.5200", "\n")

dihedral_para_info = split("X    CTL1 OHL  X        0.14    3     0.00 
    X    CTL2 OHL  X        0.14    3     0.00
    X    CTL3 OHL  X        0.14    3     0.00
    OCL  CCL  CTL1 NH3L     3.20    2   180.00
    OBL  CL   CTL2 HAL2     0.00    6   180.00
    OBL  CL   CTL3 HAL3     0.00    6   180.00
    OSL  CL   CTL2 HAL2     0.00    6   180.00
    OSL  CL   CTL3 HAL3     0.00    6   180.00
    OSLP CL   CTL2 HAL2     0.00    6   180.00
    OSLP CL   CTL3 HAL3     0.00    6   180.00
    OBL  CL   OSL  CTL1     0.965   1   180.00
    OBL  CL   OSL  CTL1     3.85    2   180.00
    OBL  CL   OSL  CTL2     0.965   1   180.00
    OBL  CL   OSL  CTL2     3.85    2   180.00
    OBL  CL   OSL  CTL3     0.965   1   180.00
    OBL  CL   OSL  CTL3     3.85    2   180.00
    X    CL   OSL  X        2.05    2   180.00
    X    CTL2 CL   X        0.05    6   180.00
    X    CTL3 CL   X        0.05    6   180.00
    X    CL   OHL  X        2.05    2   180.00
    X    CTL1 CCL  X        0.05    6   180.00
    HAL2 CTL2 CL   OHL      0.00    6   180.00
    HAL3 CTL3 CL   OHL      0.00    6   180.00
    PL   OSLP CTL2 CTL1     0.407   2     0.00
    PL   OSLP CTL2 CTL1     0.241   1   180.00
    PL   OSLP CTL2 CTL2     0.407   2     0.00
    PL   OSLP CTL2 CTL2     0.241   1   180.00
    OSL  PL   OSL  CTL1     1.20    1   180.00
    OSL  PL   OSL  CTL1     0.10    2   180.00
    OSL  PL   OSL  CTL1     0.10    3   180.00
    OSLP PL   OSLP CTL1     1.20    1   180.00
    OSLP PL   OSLP CTL1     0.10    2   180.00
    OSLP PL   OSLP CTL1     0.10    3   180.00
    OSLP PL   OSLP CTL2     1.20    1   180.00
    OSLP PL   OSLP CTL2     0.10    2   180.00
    OSLP PL   OSLP CTL2     0.10    3   180.00
    O2L  PL   OSLP CTL2     0.10    3     0.00
    O2L  PL   OSL  CTL2     0.10    3     0.00
    OSLP PL   OSLP CTL3     1.20    1   180.00
    OSLP PL   OSLP CTL3     0.10    2   180.00
    OSLP PL   OSLP CTL3     0.10    3   180.00
    O2L  PL   OSLP CTL1     0.10    3     0.00
    O2L  PL   OSL  CTL1     0.10    3     0.00
    O2L  PL   OSLP CTL3     0.10    3     0.00
    O2L  PL   OSL  CTL3     0.10    3     0.00
    OHL  PL   OSL  CTL1     0.95    2     0.00
    OHL  PL   OSL  CTL1     0.50    3     0.00
    OHL  PL   OSL  CTL2     0.95    2     0.00
    OHL  PL   OSL  CTL2     0.50    3     0.00
    OHL  PL   OSL  CTL3     0.95    2     0.00
    OHL  PL   OSL  CTL3     0.50    3     0.00
    OHL  PL   OSLP CTL2     0.95    2     0.00
    OHL  PL   OSLP CTL2     0.50    3     0.00
    OHL  PL   OSLP CTL3     0.95    2     0.00
    OHL  PL   OSLP CTL3     0.50    3     0.00
    X    OHL  PL   X        0.30    3     0.00
    X    CTL1 OSL  X        0.00    3     0.00
    X    CTL2 OSL  X        0.00    3     0.00
    X    CTL3 OSL  X        0.00    3     0.00
    X    CTL1 OSLP X        0.00    3     0.00
    X    CTL2 OSLP X        0.00    3     0.00
    X    CTL3 OSLP X        0.00    3     0.00
    CTL1 CTL2 CL   OSL     -0.15    1   180.00
    CTL1 CTL2 CL   OSL      0.53    2   180.00
    CTL2 CTL2 CL   OSL      0.000   6     0.00
    CTL2 CTL2 CL   OSL      0.030   3   180.00
    CTL2 CTL2 CL   OSL      0.432   2   180.00
    CTL2 CTL2 CL   OSL      0.332   1     0.00
    CTL3 CTL2 CL   OSL      0.000   6     0.00
    CTL3 CTL2 CL   OSL      0.030   3   180.00
    CTL3 CTL2 CL   OSL      0.432   2   180.00
    CTL3 CTL2 CL   OSL      0.332   1     0.00
    CTL3 CTL2 CTL2 CL       0.000   5   180.00
    CTL3 CTL2 CTL2 CL       0.317   3   180.00
    CTL3 CTL2 CTL2 CL       0.557   2     0.00
    CTL3 CTL2 CTL2 CL       0.753   1     0.00
    CTL2 CTL2 CTL2 CL       0.000   5   180.00
    CTL2 CTL2 CTL2 CL       0.317   3   180.00
    CTL2 CTL2 CTL2 CL       0.557   2     0.00
    CTL2 CTL2 CTL2 CL       0.753   1     0.00
    OSL  CTL2 CTL1 OSL     -0.429   4    60.00
    OSL  CTL2 CTL1 OSL      0.614   3     0.00
    OSL  CTL2 CTL1 OSL     -0.115   2    60.00
    OSL  CTL2 CTL1 OSL      0.703   1   180.00
    OSLP CTL2 CTL1 OSL      0.000   4     0.00
    OSLP CTL2 CTL1 OSL      0.607   3   180.00
    OSLP CTL2 CTL1 OSL      0.254   2    60.00
    OSLP CTL2 CTL1 OSL      2.016   1   180.00
    OSLP CTL2 CTL2 OSL      0.000   4     0.00
    OSLP CTL2 CTL2 OSL      0.607   3   180.00
    OSLP CTL2 CTL2 OSL      0.254   2    60.00
    OSLP CTL2 CTL2 OSL      2.016   1   180.00
    CTL3 CTL1 CTL2 OSL      0.000   3     0.00
    CTL2 CTL1 CTL2 OSL      0.000   3     0.00
    CTL3 CTL2 CTL2 OSL      0.000   3     0.00
    CTL2 CTL2 CTL2 OSL      0.000   3     0.00
    CL   OSL  CTL1 CTL2     0.000   4     0.00
    CL   OSL  CTL1 CTL2     0.150   3   180.00
    CL   OSL  CTL1 CTL2     1.453   2   180.00
    CL   OSL  CTL1 CTL2     0.837   1   180.00
    CL   OSL  CTL1 CTL3     0.000   4     0.00
    CL   OSL  CTL1 CTL3     0.150   3   180.00
    CL   OSL  CTL1 CTL3     1.453   2   180.00
    CL   OSL  CTL1 CTL3     0.837   1   180.00
    CL   OSL  CTL2 CTL1     0.267   3   180.00
    CL   OSL  CTL2 CTL1     0.173   2     0.00
    CL   OSL  CTL2 CTL1     0.781   1   180.00
    X    CTL2 NTL  X        0.26    3     0.00
    X    CTL5 NTL  X        0.23    3     0.00
    X    CTL1 NH3L X        0.10    3     0.00
    X    CTL2 NH3L X        0.10    3     0.00
    NH3L CTL2 CTL2 OHL      0.7     1   180.00
    NH3L CTL2 CTL2 OSLP     0.7     1   180.00
    NTL  CTL2 CTL2 OHL      4.3     1   180.00
    NTL  CTL2 CTL2 OHL     -0.4     3   180.00
    NTL  CTL2 CTL2 OSLP     3.3     1   180.00
    NTL  CTL2 CTL2 OSLP    -0.4     3   180.00
    X    CTL1 CTL1 X        0.200   3     0.00
    X    CTL1 CTL2 X        0.200   3     0.00
    X    CTL1 CTL3 X        0.200   3     0.00
    X    CTL2 CTL2 X        0.1900  3     0.00
    X    CTL2 CTL3 X        0.1600  3     0.00
    X    CTL3 CTL3 X        0.1525  3     0.00
    CTL3 CTL2 CTL2 CTL3     0.060   2     0.00
    CTL3 CTL2 CTL2 CTL3     0.035   5     0.00
    CTL2 CTL2 CTL2 CTL3     0.162   2     0.00
    CTL2 CTL2 CTL2 CTL3     0.047   3   180.00
    CTL2 CTL2 CTL2 CTL3     0.105   4     0.00
    CTL2 CTL2 CTL2 CTL3     0.177   5     0.00
    CTL2 CTL2 CTL2 CTL2     0.101   2     0.00
    CTL2 CTL2 CTL2 CTL2     0.142   3   180.00
    CTL2 CTL2 CTL2 CTL2     0.074   4     0.00
    CTL2 CTL2 CTL2 CTL2     0.097   5     0.00
    HAL3 CTL3 OSL  SL       0.00    3     0.00
    CTL2 OSL  SL   O2L      0.00    3     0.00
    CTL3 OSL  SL   O2L      0.00    3     0.00
    HEL1 CEL1 CEL1 HEL1     1.0000  2   180.00
    CTL3 CEL1 CEL1 HEL1     1.0000  2   180.00
    X    CEL1 CEL1 X        0.4500  1   180.00
    X    CEL1 CEL1 X        8.5000  2   180.00
    X    CEL2 CEL2 X        4.9000  2   180.00
    CTL2 CEL1 CEL2 HEL2     5.2000  2   180.00
    CTL3 CEL1 CEL2 HEL2     5.2000  2   180.00
    HEL1 CEL1 CEL2 HEL2     5.2000  2   180.00
    CEL1 CEL1 CTL2 HAL2     0.3000  3   180.00
    CEL1 CEL1 CTL3 HAL3     0.3000  3   180.00
    CEL1 CEL1 CTL2 CTL3     0.9100  1    180.0
    CEL1 CEL1 CTL2 CTL3     0.1800  2    180.0
    CEL1 CEL1 CTL2 CTL3     0.1700  3    180.0
    CEL1 CEL1 CTL2 CTL2     0.9100  1    180.0
    CEL1 CEL1 CTL2 CTL2     0.1800  2    180.0
    CEL1 CEL1 CTL2 CTL2     0.1700  3    180.0
    CEL1 CTL2 CTL2 CL       0.1400  1    180.0
    CEL1 CTL2 CTL2 CL       0.1700  2      0.0
    CEL1 CTL2 CTL2 CL       0.0500  3    180.0
    CEL1 CTL2 CTL2 CTL2     0.1400  1    180.0
    CEL1 CTL2 CTL2 CTL2     0.1700  2      0.0
    CEL1 CTL2 CTL2 CTL2     0.0500  3    180.0
    CEL1 CTL2 CTL2 CTL3     0.1400  1    180.0
    CEL1 CTL2 CTL2 CTL3     0.1700  2      0.0
    CEL1 CTL2 CTL2 CTL3     0.0500  3    180.0
    CEL2 CEL1 CTL2 CTL2     0.5000  1   180.00
    CEL2 CEL1 CTL2 CTL2     1.3000  3   180.00
    CEL2 CEL1 CTL2 CTL3     0.5000  1   180.00
    CEL2 CEL1 CTL2 CTL3     1.3000  3   180.00
    CEL2 CEL1 CTL2 HAL2     0.1200  3     0.00
    CEL2 CEL1 CTL3 HAL3     0.0500  3   180.00
    HEL1 CEL1 CTL2 CTL2     0.1200  3     0.00
    HEL1 CEL1 CTL2 CTL3     0.1200  3     0.00
    HEL1 CEL1 CTL2 HAL2     0.0000  3     0.00
    HEL1 CEL1 CTL3 HAL3     0.0000  3     0.00
    CEL2 CEL1 CTL2 CEL1     1.200   1   180.00
    CEL2 CEL1 CTL2 CEL1     0.400   2   180.00
    CEL2 CEL1 CTL2 CEL1     1.300   3   180.00
    CEL1 CTL2 CEL1 HEL1     0.000   2     0.00
    CEL1 CTL2 CEL1 HEL1     0.000   3     0.00
    CEL1 CEL1 CTL2 CEL1     0.850   1   180.00
    CEL1 CEL1 CTL2 CEL1     0.300   2   180.00
    CEL1 CEL1 CTL2 CEL1     0.260   3     0.00
    CEL1 CEL1 CTL2 CEL1     0.096   4     0.00
    X    CN8  ON2  X       -0.10    3     0.0
    X    CN7  CN8  X        0.20    3     0.0
    X    CN8  NN6  X        0.10    3     0.0
    CN7  ON6  CN8B HN8      0.195   1     0.0
    ON6  CN8B CN8  HN8      0.195   1     0.0
    HN7  CN7  ON6  CN8B     0.195   3     0.0
    CN8B CN8  CN7  HN7      0.195   3     0.0
    HN8  CN8B CN8  HN8      0.195   3     0.0
    HN8  CN8B CN8  CN7      0.195   3     0.0
    CN8B CN7  ON6  CN8B     0.5     5     0.0
    CN8B CN7  ON6  CN8B     0.1     3   180.0
    CN8B CN7  ON6  CN8B     0.5     1     0.0
    CN8B CN8  CN7  ON5      0.4     5     0.0
    CN8B CN8  CN7  ON5      0.4     3     0.0
    CN8B CN8  CN7  ON5      0.7     2     0.0
    CN8B CN8  CN7  ON5      0.5     1   180.0
    CN8B CN8  CN7  ON2      0.4     5     0.0
    CN8B CN8  CN7  ON2      0.4     3     0.0
    CN8B CN8  CN7  ON2      0.7     2     0.0
    CN8B CN8  CN7  ON2      0.5     1   180.0
    CN7  ON6  CN8B CN8      0.6     6   180.0
    CN7  ON6  CN8B CN8      0.6     3     0.0
    ON6  CN8B CN8  CN7      0.7     5   180.0
    ON6  CN8B CN8  CN7      0.4     4     0.0
    ON6  CN8B CN8  CN7      0.4     3   180.0
    CN7  CN7  CN8  CN8B     0.5     4     0.0
    CN7  CN7  CN8  CN8B     0.1     3     0.0
    CN8B ON6  CN7  CN7      0.5     3     0.0
    ON2  P2   ON2  CN7      0.90    1   180.0 
    ON2  P2   ON2  CN7      0.40    2   180.0 
    ON2  P2   ON2  CN7      0.20    3   180.0 
    ON2  P    ON2  CN7      1.20    1   180.0 
    ON2  P    ON2  CN7      0.10    2   180.0 
    ON2  P    ON2  CN7      0.10    3   180.0 
    ON2  P    ON2  CN7      0.00    6     0.0 
    ON2  P    ON2  CN8      1.20    1   180.0 
    ON2  P    ON2  CN8      0.10    2   180.0 
    ON2  P    ON2  CN8      0.10    3   180.0 
    ON2  P    ON2  CN8      0.00    6     0.0 
    ON2  P2   ON2  CN8      1.20    1   180.0 
    ON2  P2   ON2  CN8      0.10    2   180.0 
    ON2  P2   ON2  CN8      0.10    3   180.0 
    ON2  P2   ON2  CN8      0.00    6     0.0 
    ON2  P2   ON2  CN8B     1.20    1   180.0 
    ON2  P2   ON2  CN8B     0.10    2   180.0 
    ON2  P2   ON2  CN8B     0.10    3   180.0 
    ON2  P2   ON2  CN8B     0.00    6     0.0 
    ON2  P    ON2  CN8B     1.20    1   180.0 
    ON2  P    ON2  CN8B     0.10    2   180.0 
    ON2  P    ON2  CN8B     0.10    3   180.0 
    ON2  P    ON2  CN8B     0.00    6     0.0 
    ON2  P    ON2  CN9      1.20    1   180.0 
    ON2  P    ON2  CN9      0.10    2   180.0 
    ON2  P    ON2  CN9      0.10    3   180.0 
    ON2  P    ON2  CN9      0.00    6     0.0 
    ON2  P2   ON2  CN9      1.20    1   180.0 
    ON2  P2   ON2  CN9      0.10    2   180.0 
    ON2  P2   ON2  CN9      0.10    3   180.0 
    ON2  P2   ON2  CN9      0.00    6     0.0 
    ON3  P    ON2  CN7      0.10    3     0.0 
    ON3  P2   ON2  CN7      0.10    3     0.0 
    ON3  P    ON2  CN7B     0.10    3     0.0 
    ON3  P    ON2  CN8      0.10    3     0.0 
    ON3  P2   ON2  CN8      0.10    3     0.0 
    ON3  P    ON2  CN8B     0.10    3     0.0 
    ON3  P2   ON2  CN8B     0.10    3     0.0 
    ON3  P    ON2  CN9      0.10    3     0.0 
    ON3  P2   ON2  CN9      0.10    3     0.0 
    ON4  P    ON2  CN7      0.95    2     0.0 
    ON4  P    ON2  CN7      0.50    3     0.0 
    ON4  P2   ON2  CN7      0.95    2     0.0 
    ON4  P2   ON2  CN7      0.50    3     0.0 
    ON4  P    ON2  CN8      0.95    2     0.0 
    ON4  P    ON2  CN8      0.50    3     0.0 
    ON4  P    ON2  CN8B     0.95    2     0.0 
    ON4  P    ON2  CN8B     0.50    3     0.0 
    ON4  P2   ON2  CN8B     0.95    2     0.0 
    ON4  P2   ON2  CN8B     0.50    3     0.0 
    ON4  P    ON2  CN9      0.95    2     0.0 
    ON4  P    ON2  CN9      0.50    3     0.0 
    ON4  P2   ON2  CN9      0.95    2     0.0 
    ON4  P2   ON2  CN9      0.50    3     0.0 
    X    ON4  P    X        0.30    3     0.0 
    X    ON4  P2   X        0.30    3     0.0 
    P    ON2  CN7  HN7      0.000   3     0.0 
    P2   ON2  CN7  HN7      0.000   3     0.0 
    P    ON2  CN7B HN7      0.000   3     0.0 
    P    ON2  CN8B HN8      0.000   3     0.0 
    P2   ON2  CN8B HN8      0.000   3     0.0 
    P    ON2  CN8  HN8      0.000   3     0.0 
    P    ON2  CN9  HN9      0.000   3     0.0 
    P2   ON2  CN9  HN9      0.000   3     0.0 
    CN9  CN8  CN8  CN9      0.15    1     0.0
    CN9  CN8  CN8  CN8      0.15    1     0.0
    NN2B CN1T NN2U CN1      1.5     2   180.0 
    CN1T NN2U CN1  CN3      1.5     2   180.0 
    NN2U CN1  CN3  CN3      1.5     2   180.0 
    CN1  CN3  CN3  NN2B     6.0     2   180.0 
    CN3  CN3  NN2B CN1T     1.5     2   180.0 
    CN3  NN2B CN1T NN2U     1.5     2   180.0 
    HN3  CN3  CN3  HN3      3.0     2   180.0 
    HN3  CN3  CN1  ON1      6.0     2   180.0 
    ON1  CN1T NN2B HN2      0.0     2   180.0 
    ON1  CN1  NN2U HN2      0.0     2   180.0 
    ON1  CN1T NN2U HN2      0.0     2   180.0 
    HN2  NN2B CN3  HN3      1.5     2   180.0 
    NN2B CN1T NN2U HN2      3.8     2   180.0 
    CN3  CN1  NN2U HN2      3.8     2   180.0 
    CN3  CN3  NN2B HN2      1.6     2   180.0 
    NN2U CN1T NN2B HN2      1.6     2   180.0 
    CN1T NN2B CN3  CN3T     1.8     2   180.0 
    NN2U CN1  CN3T CN3      1.8     2   180.0 
    CN1  CN3T CN3  NN2B     3.0     2   180.0 
    NN2B CN1  CN3T CN9      5.6     2   180.0 
    NN2B CN3  CN3T CN9      5.6     2   180.0 
    CN1  CN3T CN9  HN9      0.46    3     0.0 
    CN3  CN3T CN9  HN9      0.46    3     0.0 
    CN3T CN1  NN2U HN2      4.8     2   180.0 
    CN3  NN2  CN1  NN3      0.6     2   180.0 
    NN2  CN1  NN3  CN2      0.6     2   180.0 
    CN1  NN3  CN2  CN3      6.0     2   180.0 
    NN3  CN2  CN3  CN3      0.6     2   180.0 
    CN2  CN3  CN3  NN2      6.0     2   180.0 
    CN3  CN3  NN2  CN1      0.6     2   180.0 
    NN3  CN2  NN1  HN1      1.0     2   180.0 
    CN3  CN2  NN1  HN1      1.0     2   180.0 
    NN1  CN2  NN3  CN1      2.0     2   180.0 
    NN1  CN2  CN3  CN3      2.0     2   180.0 
    NN1  CN2  CN3  HN3      2.0     2   180.0 
    ON1C CN1  NN2  HN2      3.0     2   180.0 
    ON1C CN1  NN3  CN2      1.6     2   180.0 
    ON1C CN1  NN2  CN3      1.6     2   180.0 
    NN3  CN2  CN3  HN3      3.4     2   180.0 
    NN2  CN3  CN3  HN3      3.4     2   180.0 
    CN2  CN3  CN3  HN3      4.6     2   180.0 
    CN1  NN2  CN3  HN3      4.6     2   180.0 
    X    CN2  NN3  X        2.0     2   180.0 
    CN2  NN3A CN4  NN3A     1.8     2   180.0 
    NN3A CN4  NN3A CN5      2.0     2   180.0 
    CN4  NN3A CN5  CN5      1.8     2   180.0 
    NN3A CN5  CN5  CN2      2.0     2   180.0 
    CN5  CN5  CN2  NN3A     1.8     2   180.0 
    CN5  CN2  NN3A CN4     10.0     2   180.0 
    CN5  CN5  NN4  CN4      6.0     2   180.0 
    CN5  NN4  CN4  NN2     14.0     2   180.0 
    NN4  CN4  NN2  CN5      6.0     2   180.0 
    CN4  NN2  CN5  CN5      6.0     2   180.0 
    NN2  CN5  CN5  NN4     14.0     2   180.0 
    CN2  NN3A CN4  HN3      8.5     2   180.0 
    CN5  NN3A CN4  HN3      8.5     2   180.0 
    CN5  NN4  CN4  HN3      5.2     2   180.0 
    CN5  NN2  CN4  HN3      5.2     2   180.0 
    CN5  CN5  NN2  HN2      1.2     2   180.0 
    NN4  CN4  NN2  HN2      1.2     2   180.0 
    HN2  NN2  CN4  HN3      0.0     2   180.0 
    CN4  NN3A CN2  NN1      4.0     2   180.0 
    CN5  CN5  CN2  NN1      4.0     2   180.0 
    NN4  CN5  CN2  NN1      0.0     2   180.0 
    CN5  CN2  NN1  HN1      0.5     2   180.0 
    NN3A CN2  NN1  HN1      0.5     2   180.0 
    NN3A CN5  CN5  NN4      7.0     2   180.0 
    CN2  CN5  CN5  NN2      7.0     2   180.0 
    NN3A CN2  CN5  NN4      2.0     2   180.0 
    CN2  CN5  NN4  CN4      2.0     2   180.0 
    CN4  NN3A CN5  NN2      2.0     2   180.0 
    NN3A CN5  NN2  CN4      2.0     2   180.0 
    CN1  NN2G CN2  NN3G     0.2     2   180.0 
    NN2G CN2  NN3G CN5      2.0     2   180.0 
    CN2  NN3G CN5  CN5G     0.2     2   180.0 
    NN3G CN5  CN5G CN1      2.0     2   180.0 
    CN5  CN5G CN1  NN2G     0.2     2   180.0 
    CN5G CN1  NN2G CN2      0.2     2   180.0 
    CN5  CN5G NN4  CN4      6.0     2   180.0 
    CN5G NN4  CN4  NN2B    16.0     2   180.0 
    NN4  CN4  NN2B CN5      6.0     2   180.0 
    CN4  NN2B CN5  CN5G     6.0     2   180.0 
    NN2B CN5  CN5G NN4     10.0     2   180.0 
    ON1  CN1  CN5G CN5     14.0     2   180.0 
    ON1  CN1  CN5G NN4      0.0     2   180.0 
    ON1  CN1  NN2G CN2     14.0     2   180.0 
    ON1  CN1  NN2G HN2      0.0     2   180.0 
    NN1  CN2  NN2G CN1      4.0     2   180.0 
    NN1  CN2  NN3G CN5      4.0     2   180.0 
    NN1  CN2  NN2G HN2      0.0     2   180.0 
    NN2G CN2  NN1  HN1      1.2     2   180.0 
    NN3G CN2  NN1  HN1      1.2     2   180.0 
    HN2  NN2G CN1  CN5G     3.6     2   180.0 
    HN2  NN2G CN2  NN3G     3.6     2   180.0 
    HN3  CN4  NN4  CN5G     5.6     2   180.0 
    HN3  CN4  NN2B CN5      5.6     2   180.0 
    HN3  CN4  NN2B HN2      0.0     2   180.0 
    HN2  NN2B CN5  CN5G     1.2     2   180.0 
    HN2  NN2B CN5  NN3G     1.2     2   180.0 
    HN2  NN2B CN4  NN4      1.2     2   180.0 
    NN3G CN5  CN5G NN4     10.0     2   180.0 
    CN1  CN5G CN5  NN2     10.0     2   180.0 
    NN2G CN1  CN5G NN4      2.0     2   180.0 
    CN1  CN5G NN4  CN4      2.0     2   180.0 
    CN2  NN3G CN5  NN2B     2.0     2   180.0 
    NN3G CN5  NN2B CN4      2.0     2   180.0 
    X    CN1  NN3  X        1.0     2   180.0 
    X    CN1  NN2  X        0.9     2   180.0 
    X    CN1T NN2B X        0.9     2   180.0 
    X    CN1  NN2G X        0.9     2   180.0 
    X    CN1  NN2U X        0.9     2   180.0 
    X    CN1T NN2U X        0.9     2   180.0 
    X    CN3  NN2  X        1.0     2   180.0 
    X    CN3  NN2B X        1.0     2   180.0 
    X    CN3  CN3  X        1.0     2   180.0 
    X    CN3  CN3T X        1.0     2   180.0 
    X    CN1  CN3  X        1.0     2   180.0 
    X    CN1  CN3T X        1.0     2   180.0 
    X    CN2  CN3  X        0.8     2   180.0 
    X    CN1  CN5G X        1.0     2   180.0 
    X    CN2  NN2G X        1.0     2   180.0 
    X    CN2  CN5  X        1.0     2   180.0 
    X    CN4  NN2  X        1.5     2   180.0 
    X    CN4  NN2B X        1.5     2   180.0 
    X    CN4  NN3A X        3.5     2   180.0 
    X    CN4  NN4  X        2.0     2   180.0 
    X    CN5  CN5  X        0.0     2   180.0 
    X    CN5G CN5  X        0.0     2   180.0 
    X    CN5  NN2  X        1.5     2   180.0 
    X    CN5  NN2B X        1.5     2   180.0 
    X    CN5  NN3A X        1.0     2   180.0 
    X    CN5  NN3G X        1.0     2   180.0 
    X    CN5  NN4  X        1.0     2   180.0 
    X    CN5G NN4  X        1.0     2   180.0 
    X    CN2  NN3A X        1.0     2   180.0 
    X    CN2  NN3G X        1.0     2   180.0 
    CN1  NN2  CN9  HN9      0.19    3     0.0
    CN3  NN2  CN9  HN9      0.00    3     0.0
    CN4  NN2  CN9  HN9      0.00    3     0.0
    CN5  NN2  CN9  HN9      0.19    3     0.0
    CN1  NN2B CN9  HN9      0.19    3     0.0
    CN1T NN2B CN9  HN9      0.19    3     0.0
    CN3  NN2B CN9  HN9      0.00    3     0.0
    CN4  NN2B CN9  HN9      0.00    3     0.0
    CN5  NN2B CN9  HN9      0.19    3     0.0
    CN4  NN2B CN8  HN8      0.00    3     0.0
    CN5  NN2B CN8  HN8      0.19    3     0.0
    CN4  NN2B CN8  CN9      0.00    3     0.0
    CN5  NN2B CN8  CN9      0.19    3     0.0
    X    CN8  CN8  X        0.15    3     0.0
    X    CN8  CN9  X        0.15    3     0.0
    HN7  CN7B CN7B ON2       0.195   3         0.0
    ON2  CN7B CN7B NN2       0.0     3         0.0
    CN7  CN7B ON6  CN7       0.6     6       180.0 
    CN7B CN7  CN7  CN7       0.4     6         0.0
    CN7B CN7  CN7  CN9       0.4     6         0.0
    CN7  CN7  CN7  ON6       0.6     6         0.0
    CN7  CN7  CN7B ON6       0.6     6         0.0
    ON2  CN7  CN7  CN7       0.8     6         0.0
    ON2  CN7  CN7  CN7       0.4     5         0.0
    ON2  CN7  CN7  CN7       2.0     3       180.0
    ON2  CN7B CN7  CN7       0.8     6         0.0
    ON2  CN7B CN7  CN7       0.4     5         0.0
    ON2  CN7B CN7  CN7       2.0     3       180.0
    ON5  CN7  CN7  CN7       0.8     6         0.0
    ON5  CN7  CN7  CN7       0.4     5         0.0
    ON5  CN7  CN7  CN7       2.0     3       180.0
    ON5  CN7  CN7  ON5       0.0     3         0.0
    ON5  CN7  CN7  ON2       0.0     3         0.0
    ON2  CN7  CN7B ON6       0.5     6         0.0
    ON2  CN7  CN7B ON6       0.3     5         0.0
    ON2  CN7  CN7B ON6       0.6     4       180.0
    ON2  CN7  CN7B ON6       0.2     3         0.0
    CN7  CN7  CN7  CN8B      0.5     4       180.0
    CN7B NN2  CN4  HN3       0.3     2       180.0
    CN7B NN2  CN5  CN5      11.0     2       180.0
    CN7B NN2  CN4  NN4      11.0     2       180.0
    CN7B NN2  CN4  NN3A     11.0     2       180.0
    ON6  CN7B NN2  CN5       1.1     1       180.0
    ON6  CN7B NN2  CN4       1.1     1         0.0
    ON6B CN7B NN2  CN5       1.1     1       180.0
    ON6B CN7B NN2  CN4       1.1     1         0.0
    CN8  CN7B NN2  CN5       0.3     3         0.0
    CN8  CN7B NN2  CN4       0.0     3       180.0
    CN7B CN7B NN2  CN5       0.3     3         0.0
    CN7B CN7B NN2  CN4       0.0     3       180.0
    HN7  CN7B NN2  CN5       0.0     3         0.0
    HN7  CN7B NN2  CN4       0.195   3         0.0
    CN7B NN2  CN3  HN3       0.3     2       180.0
    CN7B NN2  CN1  ON1C     11.0     2       180.0
    CN7B NN2  CN1  NN3      11.0     2       180.0
    CN7B NN2  CN3  CN3      11.0     2       180.0
    ON6  CN7B NN2  CN1       0.0     3         0.0
    ON6  CN7B NN2  CN3       1.0     1         0.0
    ON6B CN7B NN2  CN1       0.0     3         0.0
    ON6B CN7B NN2  CN3       1.0     1         0.0
    CN8  CN7B NN2  CN1       1.0     3         0.0
    CN8  CN7B NN2  CN3       0.0     3       180.0
    CN7B CN7B NN2  CN1       1.0     3         0.0
    CN7B CN7B NN2  CN3       0.0     3       180.0
    HN7  CN7B NN2  CN1       0.0     3         0.0
    HN7  CN7B NN2  CN3       0.195   3         0.0
    CN7B NN2B CN3  HN3       0.3     2       180.0
    CN7B NN2B CN1T ON1      11.0     2       180.0
    CN7B NN2B CN1T NN2U     11.0     2       180.0
    CN7B NN2B CN3  CN3T     11.0     2       180.0
    ON6  CN7B NN2B CN1       0.0     3         0.0
    ON6  CN7B NN2B CN1T      0.7     3         0.0
    ON6  CN7B NN2B CN1T      0.8     1       180.0
    ON6  CN7B NN2B CN3       0.9     1         0.0  
    ON6B CN7B NN2B CN1       0.0     3         0.0
    ON6B CN7B NN2B CN1T      0.7     3         0.0
    ON6B CN7B NN2B CN1T      0.8     1       180.0
    ON6B CN7B NN2B CN3       0.9     1         0.0 
    CN8  CN7B NN2B CN1T      0.2     3       180.0
    CN8  CN7B NN2B CN3       0.0     3       180.0
    CN7B CN7B NN2B CN1T      0.2     3       180.0
    CN7B CN7B NN2B CN3       0.0     3       180.0
    HN7  CN7B NN2B CN1T      0.0     3         0.0
    HN7  CN7B NN2B CN3       0.195   3         0.0
    CN7B NN2B CN4  HN3       0.3     2       180.0
    CN7B NN2B CN4  NN4      11.0     2       180.0
    CN7B NN2B CN5  CN5G     11.0     2       180.0
    CN7B NN2B CN5  NN3G     11.0     2       180.0
    ON6  CN7B NN2B CN5       0.2     3         0.0
    ON6  CN7B NN2B CN5       1.1     1       180.0
    ON6  CN7B NN2B CN4       1.4     1         0.0
    ON6B CN7B NN2B CN5       0.2     3         0.0
    ON6B CN7B NN2B CN5       1.1     1       180.0
    ON6B CN7B NN2B CN4       1.4     1         0.0 
    CN8  CN7B NN2B CN5       0.0     3         0.0
    CN8  CN7B NN2B CN4       0.0     3       180.0  
    CN7B CN7B NN2B CN5       0.0     3         0.0
    CN7B CN7B NN2B CN4       0.0     3       180.0
    HN7  CN7B NN2B CN5       0.0     3         0.0
    HN7  CN7B NN2B CN4       0.195   3         0.0
    CN7  ON6  CN7B NN2       0.0     3         0.0
    CN7  ON6  CN7B NN2B      0.0     3         0.0
    CN7  ON6B CN7B NN2       0.0     3         0.0
    CN7  ON6B CN7B NN2B      0.0     3         0.0
    CN7  CN8  CN7B NN2       0.0     3         0.0
    CN7  CN8  CN7B NN2B      0.0     3         0.0
    HN8  CN8  CN7B NN2       0.0     3         0.0
    HN8  CN8  CN7B NN2B      0.0     3         0.0 
    CN7  CN7B CN7B NN2       0.0     3         0.0
    CN7  CN7B CN7B NN2B      0.0     3         0.0
    HN7  CN7B CN7B NN2       0.0     3         0.0
    HN7  CN7B CN7B NN2B      0.0     3         0.0
    CN7  CN8B ON2  P        0.2     1       120.0
    CN7  CN8B ON2  P2       0.2     1       120.0
    CN7  CN8B ON5  HN5      1.3300  1         0.00
    CN7  CN8B ON5  HN5      0.1800  2         0.00
    CN7  CN8B ON5  HN5      0.3200  3         0.00
    HN8  CN8B ON5  HN5      0.0     3         0.0
    CN7  CN7  CN8B ON2      0.20    4       180.0
    CN7  CN7  CN8B ON2      0.80    3       180.0
    CN7  CN7  CN8B ON2      0.40    2         0.0
    CN7  CN7  CN8B ON2      2.50    1       180.0
    CN8  CN7  CN8B ON2      0.2     3       180.0
    CN7  CN7  CN8B ON5      0.20    4       180.0
    CN7  CN7  CN8B ON5      0.80    3       180.0
    CN7  CN7  CN8B ON5      0.40    2         0.0
    CN7  CN7  CN8B ON5      2.50    1       180.0
    ON6  CN7  CN8B ON2      3.4     1       180.0
    ON6B CN7  CN8B ON2      3.4     1       180.0
    ON6  CN7  CN8B ON5      3.4     1       180.0
    ON6B CN7  CN8B ON5      3.4     1       180.0
    HN8  CN8B CN7  CN7      0.195   3         0.0
    HN8  CN8B CN7  CN8      0.195   1         0.0
    HN8  CN8B CN7  ON6      0.195   1         0.0
    HN8  CN8B CN7  ON6B     0.195   1         0.0
    HN7  CN7  CN8B ON2      0.195	3         0.0
    HN7  CN7  CN8B ON5      0.195   3         0.0
    HN8  CN8  CN8  ON6      0.195   1         0.0
    CN9  CN7  CN7  CN8B     0.5     4       180.0
    HN7  CN7  CN9  HN9      0.195   3         0.0
    CN7  CN7  CN9  HN9      0.195   3         0.0
    ON6  CN7  CN9  HN9      0.195   3         0.0
    HN7  CN7  CN7  CN9      0.195   3         0.0
    ON2  CN7  CN7  CN9      0.2     4         0.0
    ON2  CN7  CN7  CN9      0.8     3       180.0
    CN8  CN7  CN7  CN9      0.5     4       180.0
    CN8  CN7  CN7  CN8B     0.5     4       180.0
    CN7B CN7  CN7  CN8B     0.2     4       180.0
    ON2  CN7  CN7  CN8B      0.2      4       0.0
    ON2  CN7  CN7  CN8B      0.8      3     180.0
    ON5  CN7  CN7  CN8B      0.2      4       0.0
    ON5  CN7  CN7  CN8B      0.8      3     180.0
    ON2  CN7  CN7  ON6       0.5     6        0.0
    ON2  CN7  CN7  ON6       0.3     5        0.0
    ON2  CN7  CN7  ON6       0.6     4      180.0
    ON2  CN7  CN7  ON6       0.2     3        0.0
    ON2  CN7  CN7  ON6B      0.4     6        0.0
    ON2  CN7  CN7  ON6B      0.0     5        0.0
    ON2  CN7  CN7  ON6B      0.0     4      180.0
    ON2  CN7  CN7  ON6B      1.6     3        0.0
    ON2  CN7B CN7B ON6B      0.4     6        0.0
    ON2  CN7B CN7B ON6B      0.0     5        0.0
    ON2  CN7B CN7B ON6B      0.0     4      180.0
    ON2  CN7B CN7B ON6B      1.6     3        0.0
    ON5  CN7  CN7  ON6       0.5     6        0.0
    ON5  CN7  CN7  ON6       0.3     5        0.0
    ON5  CN7  CN7  ON6       0.6     4      180.0
    ON5  CN7  CN7  ON6       0.2     3        0.0
    ON5  CN7  CN7  ON6B      0.4     6        0.0
    ON5  CN7  CN7  ON6B      0.0     5        0.0
    ON5  CN7  CN7  ON6B      0.0     4      180.0
    ON5  CN7  CN7  ON6B      1.6     3        0.0
    CN7B ON6  CN7  CN8B      0.8     3         0.0
    CN7B ON6B CN7  CN8B      2.0     3         0.0
    CN7B ON6B CN7  CN9       2.0     3         0.0
    ON2  CN7  CN8  CN7B      0.8     6         0.0
    ON2  CN7  CN8  CN7B      0.4     5         0.0
    ON2  CN7  CN8  CN7B      2.0     3       180.0
    ON2  CN7  CN7B CN7B      0.6     6         0.0
    ON2  CN7  CN7B CN7B      0.0     5         0.0
    ON2  CN7  CN7B CN7B      1.6     3       180.0
    ON5  CN7  CN8  CN7B      0.8     6         0.0
    ON5  CN7  CN8  CN7B      0.4     5         0.0
    ON5  CN7  CN8  CN7B      2.0     3       180.0
    ON5  CN7  CN7B CN7B      0.6     6         0.0
    ON5  CN7  CN7B CN7B      0.0     5         0.0
    ON5  CN7  CN7B CN7B      1.6     3       180.0
    ON2  CN7  CN8  HN8       0.195   3         0.0
    ON5  CN7  CN8  HN8       0.195   3       180.0
    ON2  CN7  CN7B HN7       0.195   3         0.0
    ON5  CN7  CN7B HN7       0.195   3       180.0
    HN7  CN7  CN7  ON2       0.195   3         0.0
    HN7  CN7  CN7  ON5       0.195   3         0.0
    CN7  CN7  ON2  P2        0.6     5         0.0
    CN7  CN7  ON2  P2        0.2     4         0.0
    CN7  CN7  ON2  P2        0.0     3       180.0
    CN7  CN7  ON2  P2        0.4     2         0.0
    CN7  CN7  ON2  P2        1.9     1       180.0
    CN7  CN7  ON2  P         0.6     5         0.0
    CN7  CN7  ON2  P         0.2     4         0.0
    CN7  CN7  ON2  P         0.0     3       180.0
    CN7  CN7  ON2  P         0.4     2         0.0
    CN7  CN7  ON2  P         1.9     1       180.0
    CN8  CN7  ON2  P         2.5     1       180.0
    CN8  CN7  ON2  P2        1.9     1       180.0
    CN7B CN7  ON2  P         2.5     1       180.0
    CN7B CN7B ON2  P         2.5     1       180.0
    CN7  CN7B ON2  P         2.5     1       180.0
    CN7  CN7  ON5  HN5       0.5     3         0.0
    CN7  CN7  ON5  HN5       0.3     2       180.0
    CN7  CN7  ON5  HN5        1.5     1        0.0
    CN8  CN7  ON5  HN5       0.5     3         0.0
    CN8  CN7  ON5  HN5       1.0     2       180.0
    CN8  CN7  ON5  HN5       0.3     1         0.0
    CN7B CN7  ON5  HN5       0.8     3         0.0
    CN7B CN7  ON5  HN5       0.5     1         0.0
    HN7  CN7  ON5  HN5       0.0     3         0.0
    HN7  CN7  CN8B HN8       0.195   3         0.0
    HN7  CN7  CN7  CN8B      0.195   3         0.0
    CN8  CN7B ON6  CN7       0.6     6       180.0
    CN8  CN7  CN7  ON6       1.0     4         0.0
    CN8  CN7  CN7  ON6       0.3     5       180.0
    CN8  CN7  CN7  ON6       0.3     6       180.0
    CN7B CN7B ON6B CN7       0.0     6         0.0
    CN7B CN7  CN7  ON6B      0.0     3         0.0
    CN7  CN8  CN7B ON6       0.6     6         0.0
    CN7  CN7B CN7B ON6B      0.4     6         0.0
    CN7B CN8  CN7  CN7       0.4      6        0.0
    CN7B CN7B CN7  CN7       0.0      6        0.0
    CN7  CN7  ON6  CN7B      0.6      6      180.0
    CN7  CN7  ON6B CN7B      0.0      6      180.0
    HN7  CN7  CN7  CN8       0.0      3        0.0
    HN7  CN7  CN8  CN7B      0.195    3       0.0
    HN7  CN7B CN8  CN7       0.195    3       0.0
    HN7  CN7  CN7  ON6       0.195    3     180.0
    HN8  CN8  CN7B ON6       0.195    3       0.0
    HN7  CN7  CN7  HN7       0.195    3       0.0
    HN7  CN7B CN8  HN8       0.195    3       0.0
    HN7  CN7  CN8  HN8       0.195    3       0.0
    HN8  CN8  CN7  CN7       0.195    3       0.0
    HN7  CN7  ON6  CN7B      0.195    3       0.0
    HN7  CN7B ON6  CN7       0.000    3       0.0
    HN7  CN7  CN7  ON6B      0.195    3     180.0
    HN9  CN9  CN7  ON6B      0.195    3     180.0
    HN8  CN8  CN7B ON6B      0.195    3       0.0
    HN7  CN7B ON6B CN7       0.000    3       0.0
    HN7  CN7  ON6B CN7B      0.195    3       0.0
    HN7  CN7  CN7B CN7B      0.195    3       0.0
    HN7  CN7B CN7B CN7       0.195    3       0.0
    HN7  CN7B CN7B ON6B      0.195    3       0.0
    NN2  CN7B CN7B ON5       0.000    3       0.0
    NN2B CN7B CN7B ON5       0.000    3       0.0
    ON5  CN7B CN7B HN7       0.000    3       0.0
    HN7  CN7B CN7B HN7       0.000    3       0.0
    CN7  CN7  CN7B ON5       0.000    3       0.0
    ON6B CN7B CN7B ON5       0.000    3       0.0
    ON5  CN7B CN7  ON2       0.000    3       0.0
    ON5  CN7  CN7B ON2       0.000    3       0.0
    ON5  CN7B CN7  ON5       0.000    3       0.0
    HN7  CN7B ON5  HN5       0.000    3       0.0
    HN5  ON5  CN7B CN7B      0.000    6     180.0
    HN5  ON5  CN7B CN7B      0.400    3       0.0
    HN5  ON5  CN7B CN7B      0.400    2       0.0
    HN5  ON5  CN7B CN7B      0.800    1       0.0
    HN5  ON5  CN7B CN7       0.200    3       0.0
    HN5  ON5  CN7B CN7       0.000    2     180.0
    HN5  ON5  CN7B CN7       2.000    1       0.0
    CN7B ON6  CN7  CN9       0.0     3         0.0
    HN7  CN7  CN7B ON5       0.195   3         0.0
    HN7  CN7B CN7  CN7       0.195   3         0.0
    HN7  CN7  CN7  CN7B      0.195   3         0.0
    HN7  CN7  CN7B HN7       0.195   3         0.0
    P3   ON2  P    ON2      0.03    2     0.0 
    P3   ON2  P    ON2      0.03    3     0.0 
    P    ON2  P3   ON2      0.03    2     0.0 
    P    ON2  P3   ON2      0.03    3     0.0 
    P3   ON2  P3   ON2      0.03    2     0.0 
    P3   ON2  P3   ON2      0.03    3     0.0 
    P    ON2  P    ON2      0.03    2     0.0 
    P    ON2  P    ON2      0.03    3     0.0 
    P    ON2  P    ON3      0.10    2     0.0 
    P    ON2  P    ON3      0.03    3     0.0 
    P    ON2  P3   ON3      0.10    2     0.0 
    P    ON2  P3   ON3      0.03    3     0.0 
    P3   ON2  P    ON3      0.10    2     0.0 
    P3   ON2  P    ON3      0.03    3     0.0 
    P3   ON2  P3   ON3      0.10    2     0.0 
    P3   ON2  P3   ON3      0.03    3     0.0 
    P    ON2  P4   ON2      0.03    2     0.0 
    P    ON2  P4   ON2      0.03    3     0.0 
    P4   ON2  P    ON2      0.03    2     0.0 
    P4   ON2  P    ON2      0.03    3     0.0 
    P    ON2  P4   ON3      0.10    2     0.0 
    P    ON2  P4   ON3      0.03    3     0.0 
    P4   ON2  P    ON3      0.10    2     0.0 
    P4   ON2  P    ON3      0.03    3     0.0 
    P    ON2  P4   ON4      0.10    2     0.0 
    P    ON2  P4   ON4      0.03    3     0.0 
    X    ON4  P4   X        0.30    3     0.0 
    NH2  CT1  C    O        0.0000  1     0.00
    NH2  CT2  C    O        0.0000  1     0.00  
    NH2  CT1  C    NH1      0.0000  1     0.00
    NH2  CT2  C    NH1      0.0000  1     0.00  
    NH2  CT2  C    N        0.0000  1     0.00  
    NH2  CT1  CT2A HA2      0.2000  3     0.00
    H    NH2  CT1  CT1      0.0000  1     0.00
    H    NH2  CT1  C        0.0000  1     0.00
    H    NH2  CT1  CD       0.0000  1     0.00
    H    NH2  CT2  C        0.0000  1     0.00  
    H    NH2  CT2  CD       0.0000  1     0.00  
    H    NH2  CT1  HB1      0.1100  3     0.00  
    H    NH2  CT2  HB2      0.1100  3     0.00  
    H    NH2  CT1  CT2      0.1100  3     0.00  
    H    NH2  CT1  CT2A     0.1100  3     0.00  
    H    NH2  CT1  CT3      0.1100  3     0.00  
    CC   CT2A CT1  NH2      0.6800  1   180.00  
    CC   CT2A CT1  NH2      0.1000  2   180.00
    CC   CT2A CT1  NH2      0.3800  3     0.00
    CD   CT1  CT2A CC       1.6100  1   180.00  
    CD   CT1  CT2A CC       1.2900  2   180.00
    CD   CT1  CT2A CC       0.5900  3   180.00
    CAI  CA   CA   CAI      3.1000  2   180.00
    CA   CPT  CPT  CA       3.0000  2   180.00
    CAI  CPT  CPT  CAI      3.0000  2   180.00
    CA   CY   CPT  CA       3.0000  2   180.00
    CA   CY   CPT  CAI      3.0000  2   180.00
    CA   NY   CPT  CA       3.0000  2   180.00
    CPT  CA   CA   CA       3.0000  2   180.00
    CPT  CPT  CA   CA       3.0000  2   180.00
    CA   NY   CPT  CAI      3.0000  2   180.00
    CPT  CAI  CA   CA       3.0000  2   180.00
    CPT  CPT  CAI  CA       3.0000  2   180.00
    CPT  CPT  CY   CA       5.0000  2   180.00
    CPT  CPT  NY   CA       6.5000  2   180.00
    CT3  CY   CPT  CA       2.5000  2   180.00
    CT3  CY   CPT  CAI      2.5000  2   180.00
    CT3  CY   CPT  CPT      3.0000  2   180.00
    CT2  CY   CPT  CA       2.5000  2   180.00
    CT2  CY   CPT  CAI      2.5000  2   180.00
    CT2  CY   CPT  CPT      3.0000  2   180.00
    CY   CA   NY   CPT      6.0000  2   180.00
    CY   CPT  CA   CA       4.0000  2   180.00
    CY   CPT  CPT  CA       4.0000  2   180.00
    CY   CPT  CAI  CA       4.0000  2   180.00
    CY   CPT  CPT  CAI      4.0000  2   180.00
    H    NY   CA   CY       0.0500  2   180.00
    H    NY   CPT  CA       0.2000  2   180.00
    H    NY   CPT  CAI      0.2000  2   180.00
    H    NY   CPT  CPT      0.8500  2   180.00
    HP   CAI  CA   CA       4.2000  2   180.00
    HP   CA   CA   CPT      3.0000  2   180.00
    HP   CA   CPT  CPT      3.0000  2   180.00
    HP   CA   CPT  CY       4.0000  2   180.00
    HP   CA   CA   CAI      4.2000  2   180.00
    HP   CA   CAI  CPT      3.0000  2   180.00
    HP   CAI  CA   HP       2.4000  2   180.00
    HP   CAI  CPT  CPT      3.0000  2   180.00
    HP   CAI  CPT  CY       4.0000  2   180.00
    HP   CA   CY   CPT      2.8000  2   180.00
    HP   CA   CY   CT3      1.2000  2   180.00
    HP   CA   CY   CT2      1.2000  2   180.00
    HP   CA   NY   CPT      2.6000  2   180.00
    HP   CA   NY   H        0.4000  2   180.00
    HP   CY   CA   HP       1.0000  2   180.00
    HP   CY   CPT  CA       2.8000  2   180.00
    HP   CY   CPT  CAI      2.8000  2   180.00
    HP   CY   CPT  CPT      2.6000  2   180.00
    NY   CA   CY   CPT      5.0000  2   180.00
    NY   CA   CY   CT3      2.5000  2   180.00
    NY   CA   CY   CT2      2.5000  2   180.00
    NY   CA   CY   HP       3.5000  2   180.00
    NY   CPT  CA   CA       3.0000  2   180.00
    NY   CPT  CA   HP       3.0000  2   180.00
    NY   CPT  CPT  CA       4.0000  2   180.00
    NY   CPT  CAI  CA       3.0000  2   180.00
    NY   CPT  CAI  HP       3.0000  2   180.00
    NY   CPT  CPT  CAI      4.0000  2   180.00
    NY   CPT  CPT  CY       6.5000  2   180.00
    CT3  CT2  CY   CA       0.3800  2     0.00
    CT3  CT2  CY   CPT      0.2500  2   180.00
    CT3  CT2  CY   CPT      0.3000  3     0.00
    HA3  CT3  CY   CA       0.0100  3     0.00
    HA3  CT3  CY   CPT      0.2000  3     0.00
    HA2  CT2  CY   CA       0.0100  3     0.00
    HA2  CT2  CY   CPT      0.2000  3     0.00
    X    CS   SS   X        0.0000  3     0.20
    C    CT1  NH1  C        0.2000  1   180.00
    C    CT2  NH1  C        0.2000  1   180.00
    C    N    CP1  C        0.8000  3     0.00
    CA   CA   CA   CA       3.1000  2   180.00
    CC   CP1  N    C        0.8000  3     0.00
    CC   CT1  CT2  CA       0.0400  3     0.00
    CC   CT1  NH1  C        0.2000  1   180.00
    CC   CT2  NH1  C        2.0000  1   180.00
    CD   CP1  N    C        0.0000  1   180.00
    CD   CT1  NH1  C        0.2000  1   180.00
    CD   CT2  NH1  C        2.0000  1   180.00
    CE1  CE1  CT3  HA3      0.0300  3     0.00
    CE2  CE1  CT2  CT3      0.5000  1   180.00
    CE2  CE1  CT2  CT3      1.3000  3   180.00	
    CE2  CE1  CT2  HA2      0.1200  3     0.00	
    CE2  CE1  CT3  HA3      0.0500  3   180.00	
    CP1  C    N    CP1      2.7500  2   180.00
    CP1  C    N    CP1      0.3000  4     0.00
    CP2  CP1  N    C        0.8000  3     0.00
    CP2  CP3  N    C        0.0000  3   180.00
    CP2  CP3  N    CP1      0.1000  3     0.00
    CP2  CP3  NP   CP1      0.0800  3     0.00
    CP3  N    C    CP1      2.7500  2   180.00
    CP3  N    C    CP1      0.3000  4     0.00
    CP3  N    CP1  C        0.1000  3     0.00
    CP3  N    CP1  CC       0.1000  3     0.00
    CP3  N    CP1  CP2      0.1000  3     0.00
    CP3  NP   CP1  C        0.0800  3     0.00
    CP3  NP   CP1  CC       0.0800  3     0.00
    CP3  NP   CP1  CD       0.0800  3     0.00
    CP3  NP   CP1  CP2      0.0800  3     0.00
    CPH2 NR1  CPH1 CPH1    14.0000  2   180.00
    CPH2 NR2  CPH1 CPH1    14.0000  2   180.00
    CPH2 NR3  CPH1 CPH1    12.0000  2   180.00
    CT1  C    N    CP1      2.7500  2   180.00
    CT1  C    N    CP1      0.3000  4     0.00
    CT1  C    N    CP3      2.7500  2   180.00
    CT1  C    N    CP3      0.3000  4     0.00
    CT1  C    NH1  CT1      1.6000  1     0.00
    CT1  C    NH1  CT1      2.5000  2   180.00
    CT1  CT1  NH1  C        1.8000  1     0.00
    CT1  NH1  C    CP1      1.6000  1     0.00
    CT1  NH1  C    CP1      2.5000  2   180.00
    CT2  C    N    CP1      2.7500  2   180.00
    CT2  C    N    CP1      0.3000  4     0.00
    CT2  C    N    CP3      2.7500  2   180.00
    CT2  C    N    CP3      0.3000  4     0.00
    CT2  C    NH1  CT1      1.6000  1     0.00
    CT2  C    NH1  CT1      2.5000  2   180.00
    CT2  C    NH1  CT2      1.6000  1     0.00
    CT2  C    NH1  CT2      2.5000  2   180.00
    CT2  C    NH1  CT3      1.6000  1     0.00
    CT2  C    NH1  CT3      2.5000  2   180.00
    CT2  CA   CA   CA       3.1000  2   180.00
    CT2  CPH1 NR1  CPH2     3.0000  2   180.00
    CT2  CPH1 NR2  CPH2     3.0000  2   180.00
    CT2  CPH1 NR3  CPH2     2.5000  2   180.00
    CT2  CT1  NH1  C        1.8000  1     0.00
    CT2  CT2  CPH1 CPH1     0.4000  1     0.00
    CT2  CT2  CT2  CT2      0.10    2   180.00
    CT2  CT2  CT2  CT2      0.15    4     0.00
    CT2  CT2  CT2  CT2      0.10    6   180.00
    CT2  CT2  CT2  CT3      0.10    2   180.00
    CT2  CT2  CT2  CT3      0.15    4     0.00
    CT2  CT2  CT2  CT3      0.10    6   180.00
    CT2  CT2  NH1  C        1.8000  1     0.00
    CT2  NH1  C    CP1      1.6000  1     0.00
    CT2  NH1  C    CP1      2.5000  2   180.00
    CT2  NH1  C    CT1      1.6000  1     0.00
    CT2  NH1  C    CT1      2.5000  2   180.00
    CT2  SM   SM   CT2      1.0000  1     0.00
    CT2  SM   SM   CT2      4.1000  2     0.00
    CT2  SM   SM   CT2      0.9000  3     0.00
    CT3  C    N    CP1      2.7500  2   180.00
    CT3  C    N    CP1      0.3000  4     0.00
    CT3  C    N    CP3      2.7500  2   180.00
    CT3  C    N    CP3      0.3000  4     0.00
    CT3  C    NH1  CT1      1.6000  1     0.00
    CT3  C    NH1  CT1      2.5000  2   180.00
    CT3  C    NH1  CT2      1.6000  1     0.00
    CT3  C    NH1  CT2      2.5000  2   180.00
    CT3  C    NH1  CT3      1.6000  1     0.00
    CT3  C    NH1  CT3      2.5000  2   180.00
    CT3  CA   CA   CA       3.1000  2   180.00
    CT3  CE1  CE2  HE2      5.2000  2   180.00	
    CT3  CPH1 NR1  CPH2     3.0000  2   180.00
    CT3  CT1  NH1  C        1.8000  1     0.00
    CT3  CT2  CA   CA       0.2300  2   180.00
    CT3  CT2  CPH1 CPH1     0.2000  1     0.00
    CT3  CT2  CPH1 CPH1     0.2700  2     0.00
    CT3  CT2  CPH1 CPH1     0.0000  3     0.00
    CT3  CT2  S    CT3      0.2400  1   180.00
    CT3  CT2  S    CT3      0.3700  3     0.00
    CT3  NH1  C    CP1      1.6000  1     0.00
    CT3  NH1  C    CP1      2.5000  2   180.00
    CT3  NH1  C    CT1      1.6000  1     0.00
    CT3  NH1  C    CT1      2.5000  2   180.00
    CT3  S    CT2  CT2      0.2400  1   180.00
    CT3  S    CT2  CT2      0.3700  3     0.00
    CT3  SM   SM   CT3      1.0000  1     0.00
    CT3  SM   SM   CT3      4.1000  2     0.00
    CT3  SM   SM   CT3      0.9000  3     0.00
    H    NH1  C    CP1      2.5000  2   180.00
    H    NH1  C    CT1      2.5000  2   180.00
    H    NH1  C    CT2      2.5000  2   180.00
    H    NH1  C    CT3      2.5000  2   180.00
    H    NH1  CT1  C        0.0000  1     0.00
    H    NH1  CT1  CC       0.0000  1     0.00
    H    NH1  CT1  CD       0.0000  1     0.00
    H    NH1  CT1  CT1      0.0000  1     0.00
    H    NH1  CT1  CT2      0.0000  1     0.00
    H    NH1  CT1  CT3      0.0000  1     0.00
    H    NH1  CT2  C        0.0000  1     0.00
    H    NH1  CT2  CC       0.0000  1     0.00
    H    NH1  CT2  CD       0.0000  1     0.00
    H    NH1  CT2  CT2      0.0000  1     0.00
    H    NH1  CT2  CT3      0.0000  1     0.00
    H    NH2  CC   CT1      1.4000  2   180.00
    H    NH2  CC   CT2      1.4000  2   180.00
    H    NH2  CC   CT3      1.4000  2   180.00
    H    NH2  CC   CP1      2.5000  2   180.00
    H    NR1  CPH1 CPH1     1.0000  2   180.00
    H    NR1  CPH1 CT2      1.0000  2   180.00
    H    NR1  CPH1 CT3      1.0000  2   180.00
    H    NR3  CPH1 CPH1     1.4000  2   180.00
    H    NR3  CPH1 CT2      3.0000  2   180.00
    H    NR3  CPH1 CT3      3.0000  2   180.00
    H    OH1  CA   CA       0.9900  2   180.00
    H    OH1  CT1  CT3      1.3300  1     0.00
    H    OH1  CT1  CT3      0.1800  2     0.00
    H    OH1  CT1  CT3      0.3200  3     0.00
    H    OH1  CT2  CT2      1.3000  1     0.00
    H    OH1  CT2  CT2      0.3000  2     0.00
    H    OH1  CT2  CT2      0.4200  3     0.00
    H    OH1  CT2  CT3      1.3000  1     0.00
    H    OH1  CT2  CT3      0.3000  2     0.00
    H    OH1  CT2  CT3      0.4200  3     0.00
    HA1  CC   NH2  H        1.4000  2   180.00
    HA2  CP3  N    C        0.0000  3   180.00
    HA2  CP3  N    CP1      0.1000  3     0.00
    HA2  CP3  NP   CP1      0.0800  3     0.00
    HA1  CT1  CT2  CA       0.0400  3     0.00
    HA2  CT2  CPH1 CPH1     0.0000  3     0.00
    HA2  CT2  NH1  C        0.0000  3     0.00
    HA2  CT2  NH1  H        0.0000  3     0.00
    HA2  CT2  S    CT3      0.2800  3     0.00
    HA3  CT3  CPH1 CPH1     0.0000  3     0.00
    HA3  CT3  CS   HA2      0.1600  3     0.00
    HA3  CT3  CS   HA3      0.1600  3     0.00
    HA3  CT3  CT2  CA       0.0400  3     0.00
    HA3  CT3  NH1  C        0.0000  3     0.00
    HA3  CT3  NH1  H        0.0000  3     0.00
    HA3  CT3  S    CT2      0.2800  3     0.00
    HE1  CE1  CE1  HE1      1.0000  2   180.00
    CT3  CE1  CE1  HE1      1.0000  2   180.00
    HE1  CE1  CE2  HE2      5.2000  2   180.00
    HE1  CE1  CT2  HA2      0.0000  3     0.00
    HE1  CE1  CT2  CT3      0.1200  3     0.00	
    HE1  CE1  CT3  HA3      0.0000  3     0.00	
    HE2  CE2  CE1  CT2      5.2000  2   180.00	
    HB1  CP1  N    C        0.8000  3     0.00
    HB1  CP1  N    CP3      0.1000  3     0.00
    HB1  CP1  NP   CP3      0.0800  3     0.00
    HB1  CT1  NH1  C        0.0000  1     0.00
    HB1  CT1  NH1  H        0.0000  1     0.00
    HB2  CT2  NH1  C        0.0000  1     0.00
    HB2  CT2  NH1  H        0.0000  1     0.00
    HC   NH2  CT2  HB2      0.1100  3     0.00
    HC   NH2  CT2  CD       0.1100  3     0.00
    HC   NH2  CT2  CT2       0.1100  3     0.00
    HC   NH2  CT2  HA2      0.1100  3     0.00
    HC   NP   CP1  C        0.0800  3     0.00
    HC   NP   CP1  CC       0.0800  3     0.00
    HC   NP   CP1  CD       0.0800  3     0.00
    HC   NP   CP1  CP2      0.0800  3     0.00
    HC   NP   CP1  HB1      0.0800  3     0.00
    HC   NP   CP3  CP2      0.0800  3     0.00
    HC   NP   CP3  HA2      0.0800  3     0.00
    HP   CA   CA   CA       4.2000  2   180.00
    HP   CA   CA   CT2      4.2000  2   180.00
    HP   CA   CA   CT3      4.2000  2   180.00
    HP   CA   CA   HP       2.4000  2   180.00
    HR1  CPH1 CPH1 CT2      1.0000  2   180.00
    HR1  CPH1 CPH1 CT3      1.0000  2   180.00
    HR1  CPH1 CPH1 HR1      1.0000  2   180.00
    HR1  CPH1 NR3  CPH2     2.5000  2   180.00
    HR1  CPH1 NR3  H        3.0000  2   180.00
    HR1  CPH2 NR1  CPH1     3.0000  2   180.00
    HR1  CPH2 NR1  H        1.0000  2   180.00
    HR1  CPH2 NR2  CPH1     3.0000  2   180.00
    HR2  CPH2 NR3  CPH1     3.0000  2   180.00
    HR2  CPH2 NR3  H        0.0000  2   180.00
    HR3  CPH1 CPH1 CT2      2.0000  2   180.00
    HR3  CPH1 CPH1 CT3      2.0000  2   180.00
    HR3  CPH1 CPH1 HR3      2.0000  2   180.00
    HR3  CPH1 NR1  CPH2     3.0000  2   180.00
    HR3  CPH1 NR1  H        1.0000  2   180.00
    HR3  CPH1 NR2  CPH2     3.0000  2   180.00
    HS   S    CT2  CT3      0.2400  1     0.00
    HS   S    CT2  CT3      0.1500  2     0.00
    HS   S    CT2  CT3      0.2700  3     0.00
    HS   S    CT2  HA2      0.2000  3     0.00
    HS   S    CT3  HA3      0.2000  3     0.00
    N    C    CP1  CP2      0.4000  1     0.00
    N    C    CP1  CP2      0.6000  2     0.00
    N    C    CP1  HB1      0.4000  1   180.00
    N    C    CP1  HB1      0.6000  2     0.00
    N    C    CP1  N        0.3000  1     0.00
    N    C    CP1  N       -0.3000  4     0.00
    N    C    CT1  CT1      0.0000  1     0.00
    N    C    CT1  CT2      0.0000  1     0.00
    N    C    CT1  CT3      0.0000  1     0.00
    N    C    CT1  HB1      0.0000  1     0.00
    N    C    CT2  HB2      0.0000  1     0.00
    N    C    CT3  HA3      0.0000  1     0.00
    N    CT1  CT2  CA       0.0400  3     0.00
    NH1  C    CP1  CP2      0.4000  1     0.00
    NH1  C    CP1  CP2      0.6000  2     0.00
    NH1  C    CP1  HB1      0.4000  1   180.00
    NH1  C    CP1  HB1      0.6000  2     0.00
    NH1  C    CP1  N        0.3000  1     0.00
    NH1  C    CP1  N       -0.3000  4     0.00
    NH1  C    CT1  CT1      0.0000  1     0.00
    NH1  C    CT1  CT2      0.0000  1     0.00
    NH1  C    CT1  CT3      0.0000  1     0.00
    NH1  C    CT1  HB1      0.0000  1     0.00
    NH1  C    CT1  NH1      0.6000  1     0.00
    NH1  C    CT2  CT2      0.0000  1     0.00
    NH1  C    CT2  HA2      0.0000  3     0.00
    NH1  C    CT2  HB2      0.0000  1     0.00
    NH1  C    CT2  NH1      0.6000  1     0.00
    NH1  C    CT3  HA3      0.0000  3     0.00
    NH1  CT1  C    N        0.4000  1     0.00
    NH1  CT2  C    N        0.4000  1     0.00
    NH2  CC   CP1  CP2      0.4000  1     0.00
    NH2  CC   CP1  CP2      0.6000  2     0.00
    NH2  CC   CP1  HB1      0.4000  1   180.00
    NH2  CC   CP1  HB1      0.6000  2     0.00
    NH2  CC   CP1  N        0.3000  1     0.00
    NH2  CC   CP1  N       -0.3000  4     0.00
    NH2  CC   CT2  HA2      0.0000  3   180.00
    NH3  CT1  C    N        0.4000  1     0.00
    NH3  CT1  C    NH1      0.6000  1     0.00
    NH3  CT1  CC   NH2      0.4000  1     0.00
    NH3  CT2  C    N        0.4000  1     0.00
    NH3  CT2  C    NH1      1.0000  1     0.00
    NH3  CT2  CC   NH2      0.4000  1     0.00
    NP   CP1  C    N        0.3000  1     0.00
    NP   CP1  C    NH1      0.3000  1     0.00
    NP   CP1  CC   NH2      0.3000  1     0.00
    NR1  CPH1 CPH1 CT2      3.0000  2   180.00
    NR1  CPH1 CPH1 CT3      3.0000  2   180.00
    NR1  CPH1 CPH1 HR3      3.0000  2   180.00
    NR1  CPH1 CT2  CT2      0.1900  3     0.00
    NR1  CPH1 CT2  CT3      0.1900  3     0.00
    NR1  CPH1 CT2  HA2      0.1900  3     0.00
    NR1  CPH1 CT3  HA3      0.1900  3     0.00
    NR1  CPH2 NR2  CPH1    14.0000  2   180.00
    NR2  CPH1 CPH1 CT2      3.0000  2   180.00
    NR2  CPH1 CPH1 CT3      3.0000  2   180.00
    NR2  CPH1 CPH1 HR3      3.0000  2   180.00
    NR2  CPH1 CPH1 NR1     14.0000  2   180.00
    NR2  CPH1 CT2  CT2      0.1900  3     0.00
    NR2  CPH1 CT2  CT3      0.1900  3     0.00
    NR2  CPH1 CT2  HA2      0.1900  3     0.00
    NR2  CPH1 CT3  HA3      0.1900  3     0.00
    NR2  CPH2 NR1  CPH1    14.0000  2   180.00
    NR2  CPH2 NR1  H        1.0000  2   180.00
    NR3  CPH1 CPH1 CT2      2.5000  2   180.00
    NR3  CPH1 CPH1 CT3      2.5000  2   180.00
    NR3  CPH1 CPH1 HR1      2.5000  2   180.00
    NR3  CPH1 CPH1 NR3     12.0000  2   180.00
    NR3  CPH1 CT2  CT2      0.1900  3     0.00
    NR3  CPH1 CT2  CT3      0.1900  3     0.00
    NR3  CPH1 CT2  HA2      0.1900  3     0.00
    NR3  CPH1 CT3  HA3      0.1900  3     0.00
    NR3  CPH2 NR3  CPH1    12.0000  2   180.00
    NR3  CPH2 NR3  H        1.4000  2   180.00
    O    C    CP1  CP2      0.4000  1   180.00
    O    C    CP1  CP2      0.6000  2     0.00
    O    C    CP1  HB1      0.4000  1     0.00
    O    C    CP1  HB1      0.6000  2     0.00
    O    C    CP1  N       -0.3000  4     0.00
    O    C    CT1  CT1      1.4000  1     0.00
    O    C    CT1  CT2      1.4000  1     0.00
    O    C    CT1  CT3      1.4000  1     0.00
    O    C    CT1  HB1      0.0000  1     0.00
    O    C    CT1  NH1      0.0000  1     0.00
    O    C    CT1  NH3      0.0000  1     0.00
    O    C    CT2  CT2      1.4000  1     0.00
    O    C    CT2  HA2      0.0000  3   180.00
    O    C    CT2  HB2      0.0000  1     0.00
    O    C    CT2  NH1      0.0000  1     0.00
    O    C    CT2  NH3      0.0000  1     0.00
    O    C    CT3  HA3      0.0000  3   180.00
    O    C    N    CP1      2.7500  2   180.00
    O    C    N    CP1      0.3000  4     0.00
    O    C    N    CP3      2.7500  2   180.00
    O    C    N    CP3      0.3000  4     0.00
    O    C    NH1  CT1      2.5000  2   180.00
    O    C    NH1  CT2      2.5000  2   180.00
    O    C    NH1  CT3      2.5000  2   180.00
    O    C    NH1  H        2.5000  2   180.00
    O    CC   CP1  CP2      0.4000  1   180.00
    O    CC   CP1  CP2      0.6000  2     0.00
    O    CC   CP1  HB1      0.4000  1     0.00
    O    CC   CP1  HB1      0.6000  2     0.00
    O    CC   CP1  N       -0.3000  4     0.00
    O    CC   CT2  HA2      0.0000  3   180.00
    O    CC   NH2  H        1.4000  2   180.00
    OB   CD   OS   CT2      0.9650  1   180.00
    OB   CD   OS   CT2      3.8500  2   180.00
    OB   CD   OS   CT3      0.9650  1   180.00
    OB   CD   OS   CT3      3.8500  2   180.00
    OC   CA   CA   CA       3.1000  2   180.00
    OC   CA   CA   HP       4.2000  2   180.00
    OC   CC   CP1  CP2      0.1600  3     0.00
    OC   CC   CP1  HB1      0.1600  3     0.00
    OC   CC   CP1  N        0.1600  3     0.00
    OC   CC   CP1  NP       0.1600  3     0.00
    OC   CC   CT1  NH3      3.2000  2   180.00
    OC   CC   CT2  NH3      3.2000  2   180.00
    OH1  CA   CA   CA       3.1000  2   180.00
    OH1  CA   CA   HP       4.2000  2   180.00
    S    CT2  CT2  HA2      0.0100  3     0.00
    SM   CT2  CT2  HA2      0.0100  3     0.00
    SM   SM   CT2  CT1      0.3100  3     0.00
    SM   SM   CT2  CT2      0.3100  3     0.00
    SM   SM   CT2  CT3      0.3100  3     0.00
    SM   SM   CT2  HA2      0.1580  3     0.00
    SM   SM   CT3  HA3      0.1580  3     0.00
    SS   CS   CT3  HA3      0.1500  3     0.00
    X    C    NC2  X        2.2500  2   180.00
    X    CD   OH1  X        2.0500  2   180.00
    X    CD   OS   X        2.0500  2   180.00
    X    CE1  CE1  X        0.1500  1     0.00
    X    CE1  CE1  X        8.5000  2   180.00
    X    CE2  CE2  X        4.9000  2   180.00	
    X    CP1  C    X        0.0000  6   180.00
    X    CP1  CC   X        0.0000  6   180.00
    X    CP1  CD   X        0.0000  6   180.00
    X    CP1  CP2  X        0.1400  3     0.00
    X    CP2  CP2  X        0.1600  3     0.00
    X    CP3  CP2  X        0.1400  3     0.00
    X    CT1  CC   X        0.0500  6   180.00
    X    CT1  CD   X        0.0000  6   180.00
    X    CT1  CT1  X        0.2000  3     0.00
    X    CT1  CT2  X        0.2000  3     0.00
    X    CT1  CT3  X        0.2000  3     0.00
    X    CT1  NH3  X        0.1000  3     0.00
    X    CT1  OH1  X        0.1400  3     0.00
    X    CT1  OS   X       -0.1000  3     0.00
    X    CT2  CA   X        0.0000  6     0.00
    X    CT2  CC   X        0.0500  6   180.00
    X    CT2  CD   X        0.0000  6   180.00
    X    CT2  CT2  X        0.1900  3     0.00
    X    CT2  CT3  X        0.1600  3     0.00
    X    CT2  NC2  X        0.0000  6   180.00
    X    CT2  NH3  X        0.1000  3     0.00
    X    CT2  OH1  X        0.1400  3     0.00
    X    CT2  OS   X       -0.1000  3     0.00
    X    CT3  CA   X        0.0000  6     0.00
    X    CT3  CC   X        0.0500  6   180.00
    X    CT3  CD   X        0.0000  6   180.00
    X    CT3  CT3  X        0.1525  3     0.00
    X    CT3  NC2  X        0.0000  6   180.00
    X    CT3  NH2  X        0.1100  3     0.00
    X    CT3  NH3  X        0.0900  3     0.00
    X    CT3  OH1  X        0.1400  3     0.00
    X    CT3  OS   X       -0.1000  3     0.00
    NH1  CT1  CT1  HA1      0.2000  3     0.00
    HB1  CT1  CT1  HA1      0.2000  3     0.00
    HB1  CT1  CT1  CT3      0.2000  3     0.00
    HA1  CT1  CT1  C        0.2000  3     0.00
    NH1  CT1  CT2  HA2      0.2000  3     0.00
    HB1  CT1  CT2  HA2      0.2000  3     0.00
    HB1  CT1  CT2  OH1      0.2000  3     0.00
    HB1  CT1  CT2  CT2      0.2000  3     0.00
    HA2  CT2  CT1  C        0.2000  3     0.00
    HA2  CT2  OH1  H        0.1400  3     0.00
    CT1  CT2  CT2  HA2      0.1900  3     0.00
    HA2  CT2  CT2  HA2      0.1900  3     0.00
    HA2  CT2  CT2  CC       0.1900  3     0.00
    HB1  CT1  CT2  S        0.2000  3     0.00
    CT2  CT2  CT2  HA2      0.1900  3     0.00
    CT2  CT2  CT2  NC2      0.1900  3     0.00
    CT2  CT2  NC2  HC       0.0000  6   180.00
    CT2  CT2  NC2  C        0.0000  6   180.00
    HA2  CT2  CT2  NC2      0.1900  3     0.00
    CT2  NC2  C    NC2      2.2500  2   180.00
    HA2  CT2  NC2  HC       0.0000  6   180.00
    HA2  CT2  NC2  C        0.0000  6   180.00
    NC2  C    NC2  HC       2.2500  2   180.00
    HB1  CT1  CT2  CC       0.2000  3     0.00
    HB1  CT1  CT2  CY       0.2000  3     0.00
    HA2  CT2  CC   OC       0.0500  6   180.00
    HB1  CT1  CT2  CPH1     0.2000  3     0.00
    CT1  CT1  CT3  HA3      0.2000  3     0.00
    CT1  CT1  CT2  HA2      0.2000  3     0.00
    HB1  CT1  CT1  CT2      0.2000  3     0.00
    CT1  CT2  CT3  HA3      0.1600  3     0.00
    HA1  CT1  CT3  HA3      0.2000  3     0.00
    HA1  CT1  CT2  HA2      0.2000  3     0.00
    HA1  CT1  CT2  CT3      0.2000  3     0.00
    CT3  CT1  CT2  HA2      0.2000  3     0.00
    CT3  CT1  CT2  CT3      0.2000  3     0.00
    HA3  CT3  CT1  CT2      0.2000  3     0.00
    HA2  CT2  CT3  HA3      0.1600  3     0.00
    CT1  CT2  CT1  HA1      0.2000  3     0.00
    HB1  CT1  CT2  CT1      0.2000  3     0.00
    CT3  CT1  CT3  HA3      0.2000  3     0.00
    CT2  CT2  CT2  NH3      0.1900  3     0.00
    CT2  CT2  NH3  HC       0.1000  3     0.00
    HA2  CT2  CT2  NH3      0.1900  3     0.00
    HA2  CT2  NH3  HC       0.1000  3     0.00
    HB1  CT1  CT2  CA       0.2000  3     0.00
    HA2  CT2  CA   CA       0.0000  6     0.00
    HB1  CT1  CT1  OH1      0.2000  3     0.00
    HA1  CT1  OH1  H        0.1400  3     0.00
    OH1  CT1  CT3  HA3      0.2000  3     0.00
    CT2  CT2  CC   O        0.0500  6   180.00
    CT2  CT2  CC   NH2      0.0500  6   180.00
    CT2  CT2  CC   OC       0.0500  6   180.00
    NH1  CT1  CT2A HA2      0.2000  3     0.00
    NH3  CT1  CT2A CT2      0.2000  3     0.00
    CT1  CT2A CT2  HA2      0.1900  3     0.00
    HB1  CT1  CT2A HA2      0.2000  3     0.00
    HB1  CT1  CT2A CT2      0.2000  3     0.00
    HA2  CT2A CT1  C        0.2000  3     0.00
    HA2  CT2A CT1  CC       0.2000  3     0.00
    HA2  CT2A CT2  HA2      0.1900  3     0.00
    HA2  CT2A CT2  CC       0.1900  3     0.00
    HB1  CT1  CT2A CPH1     0.2000  3     0.00
    C    NH1  CT1  CT2A     1.8000  1     0.00
    H    NH1  CT1  CT2A     0.0000  1     0.00
    CT2A CT1  C    O        1.4000  1     0.00
    CT2A CT1  C    NH1      0.0000  1     0.00
    CT2A CT1  C    N        0.0000  1     0.00
    CT1  CT2A CT2  CD       0.1900  3     0.00
    HA2  CT2A CT2  CD       0.1900  3     0.00
    CT2A CPH1 CPH1 HR1      1.0000  2   180.00
    CT2A CPH1 CPH1 NR3      2.5000  2   180.00
    CT2A CPH1 NR3  H        3.0000  2   180.00
    CT2A CPH1 NR3  CPH2     2.5000  2   180.00
    HA2  CT2A CPH1 CPH1     0.0000  3     0.00
    HA2  CT2A CPH1 NR3      0.1900  3     0.00
    C    CT1  CT2  CT2      0.3500  1   180.00 
    C    CT1  CT2  CT2      0.4200  2   180.00 
    C    CT1  CT2  CT2      1.9100  3   180.00 
    CT2  CT2  CT1  NH1      0.8800  1   180.00 
    CT2  CT2  CT1  NH1      0.0000  2   180.00 
    CT2  CT2  CT1  NH1      1.9000  3     0.00 
    CC   CT2  CT2  CT1      1.8400  1   180.00 
    CC   CT2  CT2  CT1      0.8400  2   180.00 
    CC   CT2  CT2  CT1      0.3900  3   180.00 
    CT1  CT2  CT2  CT2      0.6300  1   180.00 
    CT1  CT2  CT2  CT2      0.0100  2     0.00 
    CT1  CT2  CT2  CT2      0.1500  3     0.00 
    CT1  CT2  CT2  S        0.1400  1   180.00 
    CT1  CT2  CT2  S        0.5400  2     0.00 
    CT1  CT2  CT2  S        0.6900  3     0.00 
    C    CT1  CT2  CC       1.4100  1   180.00 
    C    CT1  CT2  CC       1.2900  2   180.00 
    C    CT1  CT2  CC       0.5900  3   180.00 
    CC   CT2  CT1  NH1      0.2800  1   180.00 
    CC   CT2  CT1  NH1      0.5000  2   180.00 
    CC   CT2  CT1  NH1      0.3800  3     0.00 
    CT1  CT2  CC   NH2      0.6200  1   180.00 
    CT1  CT2  CC   NH2      0.6600  2   180.00 
    CT1  CT2  CC   NH2      0.7200  3   180.00 
    CT1  CT2  CC   O        0.4200  1   180.00 
    CT1  CT2  CC   O        0.1500  2   180.00 
    CT1  CT2  CC   O        0.9500  3   180.00 
    C    CT1  CT2A CC       1.6100  1   180.00 
    C    CT1  CT2A CC       1.2900  2   180.00 
    C    CT1  CT2A CC       0.5900  3   180.00 
    CC   CT2A CT1  NH1      0.6800  1   180.00 
    CC   CT2A CT1  NH1      0.1000  2   180.00 
    CC   CT2A CT1  NH1      0.3800  3     0.00 
    CT1  CT2A CC   OC       0.8400  1     0.00 
    CT1  CT2A CC   OC       0.9800  2   180.00 
    CT1  CT2A CC   OC       1.4600  3     0.00 
    CT1  CT2  S    HS       0.2000  1     0.00
    CT1  CT2  S    HS       0.6500  2     0.00
    CT1  CT2  S    HS       0.2200  3     0.00
    C    CT1  CT2  S        0.2400  1   180.00
    C    CT1  CT2  S        0.7500  2   180.00
    C    CT1  CT2  S        1.3500  3   180.00
    NH1  CT1  CT2  S        0.3400  1     0.00
    NH1  CT1  CT2  S        0.5000  2   180.00
    NH1  CT1  CT2  S        1.4300  3     0.00
    CC   CT2  CT2A CT1      0.0000  1   180.00
    CC   CT2  CT2A CT1      0.3800  2   180.00
    CC   CT2  CT2A CT1      0.5900  3   180.00
    C    CT1  CT2A CT2      0.1100  1     0.00
    C    CT1  CT2A CT2      0.9800  2   180.00
    C    CT1  CT2A CT2      1.6000  3   180.00
    CC   CT1  CT2A CT2      1.6000  3   180.00
    CT2  CT2A CT1  NH1      0.3000  1     0.00 
    CT2  CT2A CT1  NH1      0.3500  2     0.00
    CT2  CT2A CT1  NH1      1.7600  3     0.00
    CPH1 CPH1 CT2  CT1      1.7400  1     0.00
    CPH1 CPH1 CT2  CT1      0.1500  2     0.00
    CPH1 CPH1 CT2  CT1      0.7700  3   180.00
    CT1  CT2  CPH1 NR1      1.4900  1     0.00
    CT1  CT2  CPH1 NR1      0.0900  2   180.00
    CT1  CT2  CPH1 NR1      0.7900  3   180.00
    CT1  CT2  CPH1 NR2      1.0900  1     0.00
    CT1  CT2  CPH1 NR2      0.0900  2     0.00
    CT1  CT2  CPH1 NR2      0.6700  3   180.00
    C    CT1  CT2  CPH1     0.1800  1   180.00
    C    CT1  CT2  CPH1     0.6400  2   180.00
    C    CT1  CT2  CPH1     0.8700  3   180.00
    CPH1 CT2  CT1  NH1      0.0000  1     0.00
    CPH1 CT2  CT1  NH1      0.0000  2   180.00
    CPH1 CT2  CT1  NH1      0.9000  3     0.00
    CPH1 CPH1 CT2A CT1      2.0400  1     0.00
    CPH1 CPH1 CT2A CT1      0.4400  2     0.00
    CPH1 CPH1 CT2A CT1      0.1300  3   180.00
    CT1  CT2A CPH1 NR3      0.5300  1   180.00
    CT1  CT2A CPH1 NR3      0.4200  2   180.00
    CT1  CT2A CPH1 NR3      0.3000  3   180.00
    C    CT1  CT2A CPH1     1.7500  1   180.00
    C    CT1  CT2A CPH1     0.1300  2     0.00
    C    CT1  CT2A CPH1     1.8600  3   180.00
    CPH1 CT2A CT1  NH1      1.0900  1   180.00
    CPH1 CT2A CT1  NH1      0.2200  2   180.00
    CPH1 CT2A CT1  NH1      2.3200  3     0.00
    CT1  CT1  CT2  CT3      0.3800  1   180.00
    CT1  CT1  CT2  CT3      0.1300  2   180.00
    CT1  CT1  CT2  CT3      0.2900  3   180.00
    C    CT1  CT1  CT2      0.1000  1   180.00
    C    CT1  CT1  CT2      0.5200  2   180.00
    C    CT1  CT1  CT2      0.2900  3   180.00
    CT2  CT1  CT1  NH1      0.1200  1   180.00
    CT2  CT1  CT1  NH1      0.3600  2   180.00
    CT2  CT1  CT1  NH1      0.4100  3     0.00
    CT1  CT2  CT1  CT3      0.0500  1     0.00
    CT1  CT2  CT1  CT3      0.1000  2   180.00
    CT1  CT2  CT1  CT3      0.0100  3   180.00
    C    CT1  CT2  CT1      0.3200  1   180.00
    C    CT1  CT2  CT1      0.6100  2   180.00
    C    CT1  CT2  CT1      0.7200  3   180.00
    CT1  CT2  CT1  NH1      0.4800  1   180.00
    CT1  CT2  CT1  NH1      0.4200  2   180.00
    CT1  CT2  CT1  NH1      0.6500  3     0.00
    CA   CA   CT2  CT1      1.0700  1     0.00
    CA   CA   CT2  CT1      0.2400  2   180.00
    CA   CA   CT2  CT1      0.1700  3   180.00
    C    CT1  CT2  CA       1.2800  1   180.00
    C    CT1  CT2  CA       0.9400  2   180.00
    C    CT1  CT2  CA       1.5700  3   180.00
    CA   CT2  CT1  NH1      0.5200  1   180.00
    CA   CT2  CT1  NH1      0.6200  2   180.00
    CA   CT2  CT1  NH1      1.5800  3     0.00
    CT1  CT2  OH1  H        0.0200  1     0.00
    CT1  CT2  OH1  H        0.5600  2     0.00
    CT1  CT2  OH1  H        0.4900  3     0.00
    C    CT1  CT2  OH1      0.6500  1   180.00
    C    CT1  CT2  OH1      0.2500  2   180.00
    C    CT1  CT2  OH1      1.1700  3   180.00
    NH1  CT1  CT2  OH1      0.1800  1   180.00
    NH1  CT1  CT2  OH1      0.1900  2   180.00
    NH1  CT1  CT2  OH1      1.4600  3     0.00
    CT1  CT1  OH1  H        0.1800  1     0.00
    CT1  CT1  OH1  H        0.0600  2     0.00
    CT1  CT1  OH1  H        0.2500  3     0.00
    C    CT1  CT1  OH1      0.7900  1   180.00
    C    CT1  CT1  OH1      0.3900  2   180.00
    C    CT1  CT1  OH1      0.9900  3   180.00
    NH1  CT1  CT1  OH1      0.0900  1     0.00
    NH1  CT1  CT1  OH1      0.1900  2   180.00
    NH1  CT1  CT1  OH1      0.1700  3     0.00
    CA   CY   CT2  CT1      0.0300  1     0.00
    CA   CY   CT2  CT1      0.5500  2     0.00
    CA   CY   CT2  CT1      0.3900  3   180.00
    CPT  CY   CT2  CT1      0.3600  1   180.00
    CPT  CY   CT2  CT1      0.0500  2     0.00
    CPT  CY   CT2  CT1      0.1900  3   180.00
    C    CT1  CT2  CY       1.0900  1   180.00
    C    CT1  CT2  CY       0.5000  2   180.00
    C    CT1  CT2  CY       1.1700  3   180.00
    CY   CT2  CT1  NH1      0.2900  1   180.00
    CY   CT2  CT1  NH1      0.6600  2   180.00
    CY   CT2  CT1  NH1      1.1700  3     0.00
    C    CT1  CT1  CT3      0.1400  1   180.00
    C    CT1  CT1  CT3      0.2600  2   180.00
    C    CT1  CT1  CT3      0.3300  3   180.00
    CT3  CT1  CT1  NH1      0.1800  1     0.00
    CT3  CT1  CT1  NH1      0.0600  2     0.00
    CT3  CT1  CT1  NH1      0.5900  3     0.00
    H    NH1  CT2A CC       0.0000  1     0.00
    X    CT2A CC   X        0.0500  6   180.00
    HB1  CT1  CT2A CC       0.2000  3     0.00
    HA2  CT2A CC   OC       0.0500  6   180.00
    NH3  CT1  CT2A HA2      0.2000  3     0.00
    NH3  CT1  CT2A CC       0.2000  3     0.00
    CC   CT2A CT1  CC       0.2000  3     0.00
    CPH1 CT2A CT1  CC       0.2000  3     0.00
    CPH1 CT2A CT1  NH3      0.2000  3     0.00
    CPH1 CT2A CT1  CD       0.2000  3     0.00 
    HA2  CT2A CT1  CD       0.2000  3     0.00     
    CT2  CT2A CT1  CD       0.2000  3     0.00
    H    NH2  CT1  CS       0.1100  3     0.00
    CS   CT1  NH1  C        1.8000  1     0.00
    H    NH1  CT1  CS       0.0000  1     0.00
    N    C    CT1  CS       0.0000  1     0.00
    NH1  C    CT1  CS       0.0000  1     0.00
    O    C    CT1  CS       1.4000  1     0.00
    HA2  CS   CT1  C        0.2000  3     0.00
    NH1  CT1  CS   HA2      0.2000  3     0.00
    HB1  CT1  CS   HA2      0.2000  3     0.00
    HB1  CT1  CS   SS       0.2000  3     0.00
    C    CT1  CS   SS       0.2000  3     0.00
    NH1  CT1  CS   SS       0.2000  3     0.00
    NH3  CT1  CS   HA2      0.2000  3     0.00
    NH3  CT1  CS   SS       0.2000  3     0.00
    NH2  CT1  CS   HA2      0.2000  3     0.00
    NH2  CT1  CS   SS       0.2000  3     0.00
    CC   CT1  CS   HA2      0.2000  3     0.00
    CC   CT1  CS   SS       0.2000  3     0.00
    CD   CT1  CS   HA2      0.2000  3     0.00
    CD   CT1  CS   SS       0.2000  3     0.00
    NH1  CT1  CT2  OC       0.2000  3     0.00
    NH2  CT1  CT2  OC       0.2000  3     0.00
    NH3  CT1  CT2  OC       0.2000  3     0.00
    C    CT1  CT2  OC       0.2000  3     0.00
    CC   CT1  CT2  OC       0.2000  3     0.00
    CD   CT1  CT2  OC       0.2000  3     0.00
    HB1  CT1  CT2  OC       0.2000  3     0.00
    HT  OTMG MG    OTMG     0.000   6    180.0000
    HT  OTCA CAL   OTCA     0.000   6    180.0000", "\n")

improper_para_info = split("OBL  X    X    CL         100.00       0.00 
    HEL2 HEL2 CEL2 CEL2         3.00         0.00 
    OCL  X    X    CL          96.00         0.00 
    OCL  X    X    CCL         96.00         0.00 
    HN2  X    X    NN2      1.0        0.0     
    NN2B CN4  CN5  HN2      7.0        0.0     
    HN1  X    X    NN1      4.0        0.0     
    NN1  CN2  HN1  HN1      6.0        0.0     
    CN1  X    X    ON1     90.0        0.0     
    CN1T X    X    ON1     90.0        0.0     
    CN1  NN2G CN5G ON1     90.0        0.0     
    CN1T NN2B NN2U ON1    110.0        0.0     
    CN1  NN2U CN3T ON1     90.0        0.0     
    CN1  X    X    ON1C    80.0        0.0     
    CN2  X    X    NN1     90.0        0.0     
    CN2  NN3G NN2G NN1     40.0        0.0     
    CN2  NN3A CN5  NN1     40.0        0.0     
    CN2  NN3  CN3  NN1     60.0        0.0     
    CN9  X    X    CN3T    14.0        0.0     
    HN2  X    X    NN2      1.0        0.0     
    NN2B CN4  CN5  HN2      7.0        0.0     
    HN1  X    X    NN1      4.0        0.0     
    NN1  CN2  HN1  HN1      6.0        0.0     
    CN1  X    X    ON1     90.0        0.0     
    CN1T X    X    ON1     90.0        0.0     
    CN1  NN2G CN5G ON1     90.0        0.0     
    CN1T NN2B NN2U ON1    110.0        0.0     
    CN1  NN2U CN3T ON1     90.0        0.0     
    CN1  X    X    ON1C    80.0        0.0     
    CN2  X    X    NN1     90.0        0.0     
    CN2  NN3G NN2G NN1     40.0        0.0     
    CN2  NN3A CN5  NN1     40.0        0.0     
    CN2  NN3  CN3  NN1     60.0        0.0     
    CN9  X    X    CN3T    14.0        0.0     
    HE2  HE2  CE2  CE2     3.0           0.00   		
    HR1  NR1  NR2  CPH2    0.5000        0.0000                
    HR1  NR2  NR1  CPH2    0.5000        0.0000                
    HR3  CPH1 NR1  CPH1    0.5000        0.0000                
    HR3  CPH1 NR2  CPH1    0.5000        0.0000                
    HR3  CPH1 NR3  CPH1    1.0000        0.0000                
    HR3  NR1  CPH1 CPH1    0.5000        0.0000                 
    HR3  NR2  CPH1 CPH1    0.5000        0.0000                 
    N    C    CP1  CP3     0.0000        0.0000                
    NC2  X    X    C      45.0000        0.0000                
    C   HC    HC   NC2      0.0         0.0              
    NC2  X    X    HC      -2.0         0.0              
    NH1  X    X    H      20.0000        0.0000               
    NH2  X    X    H       4.0000        0.0000                
    NR1  CPH1 CPH2 H       0.4500        0.0000                
    NR1  CPH2 CPH1 H       0.4500        0.0000                 
    NR3  CPH1 CPH2 H       1.2000        0.0000      
    NR3  CPH2 CPH1 H       1.2000        0.0000                 
    O    CP1  NH2  CC     45.0000        0.0000                
    O    CT1  NH2  CC     45.0000        0.0000                 
    O    CT2  NH2  CC     45.0000        0.0000                 
    O    CT3  NH2  CC     45.0000        0.0000                 
    O    HA1  NH2  CC     45.0000        0.0000                 
    O    N    CT2  CC    120.0000        0.0000                
    O    NH2  CP1  CC     45.0000        0.0000                 
    O    NH2  CT1  CC     45.0000        0.0000                 
    O    NH2  CT2  CC     45.0000        0.0000                
    O    NH2  CT3  CC     45.0000        0.0000                 
    O    NH2  HA1  CC     45.0000        0.0000                 
    O    X    X    C     120.0000        0.0000                 
    OB   X    X    CD    100.0000        0.0000                 
    OC   X    X    CC     96.0000        0.0000                 
    CC   X    X    CT1    96.0000        0.0000                 
    CC   X    X    CT2    96.0000        0.0000 
    CC   X    X    CT3    96.0000        0.0000", "\n")
                
# Summary
potential_charmm36 = Potential_Charmm36()