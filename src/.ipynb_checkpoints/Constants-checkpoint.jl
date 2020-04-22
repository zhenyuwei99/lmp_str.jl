function constants()
    res = Dict()
    # Physical Constants
    res["Const_k_b"] = 1.38065e-23;
    res["Const_n_a"] = 6.02214e23;
    res["Const_density_wat"] = 1e3;      # Unit: kg/m^3
    # Energy Converters
    res["Const_kcal2j"] = 4.184e3;
    res["Const_kcalm2j"] = Const_kcal2j/Const_n_a;
    res["Const_kcalm2t"] = Const_kcalm2j/Const_k_b;
    # Mass Converters
    res["Const_g2kg"] = 1e-3;
    res["Const_kg2g"] = 1e3;
    res["Const_gm2g"] = 1/Const_n_a;
    res["Const_gm2kg"] = Const_g2kg/Const_n_a;
    # Length Converters
    res["Const_an2m"] = 1e-10;
    res["Const_an2nm"] = 1e-1;
    res["Const_nm2m"] = 1e-9;
    res["Const_cm2an"] = 1e8;
    res["Const_dm2m"] = 1e-1;
    res["Const_m2dm"] = 1/Const_dm2m;
    # Time Converters
    res["Const_fs2s"] = 1e-15;
    res["Const_ps2s"] = 1e-12;
    res["Const_ns2s"] = 1e-9;
    
    return res
end

# Physical Constants
Const_k_b = 1.38065e-23;
Const_n_a = 6.02214e23;
Const_density_wat = 1e3;      # Unit: kg/m^3
# Energy Converters
Const_kcal2j = 4.184e3;
Const_kcalm2j = Const_kcal2j/Const_n_a;
Const_kcalm2t = Const_kcalm2j/Const_k_b;
# Mass Converters
Const_g2kg = 1e-3;
Const_kg2g = 1e3;
Const_gm2g = 1/Const_n_a;
Const_gm2kg = Const_g2kg/Const_n_a;
# Length Converters
Const_an2m = 1e-10;
Const_an2nm = 1e-1;
Const_nm2m = 1e-9;
Const_cm2an = 1e8;
Const_dm2m = 1e-1;
Const_m2dm = 1/Const_dm2m;
# Time Converters
Const_fs2s = 1e-15;
Const_ps2s = 1e-12;
Const_ns2s = 1e-9;