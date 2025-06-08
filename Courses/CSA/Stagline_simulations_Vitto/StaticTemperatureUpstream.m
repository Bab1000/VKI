function T1 = StaticTemperatureUpstream(H01, M1, Y)

    species_list = fieldnames(Y);
    Y_vals = cellfun(@(s) Y.(s), species_list);

    f = @(T) residual_enthalpy(T, H01, M1, species_list, Y_vals);
    T1 = fzero(f, [200, 5000]); % Kelvin
end

function res = residual_enthalpy(T, H01, M1, species_list, Y_vals)

    R_univ = 8.3145; % J/mol/K

    cp_mass = 0;
     h_mass = 0;
      R_mix = 0;

    for i = 1:length(species_list)
         sp = species_list{i};
         Yi = Y_vals(i);
        M_i = molar_mass(sp);   % kg/mol

        cp_i = JANAF_Cp(T, sp); % J/mol/K
         h_i = JANAF_h(T, sp);  % J/mol

        cp_mass = cp_mass + Yi * cp_i / M_i;   % J/kg/K
         h_mass = h_mass  + Yi * h_i / M_i;    % J/kg
          R_mix = R_mix   + Yi * R_univ / M_i; % J/kg/K
    end

    % gamma = cp_mass / (cp_mass - R_mix);
       Ma_vals = [5.487, 6.402];
    gamma_vals = [1.336, 1.314];
     Ma_interp = 6;
    gamma = interp1(Ma_vals, gamma_vals, Ma_interp);
    fprintf('Interpolated gamma at Ma = %.3f is %.4f\n', Ma_interp, gamma);

    H1 = h_mass;
    V1_squared = gamma * R_mix * T * M1^2;
    H01_computed = H1 + 0.5 * V1_squared;

    res = H01 - H01_computed;
end

function M = molar_mass(species)
    switch species
        case 'N2', M = 28.0134e-3;
        case 'O2', M = 31.9988e-3;
        case 'N',  M = 14.0067e-3;
        case 'O',  M = 15.9994e-3;
        case 'NO', M = 30.0061e-3;
        otherwise, error('Unknown species');
    end
end
