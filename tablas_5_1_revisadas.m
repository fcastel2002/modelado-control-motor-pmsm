function tablas = tablas_5_1_revisadas()
%TABLAS_5_1_REVISADAS Consolida valores y tablas numericas de 5.1.
%   El objetivo es evitar numeros tipeados a mano en el informe.

p = parametros_revisados_5_1();
mdl = modelo_lti_5_1_revisado(p);

polos_nom = eig(mdl.A_nom);
polos_aug = eig(mdl.A_aug);

tablas.parametros_nominales = table( ...
    ["J_l"; "k_l"; "J_eq"; "b_eq"; "k_t"; "tau_q"; "tau_d"; "tau_th"], ...
    [p.J_l; p.k_l; p.J_eq; p.b_eq; p.k_t; p.tau_q; p.tau_d; p.tau_th], ...
    ["kg.m^2"; "kg.m"; "kg.m^2"; "N.m.s/rad"; "N.m/A"; "s"; "s"; "s"], ...
    'VariableNames', {'Simbolo', 'Valor', 'Unidad'});

tablas.autovalores = table( ...
    (1:numel(polos_aug)).', real(polos_aug), imag(polos_aug), ...
    'VariableNames', {'Indice', 'ParteReal', 'ParteImaginaria'});

Rs_vec = linspace(0.8 * p.R_s_ref, 1.2 * p.R_s_ref, 13).';
zeta_Rs = zeros(size(Rs_vec));
wn_Rs = zeros(size(Rs_vec));
for k = 1:numel(Rs_vec)
    pk = parametros_revisados_5_1(0.0, 0.1, 20.0, 20.0);
    pk.R_s = Rs_vec(k);
    Ak = [
        0,               1,                 0;
        0,      -pk.b_eq/pk.J_eq,     pk.k_t/pk.J_eq;
        0, -pk.P_p*pk.lambda_m/pk.L_q, -pk.R_s/pk.L_q
    ];
    eig_k = eig(Ak);
    eig_em = eig_k(abs(imag(eig_k)) > 0);
    if ~isempty(eig_em)
        wn_Rs(k) = abs(eig_em(1));
        zeta_Rs(k) = -real(eig_em(1)) / wn_Rs(k);
    end
end

tablas.migracion_Rs = table(Rs_vec, wn_Rs, zeta_Rs, -Rs_vec / p.L_q, ...
    'VariableNames', {'Rs', 'omega_n', 'zeta', 'cero_G2'});

fprintf('\n=== Parametros nominales revisados 5.1 ===\n');
disp(tablas.parametros_nominales);
fprintf('\n=== Autovalores del modelo aumentado ===\n');
disp(tablas.autovalores);
fprintf('\n=== Migracion con Rs ===\n');
disp(tablas.migracion_Rs);
fprintf('Polo residual i_ds = %.6f rad/s\n', mdl.polo_ids);
fprintf('Cero de G2 = %.6f rad/s\n', mdl.cero_tld);
fprintf('Polos nominales activos = \n');
disp(polos_nom);
end
