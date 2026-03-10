function mdl = modelo_lti_5_1_revisado(p)
%MODELO_LTI_5_1_REVISADO Modelo LTI equivalente y aumentado de 5.1.
%   mdl = MODELO_LTI_5_1_REVISADO(p) devuelve matrices y funciones de
%   transferencia del modelo linealizado por retroalimentacion NL.

if nargin < 1 || isempty(p)
    p = parametros_revisados_5_1();
end

% Modelo nominal de 3 estados: x = [theta_m; omega_m; i_qs]
mdl.A_nom = [
    0,               1,                 0;
    0,      -p.b_eq/p.J_eq,     p.k_t/p.J_eq;
    0, -p.P_p*p.lambda_m/p.L_q, -p.R_s/p.L_q
];

mdl.B_vqs_nom = [0; 0; 1/p.L_q];
mdl.B_tld_nom = [0; -1/(p.r*p.J_eq); 0];
mdl.C_theta = [1, 0, 0];
mdl.D_zero = 0;

% Modelo aumentado de 4 estados: agrega i_ds residual desacoplada
mdl.A_aug = blkdiag(mdl.A_nom, -p.R_s/p.L_d);
mdl.B_vqs_aug = [mdl.B_vqs_nom; 0];
mdl.B_tld_aug = [mdl.B_tld_nom; 0];
mdl.B_vds_aug = [0; 0; 0; 1/p.L_d];
mdl.C_theta_aug = [1, 0, 0, 0];

% Denominador comun de las funciones de transferencia
mdl.den = [
    1, ...
    p.b_eq/p.J_eq + p.R_s/p.L_q, ...
    p.b_eq*p.R_s/(p.J_eq*p.L_q) + 3*p.P_p^2*p.lambda_m^2/(2*p.J_eq*p.L_q), ...
    0
];

mdl.num_vqs = [3*p.P_p*p.lambda_m/(2*p.J_eq*p.L_q)];
mdl.num_tld = [-1/(p.r*p.J_eq), -p.R_s/(p.r*p.J_eq*p.L_q)];
mdl.polo_ids = -p.R_s/p.L_d;
mdl.cero_tld = -p.R_s/p.L_q;
end
