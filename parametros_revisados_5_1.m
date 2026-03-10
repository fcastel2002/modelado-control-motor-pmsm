function p = parametros_revisados_5_1(m_l, b_l, T_s, T_amb)
%PARAMETROS_REVISADOS_5_1 Parametros nominales y derivados usados en 5.1.
%   p = PARAMETROS_REVISADOS_5_1() devuelve el conjunto nominal.
%   p = PARAMETROS_REVISADOS_5_1(m_l, b_l, T_s, T_amb) permite evaluar
%   escenarios de carga, friccion y temperatura sin pisar los scripts
%   originales del proyecto.

if nargin < 1 || isempty(m_l),  m_l = 0.0;  end
if nargin < 2 || isempty(b_l),  b_l = 0.1;  end
if nargin < 3 || isempty(T_s),  T_s = 20.0; end
if nargin < 4 || isempty(T_amb), T_amb = 20.0; end

% Datos dados por consigna
p.r = 120.0;
p.g = 9.80665;
p.m = 1.0;
p.l_cm = 0.25;
p.J_cm = 0.0208;
p.l_l = 0.50;
p.m_l = m_l;
p.b_l = b_l;
p.J_m = 14.0e-6;
p.b_m = 15.0e-6;
p.P_p = 3;
p.lambda_m = 0.016;
p.L_q = 5.8e-3;
p.L_d = 6.6e-3;
p.L_ls = 0.8e-3;
p.R_s_ref = 1.02;
p.T_s_ref = 20.0;
p.alpha_cu = 3.9e-3;
p.C_ts = 0.818;
p.R_ts_amb = 146.7;
p.T_s = T_s;
p.T_amb = T_amb;

% Parametros derivados
p.J_l = (p.m * p.l_cm^2 + p.J_cm) + p.m_l * p.l_l^2;
p.k_l = p.m * p.l_cm + p.m_l * p.l_l;
p.J_eq = p.J_m + p.J_l / p.r^2;
p.b_eq = p.b_m + p.b_l / p.r^2;
p.R_s = p.R_s_ref * (1 + p.alpha_cu * (p.T_s - p.T_s_ref));
p.k_t = 1.5 * p.P_p * p.lambda_m;
p.tau_q = p.L_q / p.R_s;
p.tau_d = p.L_d / p.R_s;
p.tau_th = p.R_ts_amb * p.C_ts;
p.tau_m = p.J_eq / p.b_eq;

% Limites de referencia de la consigna
p.n_m_nom_rpm = 6600.0;
p.omega_m_nom = p.n_m_nom_rpm * 2 * pi / 60;
p.n_l_nom_rpm = 60.0;
p.omega_l_nom = p.n_l_nom_rpm * 2 * pi / 60;
p.T_q_nom = 17.0;
p.T_q_max = 45.0;
p.T_ld_max = 5.0;

assert(p.J_eq > 0, "J_eq debe ser positivo.");
assert(p.b_eq > 0, "b_eq debe ser positivo.");
assert(p.R_s > 0, "R_s debe ser positiva.");
end
