%% CÁLCULO DE INCERTIDUMBRE Y SINTONÍA DE CONTROL
%  Este script es llamado por 'parametros_sistema_completo'.
%  Modifica las variables del workspace según el escenario y recalcula PIDs.

% --- A. GESTIÓN DE ESCENARIOS Y TOLERANCIAS ---
if ~exist('ESCENARIO', 'var'); ESCENARIO = 'NOMINAL'; end

% Definición de Factores de Tolerancia
tol_pct = 0.01;   % 1% Multiplicativo
tol_bl  = 0.03;   % Absoluto para fricción
tol_Tld = 5.0;    % Absoluto para perturbación

switch ESCENARIO
    case 'WORST'
        disp('   [LÓGICA] Aplicando tolerancias de PEOR CASO...');
        % Mecánica (Más pesado y duro)
        m_l = 1.5;                  
        b_l = b_l + tol_bl;        
        b_m = b_m * (1 + tol_pct);  
        J_m = J_m * (1 + tol_pct);  
        
        % Eléctrica (Más impedancia)
        R_sREF   = R_sREF * (1 + tol_pct);
        L_q      = L_q    * (1 + tol_pct);
        L_d      = L_d    * (1 + tol_pct);
        L_ls     = L_ls   * (1 + tol_pct);
        lambda_r = lambda_r * (1 - tol_pct); % Menos fuerza
        
        % Térmica
        C_ts     = C_ts * (1 - tol_pct);
        R_ts_amb = R_ts_amb * (1 + tol_pct);
        T_ld     = tol_Tld;

    case 'BEST'
        disp('   [LÓGICA] Aplicando tolerancias de MEJOR CASO...');
        % Mecánica (Liviano)
        m_l = 0;
        b_l = b_l - tol_bl;
        b_m = b_m * (1 - tol_pct);
        J_m = J_m * (1 - tol_pct);
        
        % Eléctrica (Ideal)
        R_sREF   = R_sREF * (1 - tol_pct);
        L_q      = L_q    * (1 - tol_pct);
        L_d      = L_d    * (1 - tol_pct);
        L_ls     = L_ls   * (1 - tol_pct);
        lambda_r = lambda_r * (1 + tol_pct);
        
        % Térmica
        C_ts     = C_ts * (1 + tol_pct);
        R_ts_amb = R_ts_amb * (1 - tol_pct);
        T_ld     = 0;

    otherwise
        disp('   [LÓGICA] Usando valores NOMINALES.');
end

% --- B. RE-CÁLCULO DE VARIABLES DERIVADAS (Obligatorio aquí) ---
b_eq = b_m + b_l/(r^2);                 
J_l  = (m*l_cm^2 + J_cm) + m_l*l_l^2;   
J_eq = J_m + J_l/r^2;                   
k_l  = m*l_cm + m_l*l_l;                

% --- C. ACTUALIZACIÓN DEL STRUCT PARA SIMULINK ---
params.Ld = L_d; params.Lq = L_q; params.Lls = L_ls;
params.Rs_ref = R_sREF; params.lambda_m = lambda_r;
params.J_total = J_eq; params.b_total = b_eq;
params.kl = k_l; params.Cts = C_ts; params.Rts_amb = R_ts_amb;

% --- D. AUTO-TUNING DE CONTROLADORES (Basado en planta actual) ---
% Anchos de banda
alpha_c = 2 * pi * 500; % 500 Hz (Corriente)
alpha_w = 2 * pi * 50;  % 50 Hz (Velocidad)

% PI Corriente
Kp_iq = alpha_c * L_q;    Ki_iq = alpha_c * R_sREF;
Kp_id = alpha_c * L_d;    Ki_id = alpha_c * R_sREF;

% PI Velocidad
Kp_w  = alpha_w * J_eq;   Ki_w  = alpha_w * b_eq;

% Ganancia Feedforward
Gain_inv_kt = 1 / (1.5 * P_p * lambda_r);

fprintf('   [CONTROL] Sintonizado para J_eq = %.2e, m_l = %.2f kg\n', J_eq, m_l);