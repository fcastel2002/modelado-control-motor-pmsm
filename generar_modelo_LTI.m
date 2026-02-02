%% GENERACIÓN DE MODELO LTI (Linealizado)
%  Proyecto Global Integrador - Automática y Máquinas Eléctricas
%  Autores: Joaquín Calderón - Francisco Castel

clearvars; clc; close all;

%% 1. Cargar Parámetros Físicos
if exist('parametros_sistema_completo.mlx', 'file')
    parametros_sistema_completo;
    disp('-> Parámetros físicos cargados exitosamente.');
else
    error('ERROR CRÍTICO: No se encuentra "parametros_sistema_completo.m".');
end

% Verificación de variables críticas calculadas en el otro script
if ~exist('J_eq', 'var') || ~exist('b_eq', 'var') || ~exist('k_l', 'var') || ~exist('r', 'var')
    error('Faltan variables derivadas (J_eq, b_eq, k_l). Revisa tu archivo de parámetros.');
end

%% 2. DEFINICIÓN SIMBÓLICA DEL SISTEMA (Jacbianos)
disp('-> Calculando Jacobianos Simbólicos...');

syms th_m w_m i_q i_d i_0 T_s real % Estados: x = [th_m; w_m; i_q; i_d; i_0; T_s]
syms v_q v_d v_0 T_ld T_amb real    % Entradas: u = [v_q; v_d; v_0], pert = [T_ld]

% --- Ecuaciones Auxiliares ---
% Resistencia variable con temperatura
Rs_var = R_sREF * (1 + a_cu*(T_s - T_sREF)); 
% Velocidad eléctrica y mecánica
w_r = P_p * w_m; 
th_l = th_m / r;

% --- ECUACIONES DE ESTADO (f) ---
% 1. Posición Mecánica
dot_th_m = w_m;

% 2. Velocidad Mecánica
% Torque Electromagnético
Tm = 1.5 * P_p * (lambda_r * i_q + (L_d - L_q)*i_d*i_q);
% Torque Gravitacional (NO LINEAL)
Tg = (g * k_l / r) * sin(th_l); 
% Ecuación Dinámica
dot_w_m = (1/J_eq) * (Tm - b_eq*w_m - Tg - T_ld/r);

% 3, 4, 5. Corrientes (Eléctrico)
dot_i_q = (1/L_q) * (v_q - Rs_var*i_q - (lambda_r + L_d*i_d)*w_r);
dot_i_d = (1/L_d) * (v_d - Rs_var*i_d + L_q*i_q*w_r);
dot_i_0 = (1/L_ls) * (v_0 - Rs_var*i_0);

% 6. Térmica (Estator)
P_loss = 1.5 * Rs_var * (i_q^2 + i_d^2 + 2*i_0^2);
dot_T_s = (1/C_ts) * (P_loss - (T_s - T_amb)/R_ts_amb);

% Vector de Funciones y Variables
f = [dot_th_m; dot_w_m; dot_i_q; dot_i_d; dot_i_0; dot_T_s];
x_vec = [th_m; w_m; i_q; i_d; i_0; T_s];
u_vec = [v_q; v_d; v_0];

% CÁLCULO DE MATRICES SIMBÓLICAS
A_sym = jacobian(f, x_vec);
B_sym = jacobian(f, u_vec);

%% 3. LINEALIZACIÓN POR CASOS (Evaluación Numérica)

% --- Condiciones Comunes ---
val_w_m = 0;        % Velocidad cero (Equilibrio estático)
val_i_d = 0;        % Control vectorial ideal
val_i_0 = 0;
val_T_s = 20;       % Temp ambiente inicial
val_T_ld = 0;       % Sin perturbación externa

% =========================================================================
% CASO A: BRAZO COLGANDO (theta_l = 0°) - Análisis de Estabilidad Natural
% =========================================================================
disp(' ');
disp('=== CASO A: Linealización en 0 grados (Péndulo Estable) ===');
val_th_l_0 = 0; 
val_th_m_0 = 0;
val_i_q_0 = 0; % No requiere torque para sostenerse

% Sustitución
A_0 = double(subs(A_sym, ...
    [th_m, w_m, i_q, i_d, i_0, T_s, v_q, v_d, v_0, T_ld, T_amb], ...
    [val_th_m_0, val_w_m, val_i_q_0, val_i_d, val_i_0, val_T_s, 0, 0, 0, 0, 20]));

B_0 = double(subs(B_sym, ...
    [th_m, w_m, i_q, i_d, i_0, T_s], ...
    [val_th_m_0, val_w_m, val_i_q_0, val_i_d, val_i_0, val_T_s]));

sys_0 = ss(A_0, B_0, eye(6), 0);

% Mostrar Autovalores Mecánicos (Filas 1 y 2)
disp('Autovalores Mecánicos (0 deg):');
eig_mec_0 = eig(A_0(1:2, 1:2));
disp(eig_mec_0);


% =========================================================================
% CASO B: BRAZO HORIZONTAL (theta_l = 90°) - Peor Caso de Carga
% =========================================================================
disp(' ');
disp('=== CASO B: Linealización en 90 grados (Máxima Carga) ===');
val_th_l_90 = pi/2; 
val_th_m_90 = val_th_l_90 * r;

% Torque requerido para sostener el peso
T_g_90 = g * k_l;       % sin(90) = 1
T_m_req = T_g_90 / r;   % Reflejado al motor
val_i_q_90 = T_m_req / (1.5 * P_p * lambda_r); % Corriente de equilibrio

% Sustitución
A_90 = double(subs(A_sym, ...
    [th_m, w_m, i_q, i_d, i_0, T_s, v_q, v_d, v_0, T_ld, T_amb], ...
    [val_th_m_90, val_w_m, val_i_q_90, val_i_d, val_i_0, val_T_s, 0, 0, 0, 0, 20]));

B_90 = B_0; % La matriz B no depende del estado en este modelo

sys_90 = ss(A_90, B_90, eye(6), 0);

% Mostrar Autovalores Mecánicos
disp('Autovalores Mecánicos (90 deg):');
eig_mec_90 = eig(A_90(1:2, 1:2));
disp(eig_mec_90);
disp('(Nota: Si aparece un 0 o valor muy bajo, indica pérdida de rigidez restauradora)');


% =========================================================================
% CASO C: MODELO DE DISEÑO (Fase 3) - Sin Gravedad
% =========================================================================
disp(' ');
disp('=== CASO C: Modelo de Diseño para Control (Sin Gravedad) ===');
% Eliminamos manualmente el término de gravedad (A21) para el diseño del PID
A_design = A_90; 
A_design(2,1) = 0; % Borramos el acople de posición (Gravedad anulada)

sys_design = ss(A_design, B_90, eye(6), 0);

disp('Modelo sys_design creado. Listo para usar en pidtune o sisotool.');

%% 4. VISUALIZACIÓN
figure('Name', 'Comparación de Polos');
hold on;
h1 = pzmap(sys_0);
h2 = pzmap(sys_90);
legend('0 Grados (Péndulo)', '90 Grados (Masa Pura)');
title('Comparación de Dinámica: 0 vs 90 grados');
grid on;

disp(' ');
disp('--- FINALIZADO ---');
disp('Variables creadas en el Workspace:');
disp('  sys_0      -> Para análisis de estabilidad en vertical');
disp('  sys_90     -> Para análisis en carga máxima');
disp('  sys_design -> Para cálculo de ganancias PID (Fase 3)');