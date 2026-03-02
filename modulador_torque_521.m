%% DISEÑO DEL MODULADOR DE TORQUE - CONSIGNA 5.2.1
%  Control vectorial interno de corriente/torque
%  Sub-ítems: a) Desacoplamiento, b) Control P de corriente,
%             c) Consigna de torque, d) Compensación de gravedad
%  Autores: Joaquín Calderón - Francisco Castel

clearvars; clc; close all;

%% 1. CARGAR PARÁMETROS FÍSICOS
parametros_sistema_completo;
R_s0 = R_sREF;
kt   = 1.5 * P_p * lambda_r;  % Constante de torque [N.m/A]

fprintf('=== DISEÑO DEL MODULADOR DE TORQUE (5.2.1) ===\n');
fprintf('Parámetros del motor:\n');
fprintf('  R_s0 = %.4f Ohm,  L_q = %.4f H,  L_d = %.4f H\n', R_s0, L_q, L_d);
fprintf('  kt   = %.6f N.m/A,  P_p = %d\n', kt, P_p);

%% =====================================================================
%  5.2.1(a): DESACOPLAMIENTO COMPLETO
%  =====================================================================
% Ecuaciones originales del estator (en qd0):
%   Lq*di_qs/dt = v_qs - Rs*i_qs - (lambda_m + Ld*i_ds)*omega_r
%   Ld*di_ds/dt = v_ds - Rs*i_ds + Lq*i_qs*omega_r
%
% Ley de control con desacoplamiento COMPLETO:
%   v_qs = v_qs' + Rs_hat*i_qs + (lambda_m + Ld*i_ds)*omega_r
%   v_ds = v_ds' + Rs_hat*i_ds - Lq*i_qs*omega_r
%
% Planta desacoplada (si Rs_hat = Rs):
%   Lq*di_qs/dt = v_qs'   -->  G_iq(s) = 1/(Lq*s)  (integrador puro)
%   Ld*di_ds/dt = v_ds'   -->  G_id(s) = 1/(Ld*s)  (integrador puro)
%
% Comparación con 5.1.2.c.VI (ley mínima + complementaria):
%   - En 5.1.2.c.VI solo se compensaron los ACOPLAMIENTOS CRUZADOS
%     (términos con omega_r), dejando Rs como parte de la planta.
%   - Aquí se compensa TAMBIÉN Rs, dejando integradores puros.
%   - Esto permite usar control P solamente (sin I) con SS error = 0.

fprintf('\n--- 5.2.1(a): Desacoplamiento Completo ---\n');
fprintf('Compensaciones aplicadas:\n');
fprintf('  Eje q: +Rs*i_qs + (lambda_m + Ld*i_ds)*omega_r  (Rs + back-EMF)\n');
fprintf('  Eje d: +Rs*i_ds - Lq*i_qs*omega_r               (Rs + cross-coupling)\n');
fprintf('Planta resultante: integradores puros 1/(L*s)\n');

%% =====================================================================
%  5.2.1(b): CONTROL PROPORCIONAL DE CORRIENTE
%  =====================================================================
% Polo deseado: p_i = -5000 rad/s  (BW ≈ 796 Hz)
alpha_i = 5000;  % [rad/s]

% Ganancias proporcionales (planta desacoplada = 1/(L*s))
% CL: Kp/(L*s) / (1 + Kp/(L*s)) = (Kp/L) / (s + Kp/L)
% Polo CL: -Kp/L = -alpha_i  =>  Kp = alpha_i * L
Kp_iq = alpha_i * L_q;
Kp_id = alpha_i * L_d;

fprintf('\n--- 5.2.1(b): Control Proporcional de Corriente ---\n');
fprintf('Polo deseado: p_i = -%.0f rad/s (BW = %.0f Hz)\n', alpha_i, alpha_i/(2*pi));
fprintf('Ganancias:\n');
fprintf('  Kp_iq = alpha_i * Lq = %.0f * %.4f = %.4f\n', alpha_i, L_q, Kp_iq);
fprintf('  Kp_id = alpha_i * Ld = %.0f * %.4f = %.4f\n', alpha_i, L_d, Kp_id);

% Funciones de transferencia de lazo cerrado
s = tf('s');
G_iq_OL = 1/(L_q*s);          % Planta desacoplada eje q
G_id_OL = 1/(L_d*s);          % Planta desacoplada eje d
G_iq_CL = feedback(Kp_iq * G_iq_OL, 1);  % = alpha_i/(s + alpha_i)
G_id_CL = feedback(Kp_id * G_id_OL, 1);

fprintf('\nFT lazo cerrado de corriente:\n');
fprintf('  G_iq_CL(s) = %.0f / (s + %.0f)    tau = %.4f ms\n', alpha_i, alpha_i, 1/alpha_i*1000);
fprintf('  G_id_CL(s) = %.0f / (s + %.0f)    tau = %.4f ms\n', alpha_i, alpha_i, 1/alpha_i*1000);

% Comparación con polos originales de la planta (5.1.6 lazo abierto)
polo_orig_q = -R_s0/L_q;
polo_orig_d = -R_s0/L_d;
fprintf('\nComparación con planta original (lazo abierto):\n');
fprintf('  Polo original eje q: %.1f rad/s (tau = %.2f ms)\n', polo_orig_q, -1/polo_orig_q*1000);
fprintf('  Polo original eje d: %.1f rad/s (tau = %.2f ms)\n', polo_orig_d, -1/polo_orig_d*1000);
fprintf('  Polo con modulador:  %.0f rad/s (tau = %.2f ms)\n', -alpha_i, 1/alpha_i*1000);
fprintf('  Factor de aceleración: %.0fx (eje q), %.0fx (eje d)\n', ...
    abs(alpha_i/polo_orig_q), abs(alpha_i/polo_orig_d));

%% --- Investigación: ¿Es necesaria acción integral? ---
fprintf('\n--- Investigación: Acción Integral ---\n');

% CASO 1: Con Rs compensado perfectamente (planta = 1/(Ls))
% P control + integrador en planta => SS error = 0 para escalón de i*
fprintf('CASO 1 (Rs compensado): Planta = 1/(Ls) [integrador]\n');
fprintf('  -> P-only da error estacionario = 0 (el integrador lo asegura)\n');
fprintf('  -> NO se necesita acción integral.\n');

% CASO 2: Sin compensar Rs (planta = 1/(Ls + Rs))
% P control sin integrador en planta => SS error ≠ 0
G_iq_nocomp = 1/(L_q*s + R_s0);  % Planta SIN compensar Rs
G_iq_CL_nocomp = feedback(Kp_iq * G_iq_nocomp, 1);
ss_error_pct = R_s0 / (R_s0 + Kp_iq) * 100;

fprintf('CASO 2 (Rs NO compensado): Planta = 1/(Ls + Rs)\n');
fprintf('  -> P-only da error estacionario = Rs/(Rs+Kp) = %.2f%%\n', ss_error_pct);
fprintf('  -> Se necesitaría acción integral para anular este error.\n');
fprintf('  -> Polo CL se desplaza a: %.1f rad/s (aún aceptable)\n', ...
    -(R_s0 + Kp_iq)/L_q);

fprintf('\nConclusión: Con desacoplamiento completo de Rs, P-only es suficiente.\n');
fprintf('Si Rs se compensa con valor nominal (no adaptativo), el error\n');
fprintf('residual es despreciable (Delta_Rs << Kp).\n');

%% =====================================================================
%  5.2.1(c): CONSIGNA DE TORQUE Y COMPENSACIÓN DE FRICCIÓN
%  =====================================================================
% El modulador recibe T* y produce:
%   i_qs* = T* / kt                           (conversión torque->corriente)
%   i_ds* = 0                                  (estrategia base)
%
% Compensación de fricción viscosa:
%   T*_total = T*_ext + b_eq * omega_m
%   i_qs* = T*_total / kt
%
% Efecto: elimina el término b_eq*omega_m de la dinámica mecánica:
%   J_eq * dw/dt = T*_ext - T_grav - T_ld/r

fprintf('\n--- 5.2.1(c): Consigna de Torque ---\n');
fprintf('Conversión: i_qs* = T* / kt = T* / %.6f\n', kt);
fprintf('Compensación de fricción: T*_total = T*_ext + b_eq * omega_m\n');
fprintf('  b_eq = %.6e N.m/(rad/s)\n', b_eq);

%% =====================================================================
%  5.2.1(d): COMPENSACIÓN DE GRAVEDAD (Feedforward)
%  =====================================================================
% T_grav = (g * k_l / r) * sin(theta_m / r)
%
% Compensación: agregar T_grav al comando de torque:
%   T*_total = T*_ext + b_eq*omega_m + (g*k_l/r)*sin(theta_m/r)
%
% Resultado: planta mecánica vista por el controlador externo:
%   J_eq * dw/dt = T*_ext - T_ld/r    (doble integrador puro + perturbación)

T_grav_max = g * k_l / r;  % Torque grav máximo (sin(90°)=1, m_l=0)
fprintf('\n--- 5.2.1(d): Compensación de Gravedad ---\n');
fprintf('Torque gravitacional máximo (m_l=0, theta_l=90°): %.4f N.m (ref. motor)\n', T_grav_max);
fprintf('Corriente i_qs equivalente: %.4f A\n', T_grav_max/kt);
fprintf('\nPlanta vista por el controlador externo (después de todas las compensaciones):\n');
fprintf('  theta_m / T*_ext = 1 / (J_eq * s^2)    [doble integrador]\n');
fprintf('  J_eq = %.6e kg.m2\n', J_eq);

%% =====================================================================
%  5.2.1(a): COMPARACIÓN Rs(Ts) REAL vs Rs NOMINAL
%  =====================================================================
fprintf('\n=== EVALUACIÓN: Rs(Ts) REAL vs Rs NOMINAL ===\n');

% Escenario: T_amb = 40°C (peor caso térmico según consigna 4.c)
T_amb_test = 40;  % °C
Rs_40 = R_sREF * (1 + a_cu * (T_amb_test - T_sREF));
Delta_Rs = Rs_40 - R_s0;

fprintf('T_amb = %d°C:  Rs(40°C) = %.4f Ohm  (nominal: %.4f Ohm)\n', T_amb_test, Rs_40, R_s0);
fprintf('Delta_Rs = %.4f Ohm (%.1f%%)\n', Delta_Rs, Delta_Rs/R_s0*100);
fprintf('\nCon P-only y desacoplamiento:\n');
fprintf('  Si Rs_hat = Rs(Ts): cancelación perfecta, polo CL = -%.0f rad/s\n', alpha_i);
fprintf('  Si Rs_hat = Rs_nom: residuo Delta_Rs*i en planta\n');
fprintf('    Polo CL efectivo ≈ -(%.0f + %.1f/%.4f) ≈ -%.0f rad/s (cambio < 0.3%%)\n', ...
    alpha_i, Delta_Rs, L_q, alpha_i + Delta_Rs/L_q);
fprintf('  => Efecto despreciable. Rs nominal es aceptable para P control.\n');

%% =====================================================================
%  SIMULACIÓN: MODULADOR DE TORQUE vs LAZO ABIERTO (5.1.6)
%  =====================================================================
fprintf('\n=== SIMULACIÓN COMPARATIVA ===\n');

% Señales de entrada (mismas que 5.1.6)
t_final = 2.2;
dt = 5e-5;
t_sim = (0:dt:t_final)';
N = length(t_sim);

V_QS_AMP = 19.596;
T_LD_AMP = 6.28;

% Señal de torque de carga (doble pulso)
tld_func = @(t) T_LD_AMP * ( ...
    (t >= 0.3 & t < 0.5) - (t >= 0.5 & t < 0.9) + ...
    (t >= 1.3 & t < 1.5) - (t >= 1.5 & t < 1.9));

% Para comparar con 5.1.6: usamos consigna de torque equivalente
% En 5.1.6, v_qs* = 19.596V producía i_qs_ss = V/Rs = 19.21A -> T = kt*19.21
i_qs_ref_amp = V_QS_AMP / R_s0;   % [A] Corriente SS equivalente
T_ref_amp = kt * i_qs_ref_amp;     % [N.m] Torque equivalente

% Consigna de torque externo (misma forma que v_qs* en 5.1.6)
Tref_func = @(t) T_ref_amp * ((t >= 0.1 & t < 0.7) - (t >= 1.1 & t < 1.7));

fprintf('Consigna de torque equivalente: T* = %.4f N.m (i_qs* = %.2f A)\n', T_ref_amp, i_qs_ref_amp);

% Opciones del integrador
ode_opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 1e-3);

% --- Caso 1: Modulador con Rs(Ts) adaptativo ---
fprintf('\nCaso 1: Rs_hat = Rs(Ts) adaptativo...\n');
x0 = [0; 0; 0; 0; 0; T_amb];
[~, x_adapt] = ode15s(@(t,x) pmsm_modulador_ode(t, x, ...
    Tref_func, tld_func, ...
    Kp_iq, Kp_id, kt, true, ...  % Rs adaptativo = true
    L_q, L_d, L_ls, R_sREF, a_cu, T_sREF, lambda_r, P_p, ...
    J_eq, b_eq, r, k_l, g, C_ts, R_ts_amb, T_amb), ...
    [0, t_final], x0, ode_opts);

% Interpolar
x_adapt = interp1(linspace(0, t_final, size(x_adapt,1)), x_adapt, t_sim, 'pchip');
fprintf('  Completado.\n');

% --- Caso 2: Modulador con Rs nominal (constante) ---
fprintf('Caso 2: Rs_hat = Rs_nominal constante...\n');
[~, x_nom] = ode15s(@(t,x) pmsm_modulador_ode(t, x, ...
    Tref_func, tld_func, ...
    Kp_iq, Kp_id, kt, false, ...  % Rs adaptativo = false
    L_q, L_d, L_ls, R_sREF, a_cu, T_sREF, lambda_r, P_p, ...
    J_eq, b_eq, r, k_l, g, C_ts, R_ts_amb, T_amb), ...
    [0, t_final], x0, ode_opts);

x_nom = interp1(linspace(0, t_final, size(x_nom,1)), x_nom, t_sim, 'pchip');
fprintf('  Completado.\n');

% --- Caso 3: Lazo abierto (reproducir 5.1.6 para comparar) ---
fprintf('Caso 3: Lazo abierto (5.1.6, v_qs directo)...\n');
vqs_func = @(t) V_QS_AMP * ((t >= 0.1 & t < 0.7) - (t >= 1.1 & t < 1.7));

[~, x_OL] = ode15s(@(t,x) pmsm_lazo_abierto_ode(t, x, ...
    vqs_func, tld_func, ...
    L_q, L_d, L_ls, R_sREF, a_cu, T_sREF, lambda_r, P_p, ...
    J_eq, b_eq, r, k_l, g, C_ts, R_ts_amb, T_amb), ...
    [0, t_final], x0, ode_opts);

x_OL = interp1(linspace(0, t_final, size(x_OL,1)), x_OL, t_sim, 'pchip');
fprintf('  Completado.\n');

% Pre-computar señales de entrada como vectores
u_tld = tld_func(t_sim);
u_Tref = Tref_func(t_sim);

%% =====================================================================
%  FIGURAS
%  =====================================================================

% --- Figura 1: Corriente i_qs (velocidad de respuesta) ---
fig1 = figure('Name', 'Corriente i_qs: Modulador vs Lazo Abierto', ...
    'Color', 'w', 'Position', [50 100 1100 500]);

subplot(2,1,1); hold on; grid on;
plot(t_sim, x_adapt(:,3), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Modulador (Rs adaptativo)');
plot(t_sim, x_nom(:,3), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Modulador (Rs nominal)');
plot(t_sim, x_OL(:,3), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Lazo abierto (5.1.6)');
yline(i_qs_ref_amp, 'g--', 'LineWidth', 0.8, 'DisplayName', 'i_{qs}^* referencia');
ylabel('i_{qs} [A]');
title('Corriente de Torque i_{qs}(t)', 'FontSize', 12);
legend('Location', 'best'); xlim([0 t_final]);

subplot(2,1,2); hold on; grid on;
% Zoom al primer escalón (t = 0.1s)
plot(t_sim, x_adapt(:,3), 'b-', 'LineWidth', 1.5);
plot(t_sim, x_nom(:,3), 'r--', 'LineWidth', 1.5);
plot(t_sim, x_OL(:,3), 'k:', 'LineWidth', 1.5);
yline(i_qs_ref_amp, 'g--', 'LineWidth', 0.8);
ylabel('i_{qs} [A]'); xlabel('Tiempo [s]');
title('Zoom al primer transitorio', 'FontSize', 12);
xlim([0.095 0.115]); % Zoom 20ms alrededor del escalón

sgtitle('Comparación: Modulador de Torque vs Lazo Abierto', ...
    'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 2: Estados mecánicos ---
fig2 = figure('Name', 'Estados Mecánicos: Modulador vs LA', ...
    'Color', 'w', 'Position', [100 50 1100 600]);

subplot(2,1,1); hold on; grid on;
plot(t_sim, x_adapt(:,1), 'b-', 'LineWidth', 1.2);
plot(t_sim, x_nom(:,1), 'r--', 'LineWidth', 1.2);
plot(t_sim, x_OL(:,1), 'k:', 'LineWidth', 1.2);
ylabel('\theta_m [rad]');
title('Posición Angular', 'FontSize', 12);
legend('Modulador (Rs adapt.)', 'Modulador (Rs nom.)', 'Lazo abierto', 'Location', 'best');
xlim([0 t_final]);

subplot(2,1,2); hold on; grid on;
plot(t_sim, x_adapt(:,2), 'b-', 'LineWidth', 1.2);
plot(t_sim, x_nom(:,2), 'r--', 'LineWidth', 1.2);
plot(t_sim, x_OL(:,2), 'k:', 'LineWidth', 1.2);
ylabel('\omega_m [rad/s]'); xlabel('Tiempo [s]');
title('Velocidad Angular', 'FontSize', 12);
legend('Modulador (Rs adapt.)', 'Modulador (Rs nom.)', 'Lazo abierto', 'Location', 'best');
xlim([0 t_final]);

sgtitle('Estados Mecánicos: Modulador de Torque vs Lazo Abierto', ...
    'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 3: Corriente i_ds ---
fig3 = figure('Name', 'Corriente i_ds', 'Color', 'w', 'Position', [150 100 1000 400]);
hold on; grid on;
plot(t_sim, x_adapt(:,4), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Modulador (Rs adapt.)');
plot(t_sim, x_nom(:,4), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Modulador (Rs nom.)');
plot(t_sim, x_OL(:,4), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Lazo abierto');
ylabel('i_{ds} [A]'); xlabel('Tiempo [s]');
title('Corriente de Campo i_{ds}(t) - Con i_{ds}^* = 0', 'FontSize', 12);
legend('Location', 'best'); xlim([0 t_final]);

% --- Figura 4: Comparación Rs(Ts) vs Rs nominal (detalle) ---
fig4 = figure('Name', 'Efecto Rs(Ts) vs Rs nominal', ...
    'Color', 'w', 'Position', [200 50 1100 700]);

subplot(2,2,1); hold on; grid on;
plot(t_sim, x_adapt(:,3), 'b-', 'LineWidth', 1.2);
plot(t_sim, x_nom(:,3), 'r--', 'LineWidth', 1.2);
ylabel('i_{qs} [A]'); title('Corriente i_{qs}', 'FontSize', 11);
legend('Rs(Ts)', 'Rs nom.', 'Location', 'best'); xlim([0 t_final]);

subplot(2,2,2); hold on; grid on;
plot(t_sim, x_adapt(:,2), 'b-', 'LineWidth', 1.2);
plot(t_sim, x_nom(:,2), 'r--', 'LineWidth', 1.2);
ylabel('\omega_m [rad/s]'); title('Velocidad \omega_m', 'FontSize', 11);
xlim([0 t_final]);

subplot(2,2,3); hold on; grid on;
plot(t_sim, x_adapt(:,6), 'b-', 'LineWidth', 1.2);
plot(t_sim, x_nom(:,6), 'r--', 'LineWidth', 1.2);
ylabel('T_s [°C]'); xlabel('Tiempo [s]');
title('Temperatura Estator', 'FontSize', 11); xlim([0 t_final]);

subplot(2,2,4); hold on; grid on;
% Error de corriente
error_adapt = i_qs_ref_amp * ((t_sim >= 0.1 & t_sim < 0.7) - (t_sim >= 1.1 & t_sim < 1.7)) - x_adapt(:,3);
error_nom   = i_qs_ref_amp * ((t_sim >= 0.1 & t_sim < 0.7) - (t_sim >= 1.1 & t_sim < 1.7)) - x_nom(:,3);
plot(t_sim, error_adapt, 'b-', 'LineWidth', 1.2);
plot(t_sim, error_nom, 'r--', 'LineWidth', 1.2);
ylabel('Error i_{qs} [A]'); xlabel('Tiempo [s]');
title('Error de Seguimiento de Corriente', 'FontSize', 11);
legend('Rs(Ts)', 'Rs nom.', 'Location', 'best'); xlim([0 t_final]);

sgtitle('Efecto de Compensación Rs(T_s) Adaptativa vs Nominal', ...
    'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 5: Diagrama de polos - Comparación ---
fig5 = figure('Name', 'Mapa de Polos: LA vs Modulador', ...
    'Color', 'w', 'Position', [250 100 800 500]);
hold on; grid on;

% Polos lazo abierto (LTI aumentado)
A_aug = [0, 1, 0, 0;
         0, -b_eq/J_eq, kt/J_eq, 0;
         0, 0, -R_s0/L_q, 0;
         0, 0, 0, -R_s0/L_d];
polos_LA = eig(A_aug);
plot(real(polos_LA), imag(polos_LA), 'rx', 'MarkerSize', 12, 'LineWidth', 2, ...
    'DisplayName', 'Lazo abierto (5.1)');

% Polos con modulador (corrientes cerradas a -5000)
% El sistema efectivo tiene: theta_m, omega_m con planta simplificada
% más los polos de corriente cerrados
polos_mod = [0; -b_eq/J_eq; -alpha_i; -alpha_i];  % Aproximación
plot(real(polos_mod), imag(polos_mod), 'bo', 'MarkerSize', 10, 'LineWidth', 2, ...
    'DisplayName', 'Con modulador (5.2.1)');

xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');
title('Migración de Polos: Lazo Abierto vs Con Modulador de Torque', 'FontSize', 12);
legend('Location', 'northwest');
xline(0, 'k--'); yline(0, 'k--');

%% RESUMEN FINAL
fprintf('\n========================================\n');
fprintf('  FIGURAS GENERADAS (guardar manualmente)\n');
fprintf('========================================\n');
fprintf('  Fig 1: corriente_iq_modulador_vs_LA.png\n');
fprintf('  Fig 2: estados_mec_modulador_vs_LA.png\n');
fprintf('  Fig 3: corriente_id_modulador.png\n');
fprintf('  Fig 4: efecto_Rs_adaptativo_vs_nominal.png\n');
fprintf('  Fig 5: polos_LA_vs_modulador.png\n');
fprintf('========================================\n');

%% =====================================================================
%  FUNCIONES LOCALES
%  =====================================================================

function dxdt = pmsm_modulador_ode(t, x, Tref_func, tld_func, ...
    Kp_iq, Kp_id, kt, Rs_adaptativo, ...
    Lq, Ld, Lls, Rs_ref, alpha_cu, Ts_ref, lambda_m, Pp, ...
    J_eq, b_eq, r_red, kl, grav, Cts, Rts_amb, Tamb)
% PMSM con MODULADOR DE TORQUE en lazo cerrado
% x = [theta_m; omega_m; i_qs; i_ds; i_0s; T_s]
% Incluye: desacoplamiento completo + control P + compensación fricción/gravedad

    theta_m = x(1); omega_m = x(2);
    i_qs = x(3); i_ds = x(4); i_0s = x(5); T_s = x(6);

    omega_r = Pp * omega_m;
    theta_l = theta_m / r_red;

    % Rs real (planta siempre tiene Rs variable)
    Rs_real = Rs_ref * (1 + alpha_cu * (T_s - Ts_ref));

    % Rs estimado (para compensación)
    if Rs_adaptativo
        Rs_hat = Rs_real;  % Usa sensor de temperatura
    else
        Rs_hat = Rs_ref;   % Usa valor nominal constante
    end

    % --- MODULADOR DE TORQUE ---
    % 1. Consigna de torque externo
    T_ref = Tref_func(t);
    T_ld  = tld_func(t);

    % 2. Compensación de fricción viscosa
    T_friccion = b_eq * omega_m;

    % 3. Compensación de gravedad (feedforward)
    T_grav_comp = (grav * kl / r_red) * sin(theta_l);

    % 4. Torque total comandado
    T_total = T_ref + T_friccion + T_grav_comp;

    % 5. Conversión a consigna de corriente
    i_qs_ref = T_total / kt;
    i_ds_ref = 0;  % Estrategia base

    % 6. Controladores P de corriente
    v_qs_ctrl = Kp_iq * (i_qs_ref - i_qs);
    v_ds_ctrl = Kp_id * (i_ds_ref - i_ds);

    % 7. Desacoplamiento completo (compensación de feedbacks)
    v_qs = v_qs_ctrl + Rs_hat * i_qs + (lambda_m + Ld * i_ds) * omega_r;
    v_ds = v_ds_ctrl + Rs_hat * i_ds - Lq * i_qs * omega_r;
    v_0s = 0;

    % --- PLANTA NL (ECUACIONES DE ESTADO) ---
    T_em = 1.5 * Pp * (lambda_m * i_qs + (Ld - Lq) * i_ds * i_qs);
    T_grav = (grav * kl / r_red) * sin(theta_l);

    d_theta_m = omega_m;
    d_omega_m = (1/J_eq) * (T_em - b_eq * omega_m - T_grav - T_ld/r_red);
    d_i_qs = (1/Lq) * (v_qs - Rs_real * i_qs - (lambda_m + Ld * i_ds) * omega_r);
    d_i_ds = (1/Ld) * (v_ds - Rs_real * i_ds + Lq * i_qs * omega_r);
    d_i_0s = (1/Lls) * (v_0s - Rs_real * i_0s);

    P_loss = 1.5 * Rs_real * (i_qs^2 + i_ds^2 + 2 * i_0s^2);
    d_T_s = (1/Cts) * (P_loss - (T_s - Tamb) / Rts_amb);

    dxdt = [d_theta_m; d_omega_m; d_i_qs; d_i_ds; d_i_0s; d_T_s];
end

function dxdt = pmsm_lazo_abierto_ode(t, x, vqs_func, tld_func, ...
    Lq, Ld, Lls, Rs_ref, alpha_cu, Ts_ref, lambda_m, Pp, ...
    J_eq, b_eq, r_red, kl, grav, Cts, Rts_amb, Tamb)
% PMSM en lazo abierto (5.1.6) con leyes de control mínima + complementaria
% x = [theta_m; omega_m; i_qs; i_ds; i_0s; T_s]

    theta_m = x(1); omega_m = x(2);
    i_qs = x(3); i_ds = x(4); i_0s = x(5); T_s = x(6);

    omega_r = Pp * omega_m;
    theta_l = theta_m / r_red;
    Rs = Rs_ref * (1 + alpha_cu * (T_s - Ts_ref));

    vqs_consigna = vqs_func(t);
    T_ld = tld_func(t);

    % Leyes de control NL (mínima + complementaria)
    v_ds = -Lq * i_qs * omega_r;
    v_qs = vqs_consigna + omega_r * Ld * i_ds;
    v_0s = 0;

    T_em = 1.5 * Pp * (lambda_m * i_qs + (Ld - Lq) * i_ds * i_qs);
    T_grav = (grav * kl / r_red) * sin(theta_l);

    d_theta_m = omega_m;
    d_omega_m = (1/J_eq) * (T_em - b_eq * omega_m - T_grav - T_ld/r_red);
    d_i_qs = (1/Lq) * (v_qs - Rs * i_qs - (lambda_m + Ld * i_ds) * omega_r);
    d_i_ds = (1/Ld) * (v_ds - Rs * i_ds + Lq * i_qs * omega_r);
    d_i_0s = (1/Lls) * (v_0s - Rs * i_0s);

    P_loss = 1.5 * Rs * (i_qs^2 + i_ds^2 + 2 * i_0s^2);
    d_T_s = (1/Cts) * (P_loss - (T_s - Tamb) / Rts_amb);

    dxdt = [d_theta_m; d_omega_m; d_i_qs; d_i_ds; d_i_0s; d_T_s];
end
