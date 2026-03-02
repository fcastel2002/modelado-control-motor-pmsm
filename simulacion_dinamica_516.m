%% SIMULACIÓN DINÁMICA - CONSIGNA 5.1.6
%  Comparación: Modelo NL completo (con ley de control NL) vs LTI aumentado
%  Sub-ítems: a) Estados vs tiempo, b) Métricas transitorias,
%             c) Efecto de i_ds(0), d) Field forcing/weakening
%  Autores: Joaquín Calderón - Francisco Castel

clearvars; clc; close all;

%% 1. CARGAR PARÁMETROS FÍSICOS
if exist('parametros_sistema_completo.m', 'file')
    parametros_sistema_completo;
    disp('-> Parámetros físicos cargados exitosamente.');
else
    error('No se encuentra "parametros_sistema_completo.m".');
end

R_s0 = R_sREF;          % Resistencia nominal (punto de operación T_s = 20°C)
kt   = 1.5 * P_p * lambda_r;  % Constante de torque electromagnético

%% 2. DEFINICIÓN DE SEÑALES DE ENTRADA
% Amplitudes según consigna
V_QS_AMP = 19.596;   % [V] Amplitud de pulso en eje q
T_LD_AMP = 6.28;     % [N.m] Amplitud de pulso de torque de carga

% --- v_qs*(t): Pulso de consigna de tensión eje q ---
%  0V -> +19.596V (t=0.1s) -> 0V (t=0.7s) -> -19.596V (t=1.1s) -> 0V (t=1.7s)
vqs_func = @(t) V_QS_AMP * ((t >= 0.1 & t < 0.7) - (t >= 1.1 & t < 1.7));

% --- T_ld(t): Doble pulso de torque de carga ---
%  0 -> +6.28 (t=0.3s) -> -6.28 (t=0.5s) -> 0 (t=0.9s)
%    -> +6.28 (t=1.3s) -> -6.28 (t=1.5s) -> 0 (t=1.9s)
tld_func = @(t) T_LD_AMP * ( ...
    (t >= 0.3 & t < 0.5)  - (t >= 0.5 & t < 0.9) + ...
    (t >= 1.3 & t < 1.5)  - (t >= 1.5 & t < 1.9));

%% 3. MODELO LTI AUMENTADO (4 estados)
% Estado: x = [theta_m; omega_m; i_qs; i_ds]
% Entradas: u = [v_qs; T_ld]

A_aug = [0,    1,          0,          0;
         0, -b_eq/J_eq,  kt/J_eq,      0;
         0,    0,       -R_s0/L_q,      0;
         0,    0,          0,    -R_s0/L_d];

B_aug = [0,         0;
         0,    -1/(r*J_eq);
         1/L_q,      0;
         0,          0];

C_aug = eye(4);
D_aug = zeros(4, 2);

sys_lti = ss(A_aug, B_aug, C_aug, D_aug);

fprintf('\n-> Modelo LTI aumentado construido (4 estados).\n');
fprintf('   Polos: '); fprintf('%.2f%+.2fj  ', [real(eig(A_aug))'; imag(eig(A_aug))']); fprintf('\n');

%% 4. CONFIGURACIÓN DE SIMULACIÓN
t_final = 2.2;        % [s] Tiempo total
dt      = 5e-5;       % [s] Paso de tiempo (50 us)
t_sim   = (0:dt:t_final)';
N       = length(t_sim);

% Pre-computar señales de entrada como vectores
u_vqs = vqs_func(t_sim);
u_tld = tld_func(t_sim);
U_lti = [u_vqs, u_tld];  % Matriz de entrada para lsim

% Casos de condición inicial i_ds(0)
ids0_cases  = [0, 0.5, -0.5];
case_names  = {'i_{ds}(0) = 0 A', 'i_{ds}(0) = +0.5 A', 'i_{ds}(0) = -0.5 A'};
case_colors = {'b', 'r', [0 0.6 0]};  % Azul, Rojo, Verde

% Opciones del integrador (sistema rígido por modos eléctricos rápidos)
ode_opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 1e-3);

%% 5. SIMULACIÓN PRINCIPAL (3 casos)
fprintf('\n=== SIMULACIÓN DINÁMICA (3 casos de i_{ds}(0)) ===\n');

% Almacenamiento de resultados
res_lti = cell(3,1);   % Resultados LTI  [N x 4]
res_nl  = cell(3,1);   % Resultados NL   [N x 6]
t_nl    = cell(3,1);   % Vector de tiempo NL (adaptativo)

for k = 1:3
    ids0 = ids0_cases(k);
    fprintf('\n--- Caso %d: i_ds(0) = %+.1f A ---\n', k, ids0);

    % --- 5.1 Simulación LTI (lsim) ---
    x0_lti = [0; 0; 0; ids0];  % [theta_m; omega_m; i_qs; i_ds]
    [y_lti, ~] = lsim(sys_lti, U_lti, t_sim, x0_lti);
    res_lti{k} = y_lti;
    fprintf('   LTI: completado.\n');

    % --- 5.2 Simulación NL (ode15s) ---
    % Estado NL: [theta_m; omega_m; i_qs; i_ds; i_0s; T_s]
    x0_nl = [0; 0; 0; ids0; 0; T_amb];

    [t_ode, x_ode] = ode15s(@(t,x) pmsm_nl_ode(t, x, ...
        vqs_func, tld_func, 0, ...  % vds_ext = 0 para sub-items a,b,c
        L_q, L_d, L_ls, R_sREF, a_cu, T_sREF, lambda_r, P_p, ...
        J_eq, b_eq, r, k_l, g, C_ts, R_ts_amb, T_amb), ...
        [0, t_final], x0_nl, ode_opts);

    % Interpolar a grilla uniforme para comparación
    x_nl_interp = interp1(t_ode, x_ode, t_sim, 'pchip');
    res_nl{k}  = x_nl_interp;
    t_nl{k}    = t_sim;
    fprintf('   NL:  completado (%d pasos adaptativos).\n', length(t_ode));
end

%% =====================================================================
%  SECCIÓN (a): RESPUESTA DEL ESTADO INTERNO
%  =====================================================================
fprintf('\n=== GENERANDO GRÁFICOS - SECCIÓN (a) ===\n');

state_labels_lti = {'\theta_m [rad]', '\omega_m [rad/s]', 'i_{qs} [A]', 'i_{ds} [A]'};
state_labels_nl  = {'\theta_m [rad]', '\omega_m [rad/s]', 'i_{qs} [A]', 'i_{ds} [A]', 'i_{0s} [A]', 'T_s [°C]'};

% --- Figura 1: Señales de Entrada ---
fig1 = figure('Name', 'Señales de Entrada', 'Color', 'w', 'Position', [50 400 900 350]);

subplot(2,1,1);
plot(t_sim, u_vqs, 'b', 'LineWidth', 1.2);
ylabel('v_{qs}^* [V]'); grid on;
title('Señales de Entrada según Consigna', 'FontSize', 13);
xlim([0 t_final]);

subplot(2,1,2);
plot(t_sim, u_tld, 'r', 'LineWidth', 1.2);
ylabel('T_{ld} [N.m]'); xlabel('Tiempo [s]'); grid on;
xlim([0 t_final]);

% --- Figura 2: Estados mecánicos (θ_m, ω_m) - NL vs LTI ---
fig2 = figure('Name', 'Estados Mecánicos - NL vs LTI', 'Color', 'w', 'Position', [50 50 1200 700]);

for k = 1:3
    % Posición angular
    subplot(3,2, 2*k-1); hold on; grid on;
    plot(t_sim, res_lti{k}(:,1), '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
    plot(t_sim, res_nl{k}(:,1),  '--', 'Color', case_colors{k}, 'LineWidth', 1.2);
    ylabel('\theta_m [rad]');
    if k == 1; title('Posición Angular', 'FontSize', 12); end
    if k == 3; xlabel('Tiempo [s]'); end
    legend('LTI', 'NL', 'Location', 'best');
    text(0.02, 0.9, case_names{k}, 'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'w');
    xlim([0 t_final]);

    % Velocidad angular
    subplot(3,2, 2*k); hold on; grid on;
    plot(t_sim, res_lti{k}(:,2), '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
    plot(t_sim, res_nl{k}(:,2),  '--', 'Color', case_colors{k}, 'LineWidth', 1.2);
    ylabel('\omega_m [rad/s]');
    if k == 1; title('Velocidad Angular', 'FontSize', 12); end
    if k == 3; xlabel('Tiempo [s]'); end
    legend('LTI', 'NL', 'Location', 'best');
    xlim([0 t_final]);
end

sgtitle('Estados Mecánicos: Modelo NL vs LTI Aumentado', 'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 3: Corrientes (i_qs, i_ds) - NL vs LTI ---
fig3 = figure('Name', 'Corrientes qd - NL vs LTI', 'Color', 'w', 'Position', [100 50 1200 700]);

for k = 1:3
    % i_qs
    subplot(3,2, 2*k-1); hold on; grid on;
    plot(t_sim, res_lti{k}(:,3), '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
    plot(t_sim, res_nl{k}(:,3),  '--', 'Color', case_colors{k}, 'LineWidth', 1.2);
    ylabel('i_{qs} [A]');
    if k == 1; title('Corriente eje q', 'FontSize', 12); end
    if k == 3; xlabel('Tiempo [s]'); end
    legend('LTI', 'NL', 'Location', 'best');
    text(0.02, 0.9, case_names{k}, 'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'w');
    xlim([0 t_final]);

    % i_ds
    subplot(3,2, 2*k); hold on; grid on;
    plot(t_sim, res_lti{k}(:,4), '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
    plot(t_sim, res_nl{k}(:,4),  '--', 'Color', case_colors{k}, 'LineWidth', 1.2);
    ylabel('i_{ds} [A]');
    if k == 1; title('Corriente eje d', 'FontSize', 12); end
    if k == 3; xlabel('Tiempo [s]'); end
    legend('LTI', 'NL', 'Location', 'best');
    xlim([0 t_final]);
end

sgtitle('Corrientes en Coordenadas qd0: Modelo NL vs LTI Aumentado', 'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 4: Estados exclusivos del NL (i_0s, T_s) ---
fig4 = figure('Name', 'Estados NL exclusivos', 'Color', 'w', 'Position', [150 100 900 500]);

subplot(2,1,1); hold on; grid on;
for k = 1:3
    plot(t_sim, res_nl{k}(:,5), '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
end
ylabel('i_{0s} [A]'); title('Corriente de secuencia cero (solo NL)', 'FontSize', 12);
legend(case_names, 'Location', 'best'); xlim([0 t_final]);

subplot(2,1,2); hold on; grid on;
for k = 1:3
    plot(t_sim, res_nl{k}(:,6), '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
end
ylabel('T_s [°C]'); xlabel('Tiempo [s]');
title('Temperatura del Estator (solo NL)', 'FontSize', 12);
legend(case_names, 'Location', 'best'); xlim([0 t_final]);

sgtitle('Variables Exclusivas del Modelo No Lineal', 'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 5: v_ds forzada (ley de control NL) ---
fig5 = figure('Name', 'Tensión v_ds forzada', 'Color', 'w', 'Position', [200 100 900 350]);

for k = 1:3
    omega_r_k = P_p * res_nl{k}(:,2);    % omega_r = P_p * omega_m
    i_qs_k    = res_nl{k}(:,3);
    vds_forced = -L_q .* i_qs_k .* omega_r_k;  % Ley de control mínima

    subplot(1,3,k); hold on; grid on;
    plot(t_sim, vds_forced, '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
    ylabel('v_{ds} [V]'); xlabel('Tiempo [s]');
    title(case_names{k}, 'FontSize', 11);
    xlim([0 t_final]);
end

sgtitle('Tensión v_{ds}(t) Forzada por Ley de Control Mínima: v_{ds} = -L_q \cdot i_{qs} \cdot \omega_r', ...
    'FontSize', 12, 'FontWeight', 'bold');

% --- Figura 6: Corrientes en coordenadas abc (Transformada Inversa de Park) ---
fig6 = figure('Name', 'Corrientes abc (NL)', 'Color', 'w', 'Position', [250 50 1200 700]);

for k = 1:3
    theta_e = P_p * res_nl{k}(:,1);  % Ángulo eléctrico
    i_q = res_nl{k}(:,3);
    i_d = res_nl{k}(:,4);
    i_0 = res_nl{k}(:,5);

    % Transformada inversa de Park (qd0 -> abc)
    i_a = i_q .* cos(theta_e)            + i_d .* sin(theta_e)            + i_0;
    i_b = i_q .* cos(theta_e - 2*pi/3)   + i_d .* sin(theta_e - 2*pi/3)  + i_0;
    i_c = i_q .* cos(theta_e + 2*pi/3)   + i_d .* sin(theta_e + 2*pi/3)  + i_0;

    subplot(3,1,k); hold on; grid on;
    plot(t_sim, i_a, 'r', 'LineWidth', 0.8);
    plot(t_sim, i_b, 'g', 'LineWidth', 0.8);
    plot(t_sim, i_c, 'b', 'LineWidth', 0.8);
    ylabel('Corriente [A]');
    if k == 3; xlabel('Tiempo [s]'); end
    title(['Corrientes abc (NL) - ', case_names{k}], 'FontSize', 11);
    legend('i_a', 'i_b', 'i_c', 'Location', 'northeast');
    xlim([0 t_final]);
end

sgtitle('Corrientes de Fase (Coordenadas abc) - Modelo NL', 'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 7: Tensiones en coordenadas abc ---
fig7 = figure('Name', 'Tensiones abc (NL)', 'Color', 'w', 'Position', [300 50 1200 700]);

for k = 1:3
    theta_e = P_p * res_nl{k}(:,1);
    omega_r = P_p * res_nl{k}(:,2);
    i_q_k = res_nl{k}(:,3);
    i_d_k = res_nl{k}(:,4);

    % Tensiones en qd0 (con leyes de control)
    v_q_applied = u_vqs + omega_r .* L_d .* i_d_k;   % Complementaria
    v_d_applied = -L_q .* i_q_k .* omega_r;            % Mínima
    v_0_applied = zeros(N, 1);

    % Inversa de Park
    v_a = v_q_applied .* cos(theta_e)            + v_d_applied .* sin(theta_e)            + v_0_applied;
    v_b = v_q_applied .* cos(theta_e - 2*pi/3)   + v_d_applied .* sin(theta_e - 2*pi/3)  + v_0_applied;
    v_c = v_q_applied .* cos(theta_e + 2*pi/3)   + v_d_applied .* sin(theta_e + 2*pi/3)  + v_0_applied;

    subplot(3,1,k); hold on; grid on;
    plot(t_sim, v_a, 'r', 'LineWidth', 0.8);
    plot(t_sim, v_b, 'g', 'LineWidth', 0.8);
    plot(t_sim, v_c, 'b', 'LineWidth', 0.8);
    ylabel('Tensión [V]');
    if k == 3; xlabel('Tiempo [s]'); end
    title(['Tensiones abc (NL) - ', case_names{k}], 'FontSize', 11);
    legend('v_a', 'v_b', 'v_c', 'Location', 'northeast');
    xlim([0 t_final]);
end

sgtitle('Tensiones de Fase (Coordenadas abc) - Modelo NL', 'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 8: Curvas paramétricas ---
fig8 = figure('Name', 'Curvas Paramétricas', 'Color', 'w', 'Position', [350 50 1200 500]);

% Torque vs Velocidad (plano de operación)
subplot(1,2,1); hold on; grid on;
for k = 1:3
    % Torque electromagnético del NL
    T_em = kt .* res_nl{k}(:,3);  % T_em = kt * i_qs (con i_ds ≈ 0)
    omega_m = res_nl{k}(:,2);
    plot(omega_m, T_em, '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
end
xlabel('\omega_m [rad/s]'); ylabel('T_{em} [N.m]');
title('Torque vs Velocidad (Cuadrantes de Operación)', 'FontSize', 12);
legend(case_names, 'Location', 'best');
xline(0, 'k--'); yline(0, 'k--');

% i_ds vs i_qs (plano de corrientes)
subplot(1,2,2); hold on; grid on;
for k = 1:3
    plot(res_nl{k}(:,3), res_nl{k}(:,4), '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
    % Marcar inicio
    plot(res_nl{k}(1,3), res_nl{k}(1,4), 'o', 'Color', case_colors{k}, 'MarkerSize', 8, 'MarkerFaceColor', case_colors{k});
end
xlabel('i_{qs} [A]'); ylabel('i_{ds} [A]');
title('Trayectoria en Plano de Corrientes (i_{ds} vs i_{qs})', 'FontSize', 12);
legend(case_names, 'Location', 'best');
xline(0, 'k--'); yline(0, 'k--');

sgtitle('Curvas Paramétricas - Modelo NL', 'FontSize', 14, 'FontWeight', 'bold');

% --- Figura 9: Ángulo de torque del rotor ---
fig9 = figure('Name', 'Ángulo de Torque', 'Color', 'w', 'Position', [400 100 900 400]);

for k = 1:3
    delta_torque = atan2(res_nl{k}(:,4), res_nl{k}(:,3));  % atan2(i_ds, i_qs)

    subplot(1,3,k); hold on; grid on;
    plot(t_sim, delta_torque * 180/pi, '-', 'Color', case_colors{k}, 'LineWidth', 1.2);
    ylabel('\delta [grados]'); xlabel('Tiempo [s]');
    title(case_names{k}, 'FontSize', 11);
    xlim([0 t_final]);
end

sgtitle('Evolución del Ángulo de Torque del Rotor (\delta = atan2(i_{ds}, i_{qs}))', ...
    'FontSize', 12, 'FontWeight', 'bold');

%% =====================================================================
%  SECCIÓN (b): MÉTRICAS DE TRANSITORIO
%  =====================================================================
fprintf('\n=== SECCIÓN (b): MÉTRICAS DE TRANSITORIO ===\n');

% Analizamos el primer transitorio (pulso positivo de v_qs: t=0.1s a t=0.7s)
% para omega_m (velocidad) e i_qs (corriente)

fprintf('\n%-20s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
    'Caso', 'w_ss [rad/s]', 'iq_ss [A]', 'tr [ms]', 'ts [ms]', 'Mp [%%]');
fprintf('%s\n', repmat('-', 1, 85));

for k = 1:3
    % Extraer segmento del primer pulso positivo
    idx_start = find(t_sim >= 0.1, 1);
    idx_end   = find(t_sim >= 0.7, 1);
    t_seg  = t_sim(idx_start:idx_end) - t_sim(idx_start);

    % Velocidad (LTI)
    w_seg = res_lti{k}(idx_start:idx_end, 2);
    w_ss  = w_seg(end);  % Valor estacionario (aproximado)

    % Corriente iq (LTI)
    iq_seg = res_lti{k}(idx_start:idx_end, 3);
    iq_ss  = iq_seg(end);

    % Tiempo de crecimiento (10% a 90% del valor final de velocidad)
    if abs(w_ss) > 1e-10
        t10 = t_seg(find(abs(w_seg) >= 0.1*abs(w_ss), 1));
        t90 = t_seg(find(abs(w_seg) >= 0.9*abs(w_ss), 1));
        if ~isempty(t10) && ~isempty(t90)
            tr = (t90 - t10) * 1000;  % [ms]
        else
            tr = NaN;
        end
    else
        tr = NaN;
    end

    % Tiempo de establecimiento (±1% del valor final)
    if abs(w_ss) > 1e-10
        settled = abs(w_seg - w_ss) <= 0.01 * abs(w_ss);
        idx_settled = find(settled, 1, 'first');
        % Verificar que se mantiene
        if ~isempty(idx_settled) && all(settled(idx_settled:end))
            ts_val = t_seg(idx_settled) * 1000;  % [ms]
        else
            ts_val = NaN;
        end
    else
        ts_val = NaN;
    end

    % Sobrepico
    if abs(w_ss) > 1e-10
        Mp = (max(abs(w_seg)) - abs(w_ss)) / abs(w_ss) * 100;
    else
        Mp = NaN;
    end

    fprintf('%-20s | %12.4f | %12.4f | %12.2f | %12.2f | %12.2f\n', ...
        case_names{k}, w_ss, iq_ss, tr, ts_val, Mp);
end

fprintf('\n-> Nota: Las métricas son para el 1er pulso positivo (t=0.1s a 0.7s) del modelo LTI.\n');
fprintf('   La velocidad crece como rampa (integrador), el transitorio rápido es de i_qs.\n');

% Influencia relativa de v_qs vs T_ld
fprintf('\n--- Influencia relativa de las acciones externas ---\n');
fprintf('v_qs* actúa directamente sobre i_qs (modo rápido eléctrico, tau_e = L_q/R_s = %.2f ms)\n', L_q/R_s0*1000);
fprintf('T_ld  actúa sobre omega_m (modo lento mecánico, tau_m = J_eq/b_eq = %.2f s)\n', J_eq/b_eq);
fprintf('=> v_qs domina la dinámica de corriente; T_ld perturba la dinámica mecánica.\n');

%% =====================================================================
%  SECCIÓN (c): COMPARACIÓN i_ds(t) PARA LOS 3 CASOS
%  =====================================================================
fprintf('\n=== SECCIÓN (c): EFECTO DE i_ds(0) ===\n');

fig10 = figure('Name', 'Comparación i_ds(t)', 'Color', 'w', 'Position', [450 100 1000 600]);

% Comparación en LTI
subplot(2,1,1); hold on; grid on;
for k = 1:3
    plot(t_sim, res_lti{k}(:,4), '-', 'Color', case_colors{k}, 'LineWidth', 1.5);
end
ylabel('i_{ds} [A]'); title('Modelo LTI Aumentado', 'FontSize', 12);
legend(case_names, 'Location', 'northeast');
xlim([0 t_final]);

% Comparación en NL
subplot(2,1,2); hold on; grid on;
for k = 1:3
    plot(t_sim, res_nl{k}(:,4), '-', 'Color', case_colors{k}, 'LineWidth', 1.5);
end
ylabel('i_{ds} [A]'); xlabel('Tiempo [s]');
title('Modelo No Lineal', 'FontSize', 12);
legend(case_names, 'Location', 'northeast');
xlim([0 t_final]);

sgtitle('Comparación de i_{ds}(t) para Distintas Condiciones Iniciales', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Constante de tiempo de decaimiento
tau_ids = L_d / R_s0;
fprintf('Constante de tiempo de i_ds: tau = L_d/R_s0 = %.4f s (%.2f ms)\n', tau_ids, tau_ids*1000);
fprintf('En el LTI, i_ds decae exponencialmente con tau = %.2f ms (desacoplado).\n', tau_ids*1000);
fprintf('En el NL, i_ds también decae pero con acoplamiento cruzado leve vía omega_r.\n');
fprintf('=> El efecto de i_ds(0) != 0 se extingue rápidamente en ambos modelos.\n');

%% =====================================================================
%  SECCIÓN (d): FIELD FORCING / WEAKENING
%  =====================================================================
fprintf('\n=== SECCIÓN (d): FIELD FORCING / WEAKENING ===\n');

% Consigna: v_ds*(t) = 0 -> ±1.9596 V en t = 0.5 s
V_DS_FF = V_QS_AMP / 10;  % 1.9596 V

% Definir señal v_ds* para forcing y weakening
vds_forcing  = @(t) V_DS_FF * (t >= 0.5);    % +1.9596V (field forcing)
vds_weakening = @(t) -V_DS_FF * (t >= 0.5);  % -1.9596V (field weakening)

vds_cases = {vds_forcing, vds_weakening};
vds_labels = {'Field Forcing (+1.96V)', 'Field Weakening (-1.96V)'};

% Solo simulamos con i_ds(0) = 0 para comparar con caso base
res_ff = cell(2,1);  % {forcing, weakening}

for j = 1:2
    x0_nl = [0; 0; 0; 0; 0; T_amb];

    [t_ode, x_ode] = ode15s(@(t,x) pmsm_nl_ode(t, x, ...
        vqs_func, tld_func, vds_cases{j}, ...
        L_q, L_d, L_ls, R_sREF, a_cu, T_sREF, lambda_r, P_p, ...
        J_eq, b_eq, r, k_l, g, C_ts, R_ts_amb, T_amb), ...
        [0, t_final], x0_nl, ode_opts);

    res_ff{j} = interp1(t_ode, x_ode, t_sim, 'pchip');
    fprintf('   %s: completado.\n', vds_labels{j});
end

% Simulación LTI con v_ds como segunda entrada
% Redefinir LTI con 3 entradas: [v_qs, T_ld, v_ds]
B_aug_3 = [B_aug, [0; 0; 0; 1/L_d]];  % Agregar columna para v_ds
D_aug_3 = zeros(4, 3);
sys_lti_3 = ss(A_aug, B_aug_3, C_aug, D_aug_3);

res_ff_lti = cell(2,1);
for j = 1:2
    if j == 1
        u_vds = V_DS_FF * (t_sim >= 0.5);
    else
        u_vds = -V_DS_FF * (t_sim >= 0.5);
    end
    U_lti_3 = [u_vqs, u_tld, u_vds];
    x0_lti = [0; 0; 0; 0];
    [y_lti_3, ~] = lsim(sys_lti_3, U_lti_3, t_sim, x0_lti);
    res_ff_lti{j} = y_lti_3;
end

% --- Figura 11: Field Forcing/Weakening - Estados ---
fig11 = figure('Name', 'Field Forcing/Weakening', 'Color', 'w', 'Position', [500 50 1300 800]);

vars_plot = {'\theta_m [rad]', '\omega_m [rad/s]', 'i_{qs} [A]', 'i_{ds} [A]'};

for v = 1:4
    subplot(2,2,v); hold on; grid on;

    % Caso base (sin v_ds*)
    plot(t_sim, res_nl{1}(:,v), 'k-', 'LineWidth', 1, 'DisplayName', 'Base (v_{ds}^*=0)');
    plot(t_sim, res_lti{1}(:,v), 'k--', 'LineWidth', 1, 'DisplayName', 'Base LTI');

    % Field forcing (NL y LTI)
    plot(t_sim, res_ff{1}(:,v), 'r-', 'LineWidth', 1.2, 'DisplayName', 'FF NL (+1.96V)');
    plot(t_sim, res_ff_lti{1}(:,v), 'r--', 'LineWidth', 1.2, 'DisplayName', 'FF LTI');

    % Field weakening (NL y LTI)
    plot(t_sim, res_ff{2}(:,v), 'b-', 'LineWidth', 1.2, 'DisplayName', 'FW NL (-1.96V)');
    plot(t_sim, res_ff_lti{2}(:,v), 'b--', 'LineWidth', 1.2, 'DisplayName', 'FW LTI');

    ylabel(vars_plot{v});
    if v >= 3; xlabel('Tiempo [s]'); end
    if v == 1; legend('Location', 'best', 'FontSize', 7); end
    xlim([0 t_final]);
end

sgtitle('Efecto del Field Forcing/Weakening (v_{ds}^* \neq 0) sobre el Sistema', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Análisis del efecto
fprintf('\n--- Análisis del Field Forcing/Weakening ---\n');
fprintf('Al aplicar v_ds* = +%.3f V (field forcing) en t=0.5s:\n', V_DS_FF);
fprintf('  -> i_ds crece (inyección de flujo en eje d)\n');
fprintf('  -> Modifica el torque electromagnético vía (L_d - L_q)*i_d*i_q\n');
fprintf('  -> En el LTI aumentado, solo afecta i_ds (desacoplado del resto)\n');
fprintf('  -> En el NL, el acoplamiento cruzado introduce efectos en i_qs y omega_m\n');
fprintf('Al aplicar v_ds* = -%.3f V (field weakening):\n', V_DS_FF);
fprintf('  -> Efecto simétrico opuesto\n');

%% =====================================================================
%  RESUMEN FINAL
%  =====================================================================
fprintf('\n========================================\n');
fprintf('  FIGURAS GENERADAS (guardar manualmente)\n');
fprintf('========================================\n');
fprintf('  Fig 1:  senales_entrada.png\n');
fprintf('  Fig 2:  estados_mecanicos_NLvsLTI.png\n');
fprintf('  Fig 3:  corrientes_qd_NLvsLTI.png\n');
fprintf('  Fig 4:  estados_NL_exclusivos.png\n');
fprintf('  Fig 5:  vds_forzada.png\n');
fprintf('  Fig 6:  corrientes_abc_NL.png\n');
fprintf('  Fig 7:  tensiones_abc_NL.png\n');
fprintf('  Fig 8:  curvas_parametricas.png\n');
fprintf('  Fig 9:  angulo_torque.png\n');
fprintf('  Fig 10: comparacion_ids.png\n');
fprintf('  Fig 11: field_forcing_weakening.png\n');
fprintf('========================================\n');
fprintf('¡Simulación completada exitosamente!\n');


%% =====================================================================
%  FUNCIONES LOCALES
%  =====================================================================

function dxdt = pmsm_nl_ode(t, x, vqs_func, tld_func, vds_ext_func, ...
    Lq, Ld, Lls, Rs_ref, alpha_cu, Ts_ref, lambda_m, Pp, ...
    J_eq, b_eq, r_red, kl, grav, Cts, Rts_amb, Tamb)
% PMSM_NL_ODE  Ecuaciones de estado del PMSM no lineal completo
%   x = [theta_m; omega_m; i_qs; i_ds; i_0s; T_s]
%   Incluye leyes de control NL (mínima + complementaria)

    % Desempaquetar estados
    theta_m = x(1);
    omega_m = x(2);
    i_qs    = x(3);
    i_ds    = x(4);
    i_0s    = x(5);
    T_s     = x(6);

    % Variables derivadas
    omega_r = Pp * omega_m;          % Velocidad eléctrica
    theta_l = theta_m / r_red;       % Posición de carga

    % Resistencia variable con temperatura
    Rs = Rs_ref * (1 + alpha_cu * (T_s - Ts_ref));

    % Entradas externas
    vqs_consigna = vqs_func(t);       % Señal de consigna v_qs*(t)
    T_ld         = tld_func(t);       % Torque de carga

    % Tensión v_ds externa (para field forcing/weakening, item d)
    if isa(vds_ext_func, 'function_handle')
        vds_ext = vds_ext_func(t);
    else
        vds_ext = 0;
    end

    % --- LEYES DE CONTROL NO LINEALES ---
    % Ley mínima (desacoplamiento de i_ds):
    v_ds = -Lq * i_qs * omega_r + vds_ext;

    % Ley complementaria (cancelación de acoplamiento cruzado en eje q):
    v_qs = vqs_consigna + omega_r * Ld * i_ds;

    % v_0s = 0 (sin secuencia cero aplicada)
    v_0s = 0;

    % --- ECUACIONES DE ESTADO ---
    % 1. Posición mecánica
    d_theta_m = omega_m;

    % 2. Velocidad mecánica
    T_em = 1.5 * Pp * (lambda_m * i_qs + (Ld - Lq) * i_ds * i_qs);  % Torque EM
    T_grav = (grav * kl / r_red) * sin(theta_l);                      % Torque gravitacional
    d_omega_m = (1/J_eq) * (T_em - b_eq * omega_m - T_grav - T_ld/r_red);

    % 3. Corriente eje q
    d_i_qs = (1/Lq) * (v_qs - Rs * i_qs - (lambda_m + Ld * i_ds) * omega_r);

    % 4. Corriente eje d
    d_i_ds = (1/Ld) * (v_ds - Rs * i_ds + Lq * i_qs * omega_r);

    % 5. Corriente secuencia cero
    d_i_0s = (1/Lls) * (v_0s - Rs * i_0s);

    % 6. Temperatura estator
    P_loss = 1.5 * Rs * (i_qs^2 + i_ds^2 + 2 * i_0s^2);
    d_T_s = (1/Cts) * (P_loss - (T_s - Tamb) / Rts_amb);

    dxdt = [d_theta_m; d_omega_m; d_i_qs; d_i_ds; d_i_0s; d_T_s];
end
