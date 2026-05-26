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

R_s0 = R_sREF * (1 + 3.9e-3*(29.5-20));   % R_s en operacion T_s=29.5C (media de la sim 5.1.6), NO a 20C. = 1.058 Ohm
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

% --- Figura 8: Torque vs Velocidad (Cuadrantes de Operación) - GRÁFICA INDEPENDIENTE ---
fig8 = figure('Name', 'Torque vs Velocidad - Cuadrantes de Operación', ...
    'Color', 'w', 'Position', [100 50 1100 780]);

hold on; grid on; grid minor; box on;
set(gca, 'FontSize', 11, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);

% === Datos del caso principal (k=1, i_ds(0) = 0 A) ===
T_em_main = kt .* res_nl{1}(:,3);
omega_m_main = res_nl{1}(:,2);

% === Casos secundarios k=2,3 como referencia en fondo ===
colors_bg = {[0.75 0.75 0.95], [0.75 0.95 0.75]};  % Azul/verde claros
for k = 2:3
    T_em_bg = kt .* res_nl{k}(:,3);
    omega_bg = res_nl{k}(:,2);
    plot(omega_bg, T_em_bg, '-', 'Color', colors_bg{k-1}, 'LineWidth', 1.0, ...
        'DisplayName', ['Ref: ', case_names{k}]);
end

% === Definición de transitorios (10 segmentos) ===
t_ev = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, t_final];
n_tr = length(t_ev) - 1;

trans_names = { ...
    'Trans. 1:  v_{qs}^* \rightarrow +19.6 V', ...
    'Trans. 2:  T_{ld} \rightarrow +6.28 Nm', ...
    'Trans. 3:  T_{ld} \rightarrow -6.28 Nm', ...
    'Trans. 4:  v_{qs}^* \rightarrow 0 V', ...
    'Trans. 5:  T_{ld} \rightarrow 0 Nm', ...
    'Trans. 6:  v_{qs}^* \rightarrow -19.6 V', ...
    'Trans. 7:  T_{ld} \rightarrow +6.28 Nm', ...
    'Trans. 8:  T_{ld} \rightarrow -6.28 Nm', ...
    'Trans. 9:  v_{qs}^* \rightarrow 0 V', ...
    'Trans. 10: T_{ld} \rightarrow 0 Nm'};

% Paleta de colores distinguibles (cálidos T1-T5, fríos T6-T10)
cmap_tr = [ ...
    0.84 0.15 0.16;   % T1  rojo
    0.99 0.55 0.24;   % T2  naranja
    0.17 0.63 0.17;   % T3  verde
    0.12 0.47 0.71;   % T4  azul
    0.58 0.40 0.74;   % T5  violeta
    0.09 0.75 0.81;   % T6  cian
    0.89 0.47 0.76;   % T7  rosa
    0.55 0.34 0.29;   % T8  marrón
    0.80 0.73 0.13;   % T9  oliva
    0.35 0.35 0.35];  % T10 gris oscuro

% Estilos de línea alternados para mayor distinción
line_styles = {'-', '-', '-', '-', '-', '--', '--', '--', '--', '--'};

% === Plotear cada transitorio ===
eq_pts = zeros(n_tr + 1, 2);  % [omega_m, T_em] puntos de equilibrio
eq_pts(1,:) = [0, 0];         % P_0: reposo inicial

for s = 1:n_tr
    if s < n_tr
        idx_s = find(t_sim >= t_ev(s) & t_sim < t_ev(s+1));
    else
        idx_s = find(t_sim >= t_ev(s) & t_sim <= t_ev(s+1));
    end

    plot(omega_m_main(idx_s), T_em_main(idx_s), line_styles{s}, ...
        'Color', cmap_tr(s,:), 'LineWidth', 2.2, 'DisplayName', trans_names{s});

    % Guardar punto de equilibrio (final del segmento)
    eq_pts(s+1,:) = [omega_m_main(idx_s(end)), T_em_main(idx_s(end))];

    % Flecha de dirección en punto medio (quiver)
    mid = round(length(idx_s) * 0.45);
    step = max(1, round(length(idx_s) * 0.02));
    if mid > step && mid + step <= length(idx_s)
        dx = omega_m_main(idx_s(mid+step)) - omega_m_main(idx_s(mid-step));
        dy = T_em_main(idx_s(mid+step)) - T_em_main(idx_s(mid-step));
        quiver(omega_m_main(idx_s(mid)), T_em_main(idx_s(mid)), dx, dy, 1.5, ...
            'Color', cmap_tr(s,:), 'LineWidth', 1.8, 'MaxHeadSize', 3, ...
            'HandleVisibility', 'off');
    end
end

% === Puntos de equilibrio (diamantes) ===
plot(eq_pts(1,1), eq_pts(1,2), 'kp', 'MarkerSize', 16, ...
    'MarkerFaceColor', [0.2 0.8 0.2], 'LineWidth', 1.5, ...
    'DisplayName', 'P_0: Inicio (reposo)');

for s = 1:n_tr
    plot(eq_pts(s+1,1), eq_pts(s+1,2), 'd', 'Color', cmap_tr(s,:), ...
        'MarkerSize', 9, 'MarkerFaceColor', cmap_tr(s,:), 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end

% Agregar un solo marcador representativo para la leyenda
plot(NaN, NaN, 'kd', 'MarkerSize', 9, 'MarkerFaceColor', [1 0.84 0], ...
    'LineWidth', 1.0, 'DisplayName', 'Puntos de equilibrio (P_1..P_{10})');

% Etiquetas de puntos de equilibrio
for s = 0:n_tr
    % Desplazamientos para evitar solapamiento
    if eq_pts(s+1,1) >= 0; ox = 1.5; else; ox = -1.5; end
    if eq_pts(s+1,2) >= 0; oy = 0.003; else; oy = -0.003; end
    text(eq_pts(s+1,1) + ox, eq_pts(s+1,2) + oy, ...
        sprintf('P_{%d}', s), 'FontSize', 8.5, 'FontWeight', 'bold', ...
        'Color', [0.1 0.1 0.1], 'BackgroundColor', [1 1 1 0.7], ...
        'Margin', 1, 'HorizontalAlignment', 'center');
end

% === Ejes de cuadrantes ===
xline(0, 'k-', 'LineWidth', 1.8, 'HandleVisibility', 'off');
yline(0, 'k-', 'LineWidth', 1.8, 'HandleVisibility', 'off');

% === Etiquetas de cuadrantes ===
ax_lims = axis;
margin_x = 0.70;  margin_y = 0.82;
text( abs(ax_lims(2))*margin_x,  abs(ax_lims(4))*margin_y, ...
    {'{\bf Cuadrante I}', 'Motor (+\omega, +T)'}, ...
    'FontSize', 10, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'center');
text(-abs(ax_lims(2))*margin_x,  abs(ax_lims(4))*margin_y, ...
    {'{\bf Cuadrante II}', 'Freno Reg. (-\omega, +T)'}, ...
    'FontSize', 10, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'center');
text(-abs(ax_lims(2))*margin_x, -abs(ax_lims(4))*margin_y, ...
    {'{\bf Cuadrante III}', 'Motor (-\omega, -T)'}, ...
    'FontSize', 10, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'center');
text( abs(ax_lims(2))*margin_x, -abs(ax_lims(4))*margin_y, ...
    {'{\bf Cuadrante IV}', 'Freno Reg. (+\omega, -T)'}, ...
    'FontSize', 10, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'center');

% === Formato final ===
xlabel('\omega_m [rad/s]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('T_{em} [N \cdot m]', 'FontSize', 13, 'FontWeight', 'bold');
title({'Torque Electromagnético vs Velocidad Angular - Cuadrantes de Operación'; ...
       'Modelo NL con Ley de Control (caso i_{ds}(0) = 0 A)'}, ...
    'FontSize', 13, 'FontWeight', 'bold');

lgd = legend('Location', 'eastoutside', 'FontSize', 8.5, 'Box', 'on');
title(lgd, 'Transitorios y Referencias', 'FontSize', 9);
exportgraphics(fig8, 'imagenes/sim_torque_vs_velocidad.png', 'Resolution', 150);

% --- Figura 8b: Trayectoria en Plano de Corrientes (3 casos en paneles separados) ---
fig8b = figure('Name', 'Plano de Corrientes i_ds vs i_qs', 'Color', 'w', ...
    'Position', [200 40 700 850]);

for k = 1:3
    subplot(3,1,k); hold on; grid on; box on;
    plot(res_nl{k}(:,3), res_nl{k}(:,4), '-', 'Color', case_colors{k}, 'LineWidth', 1.4);
    plot(res_nl{k}(1,3), res_nl{k}(1,4), 'o', 'Color', case_colors{k}, ...
        'MarkerSize', 9, 'MarkerFaceColor', case_colors{k}, 'HandleVisibility', 'off');
    xline(0, 'k--', 'HandleVisibility', 'off'); yline(0, 'k--', 'HandleVisibility', 'off');
    ylabel('i_{ds} [A]', 'FontSize', 11);
    title(case_names{k}, 'FontSize', 12, 'FontWeight', 'bold');
    ylim([-0.6 0.6]);   % i_ds es chico: fijamos escala comun para apreciar el decaimiento
end
xlabel('i_{qs} [A]', 'FontSize', 12);   % solo en el panel inferior
sgtitle('Trayectoria en Plano de Corrientes (i_{ds} vs i_{qs}) - Modelo NL', ...
    'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig8b, 'imagenes/sim_ids_vs_iqs.png', 'Resolution', 150);

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
%  SECCIÓN (b): MÉTRICAS DE TRANSITORIO (10 transitorios, caso i_ds(0)=0)
%  =====================================================================
fprintf('\n=== SECCIÓN (b): MÉTRICAS DE TRANSITORIO ===\n');

% Instantes de conmutación y descripción de cada evento
t_events = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, t_final];
event_names = { ...
    'v_{qs}^* -> +19.6V', ...
    'T_{ld} -> +6.28Nm',  ...
    'T_{ld} -> -6.28Nm',  ...
    'v_{qs}^* -> 0V',     ...
    'T_{ld} -> 0Nm',      ...
    'v_{qs}^* -> -19.6V', ...
    'T_{ld} -> +6.28Nm',  ...
    'T_{ld} -> -6.28Nm',  ...
    'v_{qs}^* -> 0V',     ...
    'T_{ld} -> 0Nm'};

% Usamos caso k=1 (i_ds(0) = 0) del modelo LTI
y_lti_nom = res_lti{1};  % [theta_m, omega_m, i_qs, i_ds]

fprintf('\n--- Métricas por transitorio (modelo LTI, i_ds(0) = 0 A) ---\n');
fprintf('%-6s | %-25s | %10s | %10s | %10s | %10s | %8s\n', ...
    'Trans', 'Evento', 'w_end', 'iq_ss', 'tr_iq', 'ts_iq', 'Mp_iq');
fprintf('%-6s | %-25s | %10s | %10s | %10s | %10s | %8s\n', ...
    '', '', '[rad/s]', '[A]', '[ms]', '[ms]', '[%]');
fprintf('%s\n', repmat('-', 1, 95));

% Almacenar para LaTeX
metricas = NaN(10, 5);  % [w_end, iq_ss, tr_iq, ts_iq, Mp_iq]
latex_rows = strings(10, 1);

for s = 1:10
    idx_ini = find(t_sim >= t_events(s), 1);
    if s < 10
        idx_fin = find(t_sim < t_events(s+1), 1, 'last');
    else
        idx_fin = length(t_sim);
    end
    if isempty(idx_fin) || idx_fin < idx_ini
        idx_fin = length(t_sim);
    end

    t_seg = t_sim(idx_ini:idx_fin) - t_sim(idx_ini);

    % --- omega_m: valor al final del segmento ---
    w_seg = y_lti_nom(idx_ini:idx_fin, 2);
    w_end = w_seg(end);

    % --- i_qs: métricas de transitorio ---
    vqs_prev = u_vqs(max(idx_ini - 1, 1));
    vqs_post = u_vqs(idx_ini);
    hubo_escalon_vqs = abs(vqs_post - vqs_prev) > 1e-9;
    iq_seg = y_lti_nom(idx_ini:idx_fin, 3);
    iq_ini = iq_seg(1);    % valor inicial del segmento
    iq_ss = vqs_post / R_s0;  % valor de régimen teórico del eje q
    delta_iq = iq_ss - iq_ini;  % cambio total teórico

    if hubo_escalon_vqs && abs(delta_iq) > 1e-6
        % Tiempo de crecimiento: 10% a 90% del cambio
        iq_10 = iq_ini + 0.1 * delta_iq;
        iq_90 = iq_ini + 0.9 * delta_iq;
        if delta_iq > 0
            idx_10 = find(iq_seg >= iq_10, 1);
            idx_90 = find(iq_seg >= iq_90, 1);
        else
            idx_10 = find(iq_seg <= iq_10, 1);
            idx_90 = find(iq_seg <= iq_90, 1);
        end
        if ~isempty(idx_10) && ~isempty(idx_90) && idx_90 > idx_10
            tr_iq = (t_seg(idx_90) - t_seg(idx_10)) * 1000;
        else
            tr_iq = NaN;
        end

        % Tiempo de establecimiento: ±1% respecto al valor de régimen teórico
        banda = 0.01 * abs(delta_iq);
        settled = abs(iq_seg - iq_ss) <= banda;
        idx_settled = find(settled, 1, 'first');
        if ~isempty(idx_settled) && all(settled(idx_settled:end))
            ts_iq = t_seg(idx_settled) * 1000;
        else
            ts_iq = NaN;
        end

        % Sobrepico respecto al cambio
        if delta_iq > 0
            Mp_iq = (max(iq_seg) - iq_ss) / abs(delta_iq) * 100;
        else
            Mp_iq = (min(iq_seg) - iq_ss) / delta_iq * 100;  % delta_iq < 0
        end
        if Mp_iq < 0.01; Mp_iq = 0; end
    else
        tr_iq = NaN;
        ts_iq = NaN;
        Mp_iq = NaN;
    end

    metricas(s,:) = [w_end, iq_ss, tr_iq, ts_iq, Mp_iq];

    tr_txt = '--';
    ts_txt = '--';
    mp_txt = '--';
    if ~isnan(tr_iq), tr_txt = sprintf('%.2f', tr_iq); end
    if ~isnan(ts_iq), ts_txt = sprintf('%.2f', ts_iq); end
    if ~isnan(Mp_iq), mp_txt = sprintf('%.2f', Mp_iq); end

    fprintf('  %2d   | %-25s | %10.2f | %10.4f | %10s | %10s | %8s\n', ...
        s, event_names{s}, w_end, iq_ss, tr_txt, ts_txt, mp_txt);

    latex_rows(s) = sprintf([ ...
        '%d & $%s$ & %6.0f & %6.2f & %s & %s & %s \\\\'], ...
        s, strrep(event_names{s}, '->', '\to '), w_end, iq_ss, ...
        strrep(tr_txt, '.', '{,}'), strrep(ts_txt, '.', '{,}'), strrep(mp_txt, '.', '{,}'));
end

fprintf('\n--- Influencia relativa de las acciones externas ---\n');
fprintf('tau_e = L_q/R_s0 = %.2f ms (dinámica eléctrica rápida de i_qs)\n', L_q/R_s0*1000);
fprintf('t_s(1%%) teorico = tau_e * ln(100) = %.2f ms\n', L_q/R_s0 * log(100) * 1000);
fprintf('\n--- Filas LaTeX sugeridas para la tabla ---\n');
for s = 1:10
    fprintf('%s\n', latex_rows(s));
end
fprintf('tau_m = J_eq/b_eq = %.2f s  (dinámica mecánica lenta de omega_m)\n', J_eq/b_eq);
fprintf('Ratio tau_m/tau_e ≈ %.0f:1\n', (J_eq/b_eq)/(L_q/R_s0));
fprintf('=> v_qs* domina la corriente i_qs (respuesta rápida, ~5*tau_e ≈ %.1f ms)\n', 5*L_q/R_s0*1000);
fprintf('=> T_ld perturba omega_m directamente (respuesta lenta, ~5*tau_m ≈ %.1f s)\n', 5*J_eq/b_eq);
fprintf('=> Los escalones de v_qs* producen cambios abruptos en i_qs; los de T_ld no afectan i_qs.\n');

%% =====================================================================
%  SECCIÓN (c): COMPARACIÓN i_ds(t) PARA LOS 3 CASOS
%  =====================================================================
fprintf('\n=== SECCIÓN (c): EFECTO DE i_ds(0) ===\n');

% --- Figura 10a: i_ds(t) modelo NL - 3 condiciones iniciales ---
fig10a = figure('Name', 'i_ds NL - 3 casos', 'Color', 'w', 'Position', [450 100 900 400]);
hold on; grid on; box on;
for k = 1:3
    plot(t_sim, res_nl{k}(:,4), '-', 'Color', case_colors{k}, 'LineWidth', 1.5);
end
ylabel('i_{ds}^r [A]', 'FontSize', 12); xlabel('Tiempo [s]', 'FontSize', 12);
title('Corriente i_{ds}^r(t) - Modelo NL: efecto de la condición inicial', 'FontSize', 13, 'FontWeight', 'bold');
legend(case_names, 'Location', 'northeast');
xlim([0 0.1]);  % Zoom a los primeros 100 ms donde ocurre el decaimiento

exportgraphics(fig10a, 'imagenes/sim_ids_NL_3casos.png', 'Resolution', 300);
fprintf('   -> Exportada: imagenes/sim_ids_NL_3casos.png\n');

% --- Figura 10b: i_ds(t) modelo LTI - 3 condiciones iniciales ---
fig10b = figure('Name', 'i_ds LTI - 3 casos', 'Color', 'w', 'Position', [500 100 900 400]);
hold on; grid on; box on;
for k = 1:3
    plot(t_sim, res_lti{k}(:,4), '-', 'Color', case_colors{k}, 'LineWidth', 1.5);
end
ylabel('i_{ds}^r [A]', 'FontSize', 12); xlabel('Tiempo [s]', 'FontSize', 12);
title('Corriente i_{ds}^r(t) - Modelo LTI: efecto de la condición inicial', 'FontSize', 13, 'FontWeight', 'bold');
legend(case_names, 'Location', 'northeast');
xlim([0 0.1]);  % Zoom a los primeros 100 ms

exportgraphics(fig10b, 'imagenes/sim_ids_LTI_3casos.png', 'Resolution', 300);
fprintf('   -> Exportada: imagenes/sim_ids_LTI_3casos.png\n');

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

% --- Figura 11: Field Forcing/Weakening - i_ds (gráfico principal para informe) ---
fig11 = figure('Name', 'Field Forcing/Weakening - i_ds', 'Color', 'w', 'Position', [500 50 900 450]);
hold on; grid on; box on;

% Caso base NL (sin v_ds*)
plot(t_sim, res_nl{1}(:,4), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Base NL (v_{ds}^*=0)');
plot(t_sim, res_lti{1}(:,4), 'k--', 'LineWidth', 1.2, 'DisplayName', 'Base LTI');

% Field forcing NL y LTI
plot(t_sim, res_ff{1}(:,4), 'r-', 'LineWidth', 1.5, 'DisplayName', 'FF NL (v_{ds}^*=+1.96V)');
plot(t_sim, res_ff_lti{1}(:,4), 'r--', 'LineWidth', 1.2, 'DisplayName', 'FF LTI');

% Field weakening NL y LTI
plot(t_sim, res_ff{2}(:,4), 'b-', 'LineWidth', 1.5, 'DisplayName', 'FW NL (v_{ds}^*=-1.96V)');
plot(t_sim, res_ff_lti{2}(:,4), 'b--', 'LineWidth', 1.2, 'DisplayName', 'FW LTI');

ylabel('i_{ds}^r [A]', 'FontSize', 12); xlabel('Tiempo [s]', 'FontSize', 12);
title({'Field Forcing / Weakening: efecto sobre i_{ds}^r(t)'; ...
       'Comparación Modelo NL vs LTI (i_{ds}^r(0) = 0 A)'}, ...
    'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
xlim([0 t_final]);

exportgraphics(fig11, 'imagenes/sim_ff_ids.png', 'Resolution', 300);
fprintf('   -> Exportada: imagenes/sim_ff_ids.png\n');

% --- Figura 11b: Field Forcing/Weakening - Todos los estados (referencia) ---
fig11b = figure('Name', 'Field Forcing/Weakening - Todos', 'Color', 'w', 'Position', [550 50 1300 800]);

vars_plot = {'\theta_m [rad]', '\omega_m [rad/s]', 'i_{qs} [A]', 'i_{ds} [A]'};

for v = 1:4
    subplot(2,2,v); hold on; grid on;

    plot(t_sim, res_nl{1}(:,v), 'k-', 'LineWidth', 1, 'DisplayName', 'Base NL');
    plot(t_sim, res_lti{1}(:,v), 'k--', 'LineWidth', 1, 'DisplayName', 'Base LTI');
    plot(t_sim, res_ff{1}(:,v), 'r-', 'LineWidth', 1.2, 'DisplayName', 'FF NL');
    plot(t_sim, res_ff_lti{1}(:,v), 'r--', 'LineWidth', 1.2, 'DisplayName', 'FF LTI');
    plot(t_sim, res_ff{2}(:,v), 'b-', 'LineWidth', 1.2, 'DisplayName', 'FW NL');
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
fprintf('  FIGURAS GENERADAS\n');
fprintf('========================================\n');
fprintf('  Fig 1:   senales_entrada (referencia)\n');
fprintf('  Fig 2:   estados_mecanicos_NLvsLTI (referencia)\n');
fprintf('  Fig 3:   corrientes_qd_NLvsLTI (referencia)\n');
fprintf('  Fig 4:   estados_NL_exclusivos (referencia)\n');
fprintf('  Fig 5:   vds_forzada (referencia)\n');
fprintf('  Fig 6:   corrientes_abc_NL (referencia)\n');
fprintf('  Fig 7:   tensiones_abc_NL (referencia)\n');
fprintf('  Fig 8:   torque_vs_velocidad_cuadrantes (referencia)\n');
fprintf('  Fig 8b:  plano_corrientes_ids_vs_iqs (referencia)\n');
fprintf('  Fig 9:   angulo_torque (referencia)\n');
fprintf('  -----------------------------------------------\n');
fprintf('  EXPORTADAS AUTOMÁTICAMENTE a imagenes/:\n');
fprintf('  Fig 10a: sim_ids_NL_3casos.png   (Sub-ítem c)\n');
fprintf('  Fig 10b: sim_ids_LTI_3casos.png  (Sub-ítem c)\n');
fprintf('  Fig 11:  sim_ff_ids.png          (Sub-ítem d)\n');
fprintf('  Fig 11b: field_forcing_todos (referencia)\n');
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
