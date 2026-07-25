%% ANÁLISIS DE ESTABILIDAD A LAZO ABIERTO - MODELO LTI AUMENTADO
%  Consigna 5.1.3: Polos, ceros, estabilidad, migración de propiedades
%  Proyecto Global Integrador - Automática y Máquinas Eléctricas
%  Autores: Joaquín Calderón - Francisco Castel

clearvars; clc; close all;

%% 1. Cargar Parámetros Físicos
if exist('parametros_sistema_completo.m', 'file')
    parametros_sistema_completo;
    disp('-> Parámetros físicos cargados exitosamente.');
else
    error('No se encuentra "parametros_sistema_completo.m".');
end

R_s0 = R_sREF * (1 + 3.9e-3*(29.5-20));   % R_s en el punto de operacion termico T_s=29.5C (media de la sim 5.1.6), NO a 20C. = 1.058 Ohm

%% 2. MODELO LTI AUMENTADO (4x4) - NOMINAL
% Estado: x = [theta_m; omega_m; i_qs; i_ds]
% Gravedad extraída, R_s constante.

A_aug = [
    0,    1,                        0,           0;
    0,   -b_eq/J_eq,  3*P_p*lambda_r/(2*J_eq),  0;
    0,   -P_p*lambda_r/L_q,        -R_s0/L_q,   0;
    0,    0,                         0,          -R_s0/L_d
];

B_c = [0; 0; 1/L_q; 0];           % Entrada: v_qs
B_d = [0; -1/(r*J_eq); 0; 0];     % Perturbación: T_ld
C_aug = [1, 0, 0, 0];             % Salida: theta_m

% Modelo nominal (3x3) para FTs
A_nom = A_aug(1:3, 1:3);
B_c_nom = B_c(1:3);
B_d_nom = B_d(1:3);
C_nom = [1, 0, 0];

% Funciones de transferencia
G1 = minreal(tf(ss(A_nom, B_c_nom, C_nom, 0)));  % theta_m / v_qs
G2 = minreal(tf(ss(A_nom, B_d_nom, C_nom, 0)));  % theta_m / T_ld

%% 3. EVALUACIÓN NUMÉRICA DE POLOS Y CEROS
fprintf('\n============================================\n');
fprintf('  POLOS Y CEROS - VALORES NOMINALES\n');
fprintf('============================================\n');

% Autovalores del sistema aumentado
eigs_aug = eig(A_aug);
fprintf('\nAutovalores del sistema aumentado (4x4):\n');
for k = 1:length(eigs_aug)
    if imag(eigs_aug(k)) == 0
        fprintf('  lambda_%d = %.4f rad/s\n', k, real(eigs_aug(k)));
    else
        fprintf('  lambda_%d = %.4f %+.4fj rad/s\n', k, real(eigs_aug(k)), imag(eigs_aug(k)));
    end
end

% Polos y ceros de G1
[z1, p1, k1] = zpkdata(G1, 'v');
fprintf('\n--- G1(s) = theta_m / v_qs ---\n');
fprintf('Polos: '); disp(p1(:).');
fprintf('Ceros: ');
if isempty(z1); fprintf('(ninguno)\n'); else; disp(z1(:).'); end

% Polos y ceros de G2
[z2, p2, k2] = zpkdata(G2, 'v');
fprintf('\n--- G2(s) = theta_m / T_ld ---\n');
fprintf('Polos: '); disp(p2(:).');
fprintf('Ceros: '); disp(z2(:).');

% Frecuencia natural y amortiguamiento (polos electromecánicos)
p_em = p1(imag(p1) > 0);  % polo con parte imaginaria positiva
wn = abs(p_em);
zeta = -real(p_em)/wn;
wd = imag(p_em);
tau_em = -1/real(p_em);

fprintf('\n--- Parámetros del par electromecánico ---\n');
fprintf('  omega_n  = %.2f rad/s (%.2f Hz)\n', wn, wn/(2*pi));
fprintf('  zeta     = %.4f\n', zeta);
fprintf('  omega_d  = %.2f rad/s (%.2f Hz)\n', wd, wd/(2*pi));
fprintf('  tau_em   = %.4f s (%.2f ms)\n', tau_em, tau_em*1000);

% Polo residual i_ds
polo_ids = -R_s0/L_d;
tau_ids = -1/polo_ids;
fprintf('\n--- Polo residual (i_ds) ---\n');
fprintf('  lambda_4 = %.2f rad/s\n', polo_ids);
fprintf('  tau_d    = %.4f s (%.2f ms)\n', tau_ids, tau_ids*1000);

% Cero de G2
z_G2 = -R_s0/L_q;
fprintf('\n--- Cero de G2 ---\n');
fprintf('  z_1      = %.2f rad/s\n', z_G2);
fprintf('  tau_z    = %.4f s (%.2f ms)\n', -1/z_G2, -1000/z_G2);

%% 4. FIGURA 1: MAPA DE POLOS Y CEROS A LAZO ABIERTO
fig1 = figure('Name', 'Mapa de Polos y Ceros - Lazo Abierto', ...
    'Color', 'w', 'Position', [50, 100, 1400, 550]);

% --- G1: theta_m / v_qs ---
subplot(1, 2, 1); hold on; grid on; box on;

% Polos del modelo nominal (aparecen en la FT) - cuadrados solidos
plot(real(p1), imag(p1), 's', 'MarkerSize', 10, 'MarkerFaceColor', [0 0.447 0.741], 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

% Polo cancelado de i_ds (no aparece en FT) - rombo solido rojo
plot(real(polo_ids), 0, 'd', 'MarkerSize', 11, 'MarkerFaceColor', [0.635 0.078 0.184], 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

% Ejes
xline(0, 'k-', 'LineWidth', 0.5);
yline(0, 'k-', 'LineWidth', 0.5);

% Anotaciones
text(real(p1(1))+3, imag(p1(1))+10, ...
    sprintf('%.1f%+.1fj', real(p1(1)), imag(p1(1))), 'FontSize', 8);
text(5, -5, 's=0', 'FontSize', 9, 'Color', 'b');
text(polo_ids-10, 12, sprintf('%.1f (cancelado)', polo_ids), ...
    'FontSize', 8, 'Color', 'r');

title('$G_1(s) = \theta_m / v_{qs}$', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');
legend('Polos activos', 'Polo cancelado (i_{ds})', 'Location', 'northwest');

% --- G2: theta_m / T_ld ---
subplot(1, 2, 2); hold on; grid on; box on;

% Polos - cuadrados solidos
plot(real(p2), imag(p2), 's', 'MarkerSize', 10, 'MarkerFaceColor', [0 0.447 0.741], 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
% Cero - circulo vacio
plot(real(z2), imag(z2), 'o', 'MarkerSize', 11, 'MarkerEdgeColor', [0 0.5 0], 'MarkerFaceColor', 'none', 'LineWidth', 1.8);
% Polo cancelado - rombo solido rojo
plot(real(polo_ids), 0, 'd', 'MarkerSize', 11, 'MarkerFaceColor', [0.635 0.078 0.184], 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

xline(0, 'k-', 'LineWidth', 0.5);
yline(0, 'k-', 'LineWidth', 0.5);

text(real(z2)+3, 10, sprintf('z=%.1f', real(z2)), 'FontSize', 8, 'Color', [0 0.5 0]);

title('$G_2(s) = \theta_m / T_{ld}$', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');
legend('Polos activos', 'Cero', 'Polo cancelado (i_{ds})', 'Location', 'northwest');

sgtitle('Mapa de Polos y Ceros a Lazo Abierto (Modelo LTI Aumentado)', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Exportar
exportgraphics(fig1, 'imagenes/polos_ceros_LA.png', 'Resolution', 200);
disp('-> Figura 1 generada y exportada: Mapa de polos y ceros');

%% 5. MIGRACIÓN DE PROPIEDADES ANTE VARIACIÓN DE Rs
% Rango: [0.8*R_sREF, 1.2*R_sREF] (variación por temperatura)

fprintf('\n============================================\n');
fprintf('  MIGRACIÓN DE POLOS - VARIACIÓN DE Rs\n');
fprintf('============================================\n');

N_Rs = 13;
Rs_vec = linspace(0.8*R_sREF, 1.2*R_sREF, N_Rs);

% Almacenamiento
polos_Rs  = zeros(3, N_Rs);   % 3 polos del modelo nominal
cero_Rs   = zeros(1, N_Rs);   % cero de G2
wn_Rs     = zeros(1, N_Rs);
zeta_Rs   = zeros(1, N_Rs);

fprintf('\n  Rs [Ohm]   |  Polos (s2,3)                  |  Cero G2    |  wn [rad/s] |  zeta\n');
fprintf('  -----------|--------------------------------|-------------|-------------|-------\n');

for i = 1:N_Rs
    Rs_i = Rs_vec(i);

    A_i = [0,  1,                         0;
           0, -b_eq/J_eq,  3*P_p*lambda_r/(2*J_eq);
           0, -P_p*lambda_r/L_q,         -Rs_i/L_q];

    p_i = eig(A_i);
    polos_Rs(:, i) = p_i;

    % Cero de G2: -Rs/Lq
    cero_Rs(i) = -Rs_i/L_q;

    % wn y zeta del par complejo
    p_cplx = p_i(imag(p_i) > 0);
    if ~isempty(p_cplx)
        wn_Rs(i) = abs(p_cplx);
        zeta_Rs(i) = -real(p_cplx)/wn_Rs(i);
    end

    fprintf('  %.3f      |  %.1f ± j%.1f         |  %.2f   |  %.1f       |  %.3f\n', ...
        Rs_i, real(p_cplx), imag(p_cplx), cero_Rs(i), wn_Rs(i), zeta_Rs(i));
end

% --- FIGURA 2: Migración de polos variando Rs ---
fig2 = figure('Name', 'Migración de Polos - Variación de Rs', ...
    'Color', 'w', 'Position', [50, 50, 1400, 600]);

subplot(1, 2, 1); hold on; grid on; box on;

% Colormap de azul (bajo Rs) a rojo (alto Rs)
cmap = [linspace(0,1,N_Rs)', zeros(N_Rs,1), linspace(1,0,N_Rs)'];

for i = 1:N_Rs
    p_i = polos_Rs(:, i);
    plot(real(p_i), imag(p_i), 's', 'MarkerFaceColor', cmap(i,:), ...
        'MarkerEdgeColor', cmap(i,:), 'MarkerSize', 8);
    % Cero de G2 - circulo vacio
    plot(cero_Rs(i), 0, 'o', 'MarkerEdgeColor', cmap(i,:), ...
        'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.5);
end

% Flecha indicando dirección de aumento de Rs
annotation('textarrow', [0.40 0.27], [0.75 0.65], ...
    'String', 'R_s \uparrow', 'FontSize', 11);

xline(0, 'k-', 'LineWidth', 0.5);
yline(0, 'k-', 'LineWidth', 0.5);
title('Migración de polos y ceros variando R_s', 'FontSize', 12);
xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');

cb = colorbar;
colormap(cmap);
clim([Rs_vec(1), Rs_vec(end)]);
cb.Label.String = 'R_s [\Omega]';

% --- Subplot derecho: wn y zeta vs Rs ---
subplot(1, 2, 2); hold on; grid on; box on;

yyaxis left
plot(Rs_vec, wn_Rs, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
ylabel('\omega_n [rad/s]', 'Color', 'b');
ylim([min(wn_Rs)-1, max(wn_Rs)+1]);

yyaxis right
plot(Rs_vec, zeta_Rs, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 5);
ylabel('\zeta', 'Color', 'r');

xlabel('R_s [\Omega]');
title('Frecuencia natural y amortiguamiento vs R_s', 'FontSize', 12);

sgtitle('Migración de Propiedades ante Variación de R_s', ...
    'FontSize', 14, 'FontWeight', 'bold');

exportgraphics(fig2, 'imagenes/migracion_Rs.png', 'Resolution', 200);
disp('-> Figura 2 generada y exportada: Migración Rs');

%% 6. MIGRACIÓN DE PROPIEDADES ANTE VARIACIÓN DE CARGA (m_l)
% Consigna ítem 4.b: m_l varía de 0 a 1.5 kg

fprintf('\n============================================\n');
fprintf('  MIGRACIÓN DE POLOS - VARIACIÓN DE m_l\n');
fprintf('============================================\n');

N_ml = 16;
ml_vec = linspace(0, 1.5, N_ml);

polos_ml = zeros(3, N_ml);
cero_ml  = zeros(1, N_ml);
wn_ml    = zeros(1, N_ml);
zeta_ml  = zeros(1, N_ml);
Jeq_ml   = zeros(1, N_ml);

fprintf('\n  m_l [kg]  |  J_eq [kg·m²]  |  Polos (s2,3)                  |  wn [rad/s] |  zeta\n');
fprintf('  ----------|----------------|--------------------------------|-------------|-------\n');

for i = 1:N_ml
    ml_i = ml_vec(i);

    % Recalcular parámetros dependientes de m_l
    J_l_i  = (m*l_cm^2 + J_cm) + ml_i * l_l^2;
    J_eq_i = J_m + J_l_i / r^2;
    % b_eq no depende de m_l (solo de b_m y b_l)

    Jeq_ml(i) = J_eq_i;

    A_i = [0,  1,                          0;
           0, -b_eq/J_eq_i,  3*P_p*lambda_r/(2*J_eq_i);
           0, -P_p*lambda_r/L_q,          -R_s0/L_q];

    p_i = eig(A_i);
    polos_ml(:, i) = p_i;
    cero_ml(i) = -R_s0/L_q;  % No depende de m_l

    p_cplx = p_i(imag(p_i) > 0);
    if ~isempty(p_cplx)
        wn_ml(i) = abs(p_cplx);
        zeta_ml(i) = -real(p_cplx)/wn_ml(i);
    else
        % Polos reales (sobreamortiguado)
        p_real = sort(p_i(p_i ~= 0));
        wn_ml(i) = sqrt(prod(abs(p_real)));
        zeta_ml(i) = -sum(p_real)/(2*wn_ml(i));
    end

    fprintf('  %.2f      |  %.4e    |  %.1f ± j%.1f         |  %.1f       |  %.3f\n', ...
        ml_i, J_eq_i, real(p_cplx), imag(p_cplx), wn_ml(i), zeta_ml(i));
end

% --- FIGURA 3: Migración de polos variando m_l ---
fig3 = figure('Name', 'Migración de Polos - Variación de m_l', ...
    'Color', 'w', 'Position', [50, 50, 1400, 600]);

subplot(1, 2, 1); hold on; grid on; box on;

cmap_ml = [linspace(0,0.8,N_ml)', linspace(0.6,0,N_ml)', linspace(0,0.8,N_ml)'];

for i = 1:N_ml
    p_i = polos_ml(:, i);
    plot(real(p_i), imag(p_i), 's', 'MarkerFaceColor', cmap_ml(i,:), ...
        'MarkerEdgeColor', cmap_ml(i,:), 'MarkerSize', 8);
end

% Marcar caso nominal (m_l = 0) - rombo solido
p_nom = polos_ml(:, 1);
plot(real(p_nom), imag(p_nom), 'd', 'MarkerSize', 14, 'MarkerFaceColor', [0 0.447 0.741], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

% Marcar peor caso (m_l = 1.5) - cuadrado solido
p_worst = polos_ml(:, end);
plot(real(p_worst), imag(p_worst), 's', 'MarkerSize', 14, 'MarkerFaceColor', [0.635 0.078 0.184], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

xline(0, 'k-', 'LineWidth', 0.5);
yline(0, 'k-', 'LineWidth', 0.5);

title('Migración de polos variando m_l', 'FontSize', 12);
xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');
legend('', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
    'm_l = 0 kg (nominal)', 'm_l = 1.5 kg (peor caso)', ...
    'Location', 'northwest');

cb = colorbar;
colormap(gca, cmap_ml);
clim([0, 1.5]);
cb.Label.String = 'm_l [kg]';

% --- Subplot derecho: wn y zeta vs m_l ---
subplot(1, 2, 2); hold on; grid on; box on;

yyaxis left
plot(ml_vec, wn_ml, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
ylabel('\omega_n [rad/s]', 'Color', 'b');

yyaxis right
plot(ml_vec, zeta_ml, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 5);
ylabel('\zeta', 'Color', 'r');

xlabel('m_l [kg]');
title('Frecuencia natural y amortiguamiento vs m_l', 'FontSize', 12);
xline(0, 'k--', 'Nominal');
xline(1.5, 'k--', 'Peor caso');

sgtitle('Migración de Propiedades ante Variación de Carga (m_l)', ...
    'FontSize', 14, 'FontWeight', 'bold');

exportgraphics(fig3, 'imagenes/migracion_ml.png', 'Resolution', 200);
disp('-> Figura 3 generada y exportada: Migración m_l');

%% 7. DIAGRAMAS DE BODE
fig4 = figure('Name', 'Diagramas de Bode - Lazo Abierto', ...
    'Color', 'w', 'Position', [100, 100, 1200, 600]);

w_bode = logspace(-1, 5, 2000);

subplot(1, 2, 1);
bode(G1, w_bode);
title('$G_1(s) = \theta_m / v_{qs}$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

subplot(1, 2, 2);
bode(G2, w_bode);
title('$G_2(s) = \theta_m / T_{ld}$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

sgtitle('Diagramas de Bode a Lazo Abierto', 'FontSize', 14, 'FontWeight', 'bold');

disp('-> Figura 4 generada: Bode');

%% 8. TABLA RESUMEN PARA CONSOLA (copiar al informe)
fprintf('\n============================================\n');
fprintf('  TABLA RESUMEN - MIGRACIÓN Rs\n');
fprintf('============================================\n');
fprintf('Rs [Ohm] & Polos & Cero & wn [rad/s] & zeta \\\\\n');
fprintf('\\hline\n');
for i = 1:N_Rs
    p_cplx = polos_Rs(imag(polos_Rs(:,i)) > 0, i);
    fprintf('%.3f & $%.1f \\pm j%.1f$ & $%.2f$ & $%.1f$ & $%.3f$ \\\\\n', ...
        Rs_vec(i), real(p_cplx), imag(p_cplx), cero_Rs(i), wn_Rs(i), zeta_Rs(i));
end

fprintf('\n============================================\n');
fprintf('  TABLA RESUMEN - MIGRACIÓN m_l\n');
fprintf('============================================\n');
fprintf('m_l [kg] & J_eq [kg·m²] & Polos & wn [rad/s] & zeta \\\\\n');
fprintf('\\hline\n');
for i = 1:N_ml
    p_cplx = polos_ml(imag(polos_ml(:,i)) > 0, i);
    if ~isempty(p_cplx)
        fprintf('%.2f & %.4e & $%.1f \\pm j%.1f$ & $%.1f$ & $%.3f$ \\\\\n', ...
            ml_vec(i), Jeq_ml(i), real(p_cplx), imag(p_cplx), wn_ml(i), zeta_ml(i));
    end
end

fprintf('\n-> Script finalizado. Figuras guardadas en Proyecto/imagenes/\n');
