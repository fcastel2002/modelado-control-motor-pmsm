%% FUNCIONES DE TRANSFERENCIA - MODELO LTI EQUIVALENTE AUMENTADO
%  Consigna 5.1.2.e: FTs desde v_qs(t) y T_l(t) hacia theta_m(t)
%  Proyecto Global Integrador - Automática y Máquinas Eléctricas
%  Autores: Joaquín Calderón - Francisco Castel
%
%  Modelo LTI aumentado (4to orden):
%    x_aug = [theta_m, omega_m, i_qs, i_ds]^T
%    Entradas: v_qs (control), T_ld (perturbación)
%    Salida:   theta_m (posición del motor)

clearvars; clc; close all;

%% 1. Cargar Parámetros Físicos
if exist('parametros_sistema_completo.m', 'file')
    parametros_sistema_completo;
    disp('-> Parámetros físicos cargados exitosamente.');
else
    error('No se encuentra "parametros_sistema_completo.m".');
end

if ~exist('J_eq', 'var') || ~exist('b_eq', 'var')
    error('Faltan variables derivadas (J_eq, b_eq). Revisa el archivo de parámetros.');
end

%% 2. CONSTRUCCIÓN DEL MODELO LTI AUMENTADO (4x4)
% =========================================================================
% Estado inicial considerado: x(0) = 0 (modelo incremental)
%   theta_m(0) = 0, omega_m(0) = 0, i_qs(0) = 0, i_ds(0) = 0
%
% La gravedad se extrae como perturbación externa (no está en A).
% R_s se considera constante a T_s = T_amb (desacoplamiento térmico).
% =========================================================================

fprintf('\n========================================\n');
fprintf('  MODELO LTI EQUIVALENTE AUMENTADO\n');
fprintf('========================================\n');

R_s0 = R_sREF * (1 + 3.9e-3*(29.5-20));  % R_s en operacion T_s=29.5C (media), NO a 20C. = 1.058 Ohm

% --- Matriz de Estado A (4x4) ---
%   x = [theta_m; omega_m; i_qs; i_ds]
A_aug = [
    0,          1,                          0,              0;
    0,         -b_eq/J_eq,     3*P_p*lambda_r/(2*J_eq),    0;
    0,         -P_p*lambda_r/L_q,          -R_s0/L_q,      0;
    0,          0,                          0,             -R_s0/L_d
];

% --- Matrices de Entrada ---
% B_c: entrada de control v_qs(t)
B_c = [0; 0; 1/L_q; 0];

% B_d: entrada de perturbación T_ld(t)
%   T_l(t) = g*k_l*sin(theta_l) + T_ld(t)
%   Pero en el modelo LTI, la gravedad se compensa por feedforward,
%   y T_ld entra reflejado al eje del motor (dividido por r).
B_d = [0; -1/(r*J_eq); 0; 0];

% --- Matriz de Salida C ---
% Salida: y = theta_m
C_aug = [1, 0, 0, 0];

% --- Matriz de transmisión directa D ---
D_c = 0;
D_d = 0;

% Mostrar matrices
fprintf('\nEstado inicial: x(0) = [0; 0; 0; 0]\n');
fprintf('  (Modelo incremental alrededor del equilibrio)\n');
fprintf('\nMatriz A (4x4):\n');
disp(A_aug);
fprintf('Matriz B_c (control, v_qs):\n');
disp(B_c);
fprintf('Matriz B_d (perturbación, T_ld):\n');
disp(B_d);
fprintf('Matriz C:\n');
disp(C_aug);

%% 3. SISTEMAS EN ESPACIO DE ESTADOS
% =========================================================================
% Se crean dos sistemas ss: uno por cada entrada
% =========================================================================

sys_vqs = ss(A_aug, B_c, C_aug, D_c);
sys_vqs.InputName  = 'v_{qs}';
sys_vqs.OutputName = '\theta_m';
sys_vqs.StateName  = {'theta_m', 'omega_m', 'i_qs', 'i_ds'};

sys_Tld = ss(A_aug, B_d, C_aug, D_d);
sys_Tld.InputName  = 'T_{ld}';
sys_Tld.OutputName = '\theta_m';
sys_Tld.StateName  = {'theta_m', 'omega_m', 'i_qs', 'i_ds'};

%% 4. FUNCIONES DE TRANSFERENCIA
% =========================================================================
% G(s) = C * (sI - A)^(-1) * B + D
% =========================================================================

fprintf('\n========================================\n');
fprintf('  FUNCIONES DE TRANSFERENCIA\n');
fprintf('========================================\n');

% --- G1(s): theta_m / v_qs ---
G_vqs = tf(sys_vqs);
G_vqs = minreal(G_vqs);  % Cancelar polos/ceros comunes

fprintf('\n--- G1(s) = theta_m(s) / v_qs(s) ---\n');
disp(G_vqs);

% --- G2(s): theta_m / T_ld ---
G_Tld = tf(sys_Tld);
G_Tld = minreal(G_Tld);

fprintf('\n--- G2(s) = theta_m(s) / T_ld(s) ---\n');
disp(G_Tld);

%% 5. ANÁLISIS DE POLOS, CEROS Y ESTADOS NO OBSERVABLES
% =========================================================================
fprintf('\n========================================\n');
fprintf('  ANÁLISIS DE POLOS Y CEROS\n');
fprintf('========================================\n');

% Autovalores del sistema completo (4x4)
eigenvalues = eig(A_aug);
fprintf('\nAutovalores del sistema aumentado (4x4):\n');
for k = 1:length(eigenvalues)
    if imag(eigenvalues(k)) == 0
        fprintf('  lambda_%d = %.4f rad/s  (tau = %.4f s)\n', ...
            k, real(eigenvalues(k)), -1/real(eigenvalues(k)));
    else
        fprintf('  lambda_%d = %.4f + j%.4f rad/s\n', ...
            k, real(eigenvalues(k)), imag(eigenvalues(k)));
    end
end

% --- Polos y ceros de G1(s) ---
fprintf('\n--- G1(s): theta_m / v_qs ---\n');
[z1, p1, k1] = zpkdata(G_vqs, 'v');
fprintf('Polos:\n');
for k = 1:length(p1)
    fprintf('  p_%d = %.4f', k, real(p1(k)));
    if imag(p1(k)) ~= 0
        fprintf(' + j%.4f', imag(p1(k)));
    end
    fprintf('\n');
end
fprintf('Ceros:\n');
if isempty(z1)
    fprintf('  (Ninguno)\n');
else
    for k = 1:length(z1)
        fprintf('  z_%d = %.4f', k, real(z1(k)));
        if imag(z1(k)) ~= 0
            fprintf(' + j%.4f', imag(z1(k)));
        end
        fprintf('\n');
    end
end
fprintf('Ganancia: K = %.6e\n', k1);

% --- Polos y ceros de G2(s) ---
fprintf('\n--- G2(s): theta_m / T_ld ---\n');
[z2, p2, k2] = zpkdata(G_Tld, 'v');
fprintf('Polos:\n');
for k = 1:length(p2)
    fprintf('  p_%d = %.4f', k, real(p2(k)));
    if imag(p2(k)) ~= 0
        fprintf(' + j%.4f', imag(p2(k)));
    end
    fprintf('\n');
end
fprintf('Ceros:\n');
if isempty(z2)
    fprintf('  (Ninguno)\n');
else
    for k = 1:length(z2)
        fprintf('  z_%d = %.4f', k, real(z2(k)));
        if imag(z2(k)) ~= 0
            fprintf(' + j%.4f', imag(z2(k)));
        end
        fprintf('\n');
    end
end
fprintf('Ganancia: K = %.6e\n', k2);

%% 6. ESTADO INTERNO NO OBSERVABLE: i_ds
% =========================================================================
% El estado i_ds (4to estado) NO aparece en las FTs.
% Esto se demuestra analizando la estructura de A y B:
%   - i_ds no es excitado por ninguna entrada (B_c(4)=0, B_d(4)=0)
%   - i_ds no acopla con ningún otro estado (fila 4 de A solo tiene -Rs0/Ld)
%   - i_ds no es medido por C (C(4)=0)
% Por tanto, i_ds es un ESTADO NO CONTROLABLE Y NO OBSERVABLE desde
% las entradas/salida consideradas.
% =========================================================================

fprintf('\n========================================\n');
fprintf('  ESTADOS INTERNOS NO REFLEJADOS EN LAS FTs\n');
fprintf('========================================\n');

fprintf('\nEl estado i_ds^r (4to estado) NO aporta ni se refleja en las FTs.\n');
fprintf('Razones:\n');
fprintf('  1. Ninguna entrada excita a i_ds: B_c(4)=0, B_d(4)=0\n');
fprintf('  2. i_ds no acopla con los demás estados (bloque diagonal)\n');
fprintf('  3. La salida theta_m no depende de i_ds: C(4)=0\n');
fprintf('\nEl polo lambda = -R_s0/L_d = %.4f rad/s asociado a i_ds\n', -R_s0/L_d);
fprintf('aparece en el denominador del ss pero se CANCELA en la FT (minreal).\n');
fprintf('\nInterpretación física:\n');
fprintf('  La dinámica residual i_ds(t) = i_ds(0)*exp(-R_s0/L_d * t)\n');
fprintf('  es una respuesta libre que decae exponencialmente.\n');
fprintf('  Solo existe si i_ds(0) ≠ 0 (condición inicial no nula).\n');
fprintf('  El modelo LTI NOMINAL (3 estados) es suficiente para las FTs.\n');

% Verificación formal: Controlabilidad y Observabilidad
fprintf('\n--- Verificación formal ---\n');

% Controlabilidad desde v_qs
Mc_vqs = ctrb(A_aug, B_c);
rank_Mc_vqs = rank(Mc_vqs);
fprintf('Rango de Controlabilidad (v_qs):  %d / %d\n', rank_Mc_vqs, size(A_aug,1));

% Controlabilidad desde T_ld
Mc_Tld = ctrb(A_aug, B_d);
rank_Mc_Tld = rank(Mc_Tld);
fprintf('Rango de Controlabilidad (T_ld):  %d / %d\n', rank_Mc_Tld, size(A_aug,1));

% Observabilidad
Mo = obsv(A_aug, C_aug);
rank_Mo = rank(Mo);
fprintf('Rango de Observabilidad:          %d / %d\n', rank_Mo, size(A_aug,1));

if rank_Mc_vqs < size(A_aug,1)
    fprintf('\n=> El sistema NO es completamente controlable desde v_qs.\n');
    fprintf('   El estado no controlable es i_ds (modo autónomo).\n');
end
if rank_Mo < size(A_aug,1)
    fprintf('=> El sistema NO es completamente observable desde theta_m.\n');
    fprintf('   El estado no observable es i_ds.\n');
end

%% 7. COMPARACIÓN: MODELO AUMENTADO (4) vs NOMINAL (3)
% =========================================================================
% Se verifica que las FTs del modelo aumentado coinciden con las del
% modelo nominal de 3 estados (sin i_ds).
% =========================================================================

fprintf('\n========================================\n');
fprintf('  COMPARACIÓN: AUMENTADO (4) vs NOMINAL (3)\n');
fprintf('========================================\n');

% Modelo nominal (3x3): x = [theta_m; omega_m; i_qs]
A_nom = A_aug(1:3, 1:3);
B_c_nom = B_c(1:3);
B_d_nom = B_d(1:3);
C_nom = [1, 0, 0];

sys_vqs_nom = ss(A_nom, B_c_nom, C_nom, 0);
sys_Tld_nom = ss(A_nom, B_d_nom, C_nom, 0);

G_vqs_nom = tf(sys_vqs_nom);
G_Tld_nom = tf(sys_Tld_nom);

fprintf('\nG1_nominal(s) = theta_m / v_qs (3 estados):\n');
disp(G_vqs_nom);

fprintf('G2_nominal(s) = theta_m / T_ld (3 estados):\n');
disp(G_Tld_nom);

% Verificación numérica
w_test = logspace(-2, 5, 1000);
[mag1_aug, ~] = bode(G_vqs, w_test);
[mag1_nom, ~] = bode(G_vqs_nom, w_test);
error_max_G1 = max(abs(squeeze(mag1_aug) - squeeze(mag1_nom)));

[mag2_aug, ~] = bode(G_Tld, w_test);
[mag2_nom, ~] = bode(G_Tld_nom, w_test);
error_max_G2 = max(abs(squeeze(mag2_aug) - squeeze(mag2_nom)));

fprintf('Error máximo en magnitud (G1): %.2e\n', error_max_G1);
fprintf('Error máximo en magnitud (G2): %.2e\n', error_max_G2);

if error_max_G1 < 1e-10 && error_max_G2 < 1e-10
    fprintf('=> CONFIRMADO: Las FTs del modelo aumentado y nominal son idénticas.\n');
    fprintf('   El estado i_ds no aporta a las funciones de transferencia.\n');
end

%% 8. EXPRESIONES SIMBÓLICAS DE LAS FTs
% =========================================================================
% Se derivan las FTs en forma simbólica para incluir en el informe.
% =========================================================================

fprintf('\n========================================\n');
fprintf('  EXPRESIONES SIMBÓLICAS\n');
fprintf('========================================\n');

syms s Jeq beq Pp lam Lq Ld Rs0 rr positive

A_sym = [0,  1,              0;
         0, -beq/Jeq,  3*Pp*lam/(2*Jeq);
         0, -Pp*lam/Lq,    -Rs0/Lq];

B_c_sym = [0; 0; 1/Lq];
B_d_sym = [0; -1/(rr*Jeq); 0];
C_sym   = [1, 0, 0];

% G1(s) = C * (sI - A)^-1 * B_c
resolvent = inv(s*eye(3) - A_sym);
G1_sym = simplify(C_sym * resolvent * B_c_sym);
G2_sym = simplify(C_sym * resolvent * B_d_sym);

fprintf('\nG1(s) = theta_m / v_qs =\n');
[num1, den1] = numden(G1_sym);
fprintf('  Numerador:   '); disp(num1);
fprintf('  Denominador: '); disp(den1);

fprintf('\nG2(s) = theta_m / T_ld =\n');
[num2, den2] = numden(G2_sym);
fprintf('  Numerador:   '); disp(num2);
fprintf('  Denominador: '); disp(den2);

% Denominador común (polinomio característico)
char_poly = simplify(det(s*eye(3) - A_sym));
fprintf('\nPolinomio Característico del modelo nominal:\n');
fprintf('  det(sI - A) = '); disp(char_poly);

%% 9. VISUALIZACIÓN: MAPA DE POLOS-CEROS Y BODE
% =========================================================================

figure('Name', 'Funciones de Transferencia - LTI Aumentado', ...
       'Color', 'w', 'Position', [50, 50, 1400, 800]);

% --- Mapa de Polos y Ceros ---
subplot(2, 2, 1);
pzmap(G_vqs);
hold on;
% Marcar el polo cancelado de i_ds
plot(real(-R_s0/L_d), 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
title('Polos y Ceros: G_1(s) = \theta_m / v_{qs}');
legend('Polos/Ceros activos', 'Polo cancelado (i_{ds})', 'Location', 'best');
grid on;

subplot(2, 2, 2);
pzmap(G_Tld);
hold on;
plot(real(-R_s0/L_d), 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
title('Polos y Ceros: G_2(s) = \theta_m / T_{ld}');
legend('Polos/Ceros activos', 'Polo cancelado (i_{ds})', 'Location', 'best');
grid on;

% --- Diagramas de Bode ---
subplot(2, 2, 3);
bode(G_vqs, w_test);
title('Bode: G_1(s) = \theta_m / v_{qs}');
grid on;

subplot(2, 2, 4);
bode(G_Tld, w_test);
title('Bode: G_2(s) = \theta_m / T_{ld}');
grid on;

sgtitle('Consigna 5.1.2.e: Funciones de Transferencia del Modelo LTI Aumentado', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% 10. RESPUESTA AL ESCALÓN
% =========================================================================

figure('Name', 'Respuesta al Escalón - LTI Aumentado', ...
       'Color', 'w', 'Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
step(G_vqs);
title('Respuesta al escalón: \theta_m ante escalón en v_{qs}');
xlabel('Tiempo [s]'); ylabel('\theta_m [rad]');
grid on;

subplot(1, 2, 2);
step(G_Tld);
title('Respuesta al escalón: \theta_m ante escalón en T_{ld}');
xlabel('Tiempo [s]'); ylabel('\theta_m [rad]');
grid on;

sgtitle('Respuestas al Escalón Unitario', 'FontSize', 14);

%% 11. RESUMEN FINAL
% =========================================================================

fprintf('\n========================================\n');
fprintf('  RESUMEN - CONSIGNA 5.1.2.e\n');
fprintf('========================================\n');
fprintf('\n1) ESTADO INICIAL CONSIDERADO:\n');
fprintf('   x(0) = [theta_m; omega_m; i_qs; i_ds] = [0; 0; 0; 0]\n');
fprintf('   (Modelo incremental alrededor del punto de equilibrio)\n');
fprintf('   La gravedad se trata como perturbación externa compensada\n');
fprintf('   por feedforward => no aparece en A.\n');

fprintf('\n2) FUNCIONES DE TRANSFERENCIA (Orden 3, modelo nominal):\n');
fprintf('   G1(s) = theta_m(s) / v_qs(s)\n');
fprintf('   G2(s) = theta_m(s) / T_ld(s)\n');
fprintf('   Ambas comparten el mismo denominador (polinomio característico).\n');

fprintf('\n3) ESTADO INTERNO QUE NO APORTA:\n');
fprintf('   El estado i_ds^r (corriente de eje directo) NO se refleja\n');
fprintf('   en las funciones de transferencia.\n');
fprintf('   - No es controlable desde v_qs ni desde T_ld.\n');
fprintf('   - No es observable desde theta_m.\n');
fprintf('   - Su polo (-R_s0/L_d = %.2f rad/s) se cancela al\n', -R_s0/L_d);
fprintf('     pasar de representación ss a tf.\n');
fprintf('   - Físicamente: solo responde a su condición inicial\n');
fprintf('     (dinámica residual libre, decae exponencialmente).\n');

fprintf('\n========================================\n');
fprintf('  FIN DEL ANÁLISIS\n');
fprintf('========================================\n');
