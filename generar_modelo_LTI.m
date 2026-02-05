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


%% =========================================================================
%  CASO C: MODELO NOMINAL DE DISEÑO (SIN GRAVEDAD)
%  Este es el modelo que "verá" el PID asumiendo Feedforward perfecto.
% =========================================================================

% 1. Tomamos la matriz del PEOR CASO DE INERCIA (90 grados)
%    (Esto asegura que usemos J_eq máximo, lo cual es conservador y seguro)
A_design = A_90; 

% 2. ELIMINACIÓN DE LA GRAVEDAD (Linearization by Feedback)
%    La gravedad es el término A(2,1) que vincula Posición -> Velocidad.
%    Al hacerlo cero, desacoplamos la posición de la aceleración.
A_design(2,1) = 0; 

% 3. CREACIÓN DEL OBJETO LTI
%    Nota: La matriz B es la misma (la potencia del motor no cambia)
sys_design = ss(A_design, B_90, eye(6), 0);

% 4. MATRICES PARA SIMULINK
A_design_sim = A_design;
B_design_sim = B_90;

disp(' ');
disp('=== MODELO DE DISEÑO CREADO (sys_design) ===');
disp('-> Se ha eliminado el término de rigidez gravitacional (A21=0).');
disp('-> Este modelo representa una Inercia Pura + Fricción.');
disp('-> Listo para usar en Simulink y PID Tuner.');

%% 4. VISUALIZACIÓN TOTAL (Matriz Panorámica vs Zoom)
%  Generamos 6 gráficos: 
%  - Arriba: Escala completa para ver polos Eléctricos lejanos.
%  - Abajo: Zoom para ver polos Mecánicos cercanos.

figure('Name', 'Analisis de Polos: Macro vs Micro', 'Color', 'w', 'Position', [50, 50, 1400, 700]);

% --- DATOS ---
p0   = pole(sys_0);      % Azul
p90  = pole(sys_90);     % Rojo
pDes = pole(sys_design); % Magenta

% Estilo
msize = 8; lw = 1.5;

% Límites para el ZOOM (Fila de abajo)
zoom_x = [-4, 1];   % Eje Real cerca del origen
zoom_y = [-5, 5];   % Eje Imag (Oscilación)

% =========================================================================
% FILA 1: VISTA PANORÁMICA (Para ver los Polos Eléctricos Lejanos)
% =========================================================================

% CASO 0° (PANORAMA)
subplot(2, 3, 1); hold on; grid on; box on;
plot(real(p0), imag(p0), 'bx', 'MarkerSize', msize, 'LineWidth', lw);
title('0° (VISTA COMPLETA)', 'FontWeight', 'bold');
xlabel('Eje Real (Grandes valores)'); ylabel('Eje Imag');
xlim('auto'); ylim('auto'); % Deja que Matlab muestre todo el rango

% CASO 90° (PANORAMA)
subplot(2, 3, 2); hold on; grid on; box on;
plot(real(p90), imag(p90), 'ro', 'MarkerSize', msize, 'LineWidth', lw);
title('90° (VISTA COMPLETA)', 'FontWeight', 'bold');
xlim('auto'); ylim('auto');

% CASO DISEÑO (PANORAMA)
subplot(2, 3, 3); hold on; grid on; box on;
plot(real(pDes), imag(pDes), 'm+', 'MarkerSize', msize+2, 'LineWidth', lw);
title('DISEÑO (VISTA COMPLETA)', 'FontWeight', 'bold');
xlim('auto'); ylim('auto');

% =========================================================================
% FILA 2: VISTA ZOOM (Para ver la Dinámica Mecánica)
% =========================================================================

% CASO 0° (ZOOM)
subplot(2, 3, 4); hold on; grid on; box on;
plot(real(p0), imag(p0), 'bx', 'MarkerSize', msize+2, 'LineWidth', lw);
xline(0, 'k--'); yline(0, 'k--');
title('0° (ZOOM MECÁNICO)', 'Color', 'b');
xlabel('Eje Real (Cerca de 0)'); ylabel('Eje Imag');
xlim(zoom_x); ylim(zoom_y); % FORZAMOS ZOOM
text(-2, 2, 'Oscilación Visible', 'FontSize', 8, 'BackgroundColor', 'w');

% CASO 90° (ZOOM)
subplot(2, 3, 5); hold on; grid on; box on;
plot(real(p90), imag(p90), 'ro', 'MarkerSize', msize+2, 'LineWidth', lw);
xline(0, 'k--'); yline(0, 'k--');
title('90° (ZOOM MECÁNICO)', 'Color', 'r');
xlabel('Eje Real'); 
xlim(zoom_x); ylim(zoom_y); 

% CASO DISEÑO (ZOOM)
subplot(2, 3, 6); hold on; grid on; box on;
plot(real(pDes), imag(pDes), 'm+', 'MarkerSize', msize+4, 'LineWidth', lw);
xline(0, 'k--'); yline(0, 'k--');
title('DISEÑO (ZOOM MECÁNICO)', 'Color', 'm');
xlabel('Eje Real'); 
xlim(zoom_x); ylim(zoom_y);
text(-2, 1, 'Integrador (0,0)', 'FontSize', 8, 'BackgroundColor', 'w');

sgtitle('MAPA COMPLETO DE POLOS: Arriba (Eléctricos) vs. Abajo (Mecánicos)', 'FontSize', 16);

disp(' ');
disp('--- ANÁLISIS GENERADO ---');
disp('-> FILA SUPERIOR: Mira el Eje X negativo. Ahí verás los polos eléctricos lejanos.');
disp('   (Nota: Los polos mecánicos se ven amontonados en el cero, es normal).');
disp('-> FILA INFERIOR: Es un microscopio al origen. Ahí verás la estabilidad real.');

%% 5. CONFIGURACIÓN DE VALIDACIÓN Y ENTORNO DE SIMULACIÓN
%  Selección del escenario de validación para el diagrama de Simulink.
%  Permite alternar entre la validación del punto de operación (física)
%  y la validación del modelo nominal de diseño (control).

% --- Gestión de variables para Simulink (Corrección de Tipos) ---
if exist('parametros_sistema_completo.mlx', 'file')
    clear T_amb; T_amb = 20; % Asegura valor numérico (double) en vez de simbólico
end
T_ld = 0; % Inicializamos perturbación externa en cero

fprintf('\n================================================================\n');
fprintf('   SELECCIÓN DE ESCENARIO DE VALIDACIÓN (MODELO VS PLANTA)\n');
fprintf('================================================================\n');
fprintf(' [1] ESCENARIO A: Modelo Completo (Gravedad como Dinámica Interna)\n');
fprintf('     -> Punto de Operación: 90 grados (Equilibrio Inestable).\n');
fprintf('     -> Objetivo: Validar linealización de la física (Jacobianos).\n');
fprintf('     -> Resultado esperado: Divergencia leve por no-linealidad.\n');
fprintf('\n');
fprintf(' [2] ESCENARIO B: Modelo Nominal (Gravedad como Perturbación)\n');
fprintf('     -> Punto de Operación: Dinámica Incremental (Inercial).\n');
fprintf('     -> Objetivo: Validar modelo para diseño de control (Feedforward).\n');
fprintf('     -> Resultado esperado: Coincidencia exacta (Rampas superpuestas).\n');
fprintf('----------------------------------------------------------------\n');

% Validación de entrada para evitar errores si das Enter sin querer
opcion_valida = false;
while ~opcion_valida
    opcion_validacion = input('-> Seleccione el escenario (1 o 2): ');
    if opcion_validacion == 1 || opcion_validacion == 2
        opcion_valida = true;
    else
        disp('Error: Por favor ingrese 1 o 2.');
    end
end

if opcion_validacion == 1
    %% CONFIGURACIÓN ESCENARIO A: GRAVEDAD INCLUIDA
    disp('>> Configurando simulación: MODELO COMPLETO (90°) ...');
    
    % 1. Matrices del Modelo LTI (Incluye término A21 de gravedad)
    A_sim = sys_90.A;
    B_sim = sys_90.B;
    
    % 2. Configuración de la Planta No Lineal
    k_g_interna = 1;           % HABILITA la realimentación interna de gravedad
    th_init     = val_th_m_90; % Condición Inicial: 90 grados
    
    % 3. Condiciones de Equilibrio (Bias)
    % Se inyecta el voltaje necesario para sostener la carga
    v_q_eq = R_sREF * val_i_q_90;
    U_bias = [v_q_eq; 0; 0];  
    
    % 4. Ajuste de Observación (Offset)
    % Resta el equilibrio para que la gráfica empiece en 0
    X_bias = [val_th_m_90; 0; val_i_q_90; 0; 0; 20]; 

else
    %% CONFIGURACIÓN ESCENARIO B: GRAVEDAD COMO PERTURBACIÓN (DESACOPLADA)
    disp('>> Configurando simulación: MODELO NOMINAL DE DISEÑO ...');
    disp('   (Se asume compensación perfecta de gravedad para validar J_eq)');
    
    % 1. Matrices del Modelo LTI (Sin término A21)
    % Usamos el modelo sys_design que creaste en la Sección C
    A_sim = sys_design.A;
    B_sim = sys_design.B;         
    
    % 2. Configuración de la Planta No Lineal
    k_g_interna = 0;          % DESACOPLA la gravedad interna (se trata como d(t)=0)
    th_init     = 0;          % Condición Inicial: 0 (Incremental)
    
    % 3. Condiciones de Equilibrio
    U_bias = [0; 0; 0];       % Sin esfuerzo de sostenimiento (Equilibrio natural)
    
    % 4. Ajuste de Observación
    X_bias = zeros(6,1);      % Sin offset
end

% Matrices de Salida y Paso Directo (Comunes para ambos casos)
C_sim = eye(6);
D_sim = zeros(6,3);

disp(' ');
disp('-> Entorno configurado exitosamente.');
disp('-> Variables listas: A_sim, B_sim, U_bias, X_bias, k_g_interna, th_init');