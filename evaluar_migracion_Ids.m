%% SCRIPT: EVALUACIÓN DE MIGRACIÓN DE POLOS AL VARIAR I_ds (LPV)
% Objetivo: Evaluar la migración de autovalores de la matriz AO al variar I_ds.
% Este script es independiente y se puede usar en el enfoque LPV.


%% 1. Configuración de Variación de I_ds
disp('=== Configuración de Variación de I_ds ===');
I_ds_min = input('-> Ingrese el valor mínimo de I_ds [A] (ej. -10) o Enter para -10: ');
if isempty(I_ds_min), I_ds_min = -10; end

I_ds_max = input('-> Ingrese el valor máximo de I_ds [A] (ej. 10) o Enter para 10: ');
if isempty(I_ds_max), I_ds_max = 10; end

I_ds_step = input('-> Ingrese el paso de corriente [A] (ej. 0.5) o Enter para 0.5: ');
if isempty(I_ds_step), I_ds_step = 0.5; end

Ids_vals = I_ds_min:I_ds_step:I_ds_max;
num_points = length(Ids_vals);

%% 2. Definición Arbitraría de Parámetros (Ajustar según se necesite)
% Como estamos en un contexto general LPV, definimos parámetros que pueden
% ser modificados a gusto por el usuario sin depender de otros scripts.

params = struct();
% --- Subsistema electromagnético ---
params.Ld        = L_d;
params.Lq        = L_q;
params.Lls       = L_ls;
params.Rs_ref    = R_sREF;
params.alpha_cu  = a_cu;
params.lambda_m  = lambda_r;
params.Pp        = P_p;

% --- Subsistema mecánico ---
params.J_total   = J_eq;
params.b_total   = b_eq;
params.r_red     = r;
params.kl        = k_l;
params.g         = g;

% --- Subsistema térmico ---
params.Cts       = C_ts;
params.Rts_amb   = R_ts_amb;
params.Ts_ref    = T_sREF;      % Temperatura de referencia

%% 3. Definición del Punto de Operación Base
% Vector x_op = [theta_mO; w_mO; i_qsO; i_dsO; i_0sO; TsO]
theta_mO = 0;              
w_mO = 6.28;                % Una velocidad mecánica arbitraria
i_qsO = 5;                 % Corriente en q arbitraria
i_0sO = 0;                 % Secuencia cero 
TsO = 20;                  % Temperatura de operación

%% 4. Evaluación de Autovalores en Bucle
disp('-> Calculando autovalores de la matriz AO_LPV...');
polos_totales = zeros(6, num_points); % Almacena 6 polos por iteración

for k = 1:num_points
    i_dsO_actual = Ids_vals(k);
    
    % Actualizamos el x_op particular para esta evaluación
    x_op = [theta_mO; w_mO; i_qsO; i_dsO_actual; i_0sO; TsO];
    
    % Obtener la matriz AO para este conjunto LPV
    AO_actual = calcular_matriz_AO(x_op, params);
    
    % Calcular los autovalores y organizarlos para seguir ramas de forma más estable
    autovalores = sort(eig(AO_actual), 'ComparisonMethod', 'real');
    
    polos_totales(:, k) = autovalores;
end

%% 5. Gráficas y Resultados
disp('-> Generando diagrama de polos...');

% Colores y marcadores distintos para cada rama (autovalor)
colores_rama = lines(6);
marcadores   = {'o', 's', 'd', '^', 'v', 'p'};
nombres_polo = {'\lambda_1', '\lambda_2', '\lambda_3', ...
                '\lambda_4', '\lambda_5', '\lambda_6'};

figure('Name', 'Migración de Polos vs I_{ds} (LPV)', 'Color', 'w', 'Position', [100, 100, 900, 650]);
hold on; grid on; box on;

title('Migración de Polos al Variar la Corriente I_{ds}', 'FontSize', 14);
xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');

h_leg = gobjects(6, 1); % Handles para la leyenda

for polo_idx = 1:6
    % Traza continua de la rama
    h_leg(polo_idx) = plot(real(polos_totales(polo_idx, :)), ...
        imag(polos_totales(polo_idx, :)), ...
        '-', 'Color', colores_rama(polo_idx,:), 'LineWidth', 1.2, ...
        'Marker', marcadores{polo_idx}, 'MarkerSize', 5, ...
        'DisplayName', nombres_polo{polo_idx});

    % Marcador de INICIO (I_ds mínimo) -> círculo relleno grande
    plot(real(polos_totales(polo_idx, 1)), imag(polos_totales(polo_idx, 1)), ...
        marcadores{polo_idx}, 'Color', colores_rama(polo_idx,:), ...
        'MarkerFaceColor', colores_rama(polo_idx,:), 'MarkerSize', 10, ...
        'HandleVisibility', 'off');

    % Marcador de FIN (I_ds máximo) -> marcador con borde negro
    plot(real(polos_totales(polo_idx, end)), imag(polos_totales(polo_idx, end)), ...
        marcadores{polo_idx}, 'Color', 'k', ...
        'MarkerFaceColor', colores_rama(polo_idx,:), 'MarkerSize', 10, ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% Marcador especial para I_ds = 0
idx_cero = find(abs(Ids_vals) < eps);
if ~isempty(idx_cero)
    h_ids0 = plot(real(polos_totales(:, idx_cero)), imag(polos_totales(:, idx_cero)), ...
        'r*', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', 'I_{ds} = 0');
else
    warning('El valor I_{ds} = 0 no está en el rango definido.');
    h_ids0 = [];
end

% Líneas de referencia (ejes 0)
xline(0, 'k--', 'LineWidth', 1);
yline(0, 'k--', 'LineWidth', 1);

% Leyenda con las ramas + marcador I_ds=0
if ~isempty(h_ids0)
    legend([h_leg; h_ids0(1)], [nombres_polo, {'I_{ds} = 0'}], 'Location', 'best', 'FontSize', 11);
else
    legend(h_leg, nombres_polo, 'Location', 'best', 'FontSize', 11);
end

% Subtítulo indicando rango y dirección de migración
subtitle(sprintf('Relleno = I_{ds,min} (%.2g A)  |  Borde negro = I_{ds,max} (%.2g A)', ...
    I_ds_min, I_ds_max), 'FontSize', 10);

%% 6. Salida en formato CSV para LaTeX
disp(' ');
disp('=== Autovalores en formato CSV (para tabla LaTeX) ===');
disp('I_ds [A], Re(lambda_1), Im(lambda_1), Re(lambda_2), Im(lambda_2), Re(lambda_3), Im(lambda_3), Re(lambda_4), Im(lambda_4), Re(lambda_5), Im(lambda_5), Re(lambda_6), Im(lambda_6)');
for k = 1:num_points
    fprintf('%.4f', Ids_vals(k));
    for polo_idx = 1:6
        fprintf(', %.6f, %.6f', real(polos_totales(polo_idx, k)), imag(polos_totales(polo_idx, k)));
    end
    fprintf('\n');
end

disp(' ');
disp('=== Análisis Terminado ===');
