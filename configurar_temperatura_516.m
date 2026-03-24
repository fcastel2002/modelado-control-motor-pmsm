%% CONFIGURAR_TEMPERATURA_516
% Selecciona una temperatura de estator para la comparacion NL vs LTI de
% la seccion 5.1.6 y calcula la resistencia equivalente R_s_sim.
%
% Uso recomendado:
%   1) Ejecutar este script.
%   2) Elegir la temperatura T_s_sim deseada.
%   3) Correr la simulacion en Simulink usando R_s_sim en ambos modelos.
%   4) Repetir el script para otra temperatura.

clc;

if exist('ESCENARIO', 'var')
    escenario_activo = ESCENARIO;
else
    escenario_activo = 'NOMINAL';
    ESCENARIO = escenario_activo;
end

parametros_sistema_completo;

fprintf('=======================================================\n');
fprintf('   CONFIGURACION TERMICA PARA SIMULACION 5.1.6\n');
fprintf('=======================================================\n');
fprintf('Escenario activo: %s\n', escenario_activo);
fprintf('Referencia: T_sREF = %.1f degC, R_sREF = %.4f Ohm\n', T_sREF, R_sREF);
fprintf('Temperaturas sugeridas [degC]: 20, 40, 60, 80, 100, 115\n');
fprintf('\n');

T_default = T_sREF;
if exist('T_s_sim', 'var') && isnumeric(T_s_sim) && isscalar(T_s_sim) && isfinite(T_s_sim)
    T_default = T_s_sim;
    fprintf('Ultima temperatura detectada en workspace: %.2f degC\n', T_s_sim);
end

T_s_sim = pedir_temperatura(T_default);

R_s_sim = R_sREF * (1 + a_cu * (T_s_sim - T_sREF));
delta_R_s_sim = R_s_sim - R_sREF;
delta_R_s_sim_pct = 100 * delta_R_s_sim / R_sREF;
tau_q_sim = L_q / R_s_sim;
tau_d_sim = L_d / R_s_sim;

cfg_rs_sim = struct( ...
    'ESCENARIO', escenario_activo, ...
    'T_s_sim', T_s_sim, ...
    'R_s_sim', R_s_sim, ...
    'R_sREF', R_sREF, ...
    'a_cu', a_cu, ...
    'T_sREF', T_sREF, ...
    'delta_R_s_sim', delta_R_s_sim, ...
    'delta_R_s_sim_pct', delta_R_s_sim_pct, ...
    'tau_q_sim', tau_q_sim, ...
    'tau_d_sim', tau_d_sim);

fprintf('\n');
fprintf('--------------- RESULTADOS DE LA CONFIGURACION ---------\n');
fprintf('Temperatura elegida:      T_s_sim  = %8.3f degC\n', T_s_sim);
fprintf('Resistencia calculada:    R_s_sim  = %8.5f Ohm\n', R_s_sim);
fprintf('Variacion respecto ref.:  Delta Rs = %8.5f Ohm (%+7.3f %%)\n', ...
    delta_R_s_sim, delta_R_s_sim_pct);
fprintf('Constante de tiempo q:    tau_q    = %8.5f ms\n', tau_q_sim * 1e3);
fprintf('Constante de tiempo d:    tau_d    = %8.5f ms\n', tau_d_sim * 1e3);

if T_s_sim > 115
    fprintf('\nADVERTENCIA: T_s_sim supera el limite de 115 degC indicado en la consigna.\n');
elseif T_s_sim < -15
    fprintf('\nADVERTENCIA: T_s_sim queda por debajo del rango ambiental minimo de la consigna.\n');
end

fprintf('\nVariables listas en workspace para Simulink:\n');
fprintf('  - T_s_sim\n');
fprintf('  - R_s_sim\n');
fprintf('  - cfg_rs_sim\n');
fprintf('\nEjecute nuevamente este script cada vez que quiera cambiar la temperatura.\n');

function T_sel = pedir_temperatura(T_default)
prompt = sprintf('Ingrese T_s_sim [degC] (Enter = %.1f): ', T_default);

while true
    entrada = input(prompt, 's');

    if isempty(strtrim(entrada))
        T_sel = T_default;
        return;
    end

    T_sel = str2double(strrep(strtrim(entrada), ',', '.'));
    if isscalar(T_sel) && isfinite(T_sel)
        return;
    end

    fprintf('Entrada invalida. Ingrese un numero real en degC.\n');
end
end
