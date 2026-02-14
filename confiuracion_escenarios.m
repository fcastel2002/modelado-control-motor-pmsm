%% CONFIGURACIÓN DE ESCENARIOS - PROYECTO GLOBAL INTEGRADOR
%  Este script permite seleccionar el caso de incertidumbre paramétrica
%  y carga automáticamente las variables necesarias para Simulink.
%  Autores: J. Calderón - F. Castel

clearvars; clc;

fprintf('=======================================================\n');
fprintf('   INICIALIZACIÓN DEL SISTEMA Y ESCENARIOS DE CARGA    \n');
fprintf('=======================================================\n');
fprintf('Seleccione el escenario de operación e incertidumbre:\n\n');
fprintf('   [0] NOMINAL    (Valores medios / Hoja de datos)\n');
fprintf('   [1] PEOR CASO  (Carga Máxima, Max Inercia, Max Resistencia)\n');
fprintf('   [2] MEJOR CASO (Sin Carga, Min Inercia, Min Resistencia)\n');
fprintf('\n');

% Validación de entrada
opcion_valida = false;
while ~opcion_valida
    seleccion = input('-> Ingrese opción (0/1/2): ');
    if ismember(seleccion, [0, 1, 2])
        opcion_valida = true;
    else
        fprintf('   Error: Opción no válida. Intente nuevamente.\n');
    end
end

% Asignación de variable de control
switch seleccion
    case 1
        ESCENARIO = 'WORST';  % Peor Caso
        nombre_esc = 'PEOR CASO (Carga Máxima)';
    case 2
        ESCENARIO = 'BEST';   % Mejor Caso
        nombre_esc = 'MEJOR CASO (Liviano)';
    otherwise
        ESCENARIO = 'NOMINAL';
        nombre_esc = 'NOMINAL (Estándar)';
end

fprintf('\n-------------------------------------------------------\n');
fprintf('-> Configurando parámetros para: [%s]\n', nombre_esc);

% --- EJECUCIÓN DEL ARCHIVO DE PARÁMETROS ---
% Al correr esto, el script 'parametros_sistema_completo' detectará
% la variable 'ESCENARIO' y actuará en consecuencia.
run('parametros_sistema_completo.mlx'); 

% --- REPORTE DE VALORES CARGADOS ---
fprintf('-> Parámetros físicos actualizados exitosamente.\n');
fprintf('\n   RESUMEN DE VARIABLES CLAVE:\n');
fprintf('   > Masa Carga Útil (m_l):  %6.3f kg\n', m_l);
fprintf('   > Inercia Total (J_eq):   %6.3e kg.m2\n', J_eq);
fprintf('   > Resistencia (Rs_ref):   %6.4f Ohm\n', R_sREF);
fprintf('   > Fricción Total (b_eq):  %6.3e N.m/(rad/s)\n', b_eq);

fprintf('-------------------------------------------------------\n');
fprintf('¡Listo! Puede ejecutar la simulación en Simulink ahora.\n');