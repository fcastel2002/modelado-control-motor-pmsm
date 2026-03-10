%% SIMULACION_5_1_REVISADA
% Base corregida para rehacer 5.1.6 sin mezclar pares, angulos ni unidades.
% No pisa ningun archivo original.

clearvars; clc;

p = parametros_revisados_5_1();
mdl = modelo_lti_5_1_revisado(p);

fprintf('Modelo nominal cargado.\n');
fprintf('omega_m,nom = %.3f rad/s\n', p.omega_m_nom);
fprintf('tau_q = %.6f s, tau_d = %.6f s, tau_m = %.6f s\n', ...
    p.tau_q, p.tau_d, p.tau_m);

fprintf('\nSeñales exigidas por consigna 5.1.6:\n');
fprintf('v_qs*: +19.596 V [0.1,0.7), -19.596 V [1.1,1.7)\n');
fprintf('T_ld:  +6.28 Nm [0.3,0.5) y [1.3,1.5), -6.28 Nm [0.5,0.9) y [1.5,1.9)\n');

fprintf('\nChequeos previos obligatorios antes de correr en Simulink/MATLAB:\n');
fprintf('1) Mantener theta_m mecanico y theta_r = P_p*theta_m separados.\n');
fprintf('2) Usar T_l para torque total de carga y T_ld solo para disturbio adicional.\n');
fprintf('3) No inferir delta(t) a partir de atan2(i_d, i_q); definir theta_ev fisico.\n');
fprintf('4) Verificar que omega_m no exceda %.1f rad/s sin justificar saturacion.\n', p.omega_m_nom);

fprintf('\nMatrices del modelo aumentado:\n');
disp(mdl.A_aug);
