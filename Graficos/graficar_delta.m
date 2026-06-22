% graficar_delta.m
% Grafica el angulo de carga (arriba) y las consignas v_qs* y T_ld (abajo)
% de la ULTIMA corrida del dataset 'datos' exportado del Data Inspector.
%
% Acceso por INDICE (los nombres tienen caracteres especiales y duplicados):
%   14: NL delta [deg]   (angulo de carga, ya en grados)
%   19: V_qs*
%   18: T_l

run = datos{end};        % ultima corrida

delta = run{14}.Values;
vqs   = run{19}.Values;
tld   = run{18}.Values;

% Extraer tiempo y datos como vectores columna.
% squeeze() colapsa dimensiones sobrantes (algunas senales vienen en 3D).
td = delta.Time(:);   yd = squeeze(delta.Data); yd = yd(:);
tv = vqs.Time(:);     yv = squeeze(vqs.Data);   yv = yv(:);
tt = tld.Time(:);     yt = squeeze(tld.Data);   yt = yt(:);

% --- Graficar ---
figure('Color', 'w');

ax1 = subplot(2,1,1);
plot(td, yd, 'Color', '#ff13a6', 'LineWidth', 1.6);
grid on;
ylabel('\delta [\circ]', 'FontSize', 14);
title('Ángulo de carga', 'FontSize', 15);
set(gca, 'FontSize', 13);

ax2 = subplot(2,1,2);
plot(tv, yv, 'Color', '#b746ff', 'LineWidth', 1.6); hold on;
plot(tt, yt, 'Color', '#64d413', 'LineWidth', 1.6); hold off;
grid on;
xlabel('Tiempo [s]', 'FontSize', 14);
ylabel('Consignas', 'FontSize', 14);
legend('v_{qs}^* [V]', 'T_{ld} [N\cdotm]', 'FontSize', 12, 'Location', 'best');
set(gca, 'FontSize', 13);

linkaxes([ax1 ax2], 'x');
xlim([0 2]);
