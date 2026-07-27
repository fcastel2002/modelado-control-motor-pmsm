%% graficar_comp_codcont.m
% Comparacion controlador CONTINUO vs CODIFICADO (discreto, Tustin).
%   Run 42 = CONTINUO ; Run 41 = CODIFICADO
% Sensores no ideales N=3, Tld=0, trayectoria curva S 0->2pi->0 (20 s).
% La diferencia entre ambos refleja SOLO el error de discretizacion.
% Paleta coherente: continuo azul (#0072BD), codificado naranja (#D95319),
% error gris (#7f7f7f). Ambas curvas solidas.

close all;
outdir = fullfile(fileparts(mfilename('fullpath')), '..', ...
                  'Trabajo_Final___AyME___2026', 'imagenes');
assert(exist(outdir,'dir')==7, 'No existe la carpeta de imagenes');

cCon = '#0072BD';  cCod = '#D95319';  cErr = '#7f7f7f';
FS = 13; LW = 1.7; DPI = 300;
OM = char(969); TH = char(952);
NM_POS = [TH 'l,med(t) [rad]'];
NM_VEL = [OM '~_l(t) [rad/s]'];

r_con = getrun(42);   % continuo
r_cod = getrun(41);   % codificado
tt = (0:2e-4:20)';    % grilla uniforme comun

% ================= FIG 1: posicion =================
yc = rs(r_con, NM_POS, tt);
yd = rs(r_cod, NM_POS, tt);
figure('Color','w','Position',[100 100 1000 680]);
ax1 = subplot(2,1,1);
plot(tt,yc,'Color',cCon,'LineWidth',1.8,'DisplayName','Continuo'); hold on;
plot(tt,yd,'Color',cCod,'LineStyle','--','LineWidth',1.6,'DisplayName','Codificado'); hold off;
title('Seguimiento de posición \theta_l: continuo vs. codificado','FontSize',15);
ylabel('\theta_l [rad]','FontSize',14);
legend('Location','south','Orientation','horizontal','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
h0 = yline(0,'LineStyle',':','Color',[.6 .6 .6],'LineWidth',1);
h0.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on;
plot(tt,yd-yc,'Color',cErr,'LineWidth',1.3); hold off;
title('Error de discretización (codificado − continuo)','FontSize',14);
ylabel('\Delta\theta_l [rad]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
fmtejes(FS); linkaxes([ax1 ax2],'x'); xlim([0 20]);
exportgraphics(gcf, fullfile(outdir,'comp_codcont_pos.png'), 'Resolution', DPI);
fprintf('  comp_codcont_pos.png\n');

% ================= FIG 2: velocidad =================
yc = rs(r_con, NM_VEL, tt);
yd = rs(r_cod, NM_VEL, tt);
figure('Color','w','Position',[100 100 1000 680]);
ax1 = subplot(2,1,1);
plot(tt,yc,'Color',cCon,'LineWidth',1.8,'DisplayName','Continuo'); hold on;
plot(tt,yd,'Color',cCod,'LineStyle','--','LineWidth',1.6,'DisplayName','Codificado'); hold off;
title('Seguimiento de velocidad \omega_l: continuo vs. codificado','FontSize',15);
ylabel('\omega_l [rad/s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
h0 = yline(0,'LineStyle',':','Color',[.6 .6 .6],'LineWidth',1);
h0.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on;
plot(tt,yd-yc,'Color',cErr,'LineWidth',1.3); hold off;
title('Error de discretización (codificado − continuo)','FontSize',14);
ylabel('\Delta\omega_l [rad/s]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
fmtejes(FS); linkaxes([ax1 ax2],'x'); xlim([0 20]);
exportgraphics(gcf, fullfile(outdir,'comp_codcont_vel.png'), 'Resolution', DPI);
fprintf('  comp_codcont_vel.png\n');
fprintf('OK: figuras comp continuo/codificado generadas\n');

%% ===================== helpers =====================
function yy = rs(r, name, tt)
  [t,y] = sig(r,name);
  [tu,ia] = unique(t);
  yy = interp1(tu, y(ia), tt, 'linear', 'extrap');
end
function r = getrun(N)
  ids = Simulink.sdi.getAllRunIDs; r = [];
  for k = 1:numel(ids)
    rr = Simulink.sdi.getRun(ids(k)); tok = regexp(rr.Name,'Run\s+(\d+)','tokens','once');
    if ~isempty(tok) && str2double(tok{1})==N, r = rr; return; end
  end
  error('No se encontro Run %d', N);
end
function [t,Y] = sig(r,name)
  idx = [];
  for i = 1:r.SignalCount
    if strcmp(r.getSignalByIndex(i).Name, name), idx = i; break; end
  end
  if isempty(idx), error('Senal "%s" no encontrada en %s', name, r.Name); end
  v = r.getSignalByIndex(idx).Values; t = v.Time(:); Y = squeeze(v.Data);
  if isrow(Y), Y = Y(:); end
  if size(Y,1) ~= numel(t) && size(Y,2) == numel(t), Y = Y.'; end
end
function fmtejes(FS)
  grid on; grid minor;
  set(gca,'FontSize',FS,'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.12);
end
