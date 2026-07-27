%% graficar_sensores_526.m
% Figuras del punto 5.2.6 (sensores no ideales) desde el Simulink Data Inspector.
% Variante de planta "perfil trapezoidal aceleracion1" (con sensores en espacio de estados).
%   Run 30 -> N=1 (polos lentos): el sistema oscila     -> sensores_no_ideales_1
%   Run 29 -> N=3 (polos rapidos): sigue la consigna     -> sensores_no_ideales_2
% Layout: arriba posicion medida + consigna + velocidad; abajo corriente i_qs.
% Misma paleta coherente que graficar_lazo_cerrado.m.

close all;
outdir = fullfile(fileparts(mfilename('fullpath')), '..', ...
                  'Trabajo_Final___AyME___2026', 'imagenes');
assert(exist(outdir,'dir')==7, 'No existe la carpeta de imagenes:\n  %s', outdir);

% paleta
cMed  = '#0072BD';   % medida -> azul
cCons = '#A2142F';   % consigna -> rojo vino, '--'
cVel  = '#7E2F8E';   % velocidad -> violeta
FS = 13; LWm = 1.7; LWc = 2.1;

% nombres de senal (ASCII-safe con char())
TH = char(952); OM = char(969);
NM_THMED = [TH 'l,med(t) [rad]'];
NM_CONS  = [TH 'l*(t) [rad]'];
NM_WL    = [OM '~_l(t) [rad/s]'];
NM_IQS   = 'I_qd0s [A](1)';

casos = { 30, 'sensores_no_ideales_1.png',  'Sensores no ideales, N = 1 (polos lentos): el sistema oscila'; ...
          29, 'sensores_no_ideales_2.png',  'Sensores no ideales, N = 3 (polos rápidos): sigue la consigna'; ...
          31, 'inversor_no_ideal_6000.png', 'Inversor no ideal, \omega_n = 6000 rad/s (lento): oscila'; ...
          29, 'inversor_no_ideal_12000.png','Inversor no ideal, \omega_n = 12000 rad/s: desempeño satisfactorio' };

for c = 1:size(casos,1)
  r = getrun(casos{c,1});
  [t,thm] = sig(r, NM_THMED);
  [tc,cn] = sig(r, NM_CONS);
  [tw,wl] = sig(r, NM_WL);
  [ti,iq] = sig(r, NM_IQS);
  figure('Color','w','Position',[100 100 1000 720]);
  ax1 = subplot(2,1,1);
  plot(t,thm,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l medida'); hold on;
  plot(tw,wl,'Color',cVel,'LineWidth',LWm,'DisplayName','\omega_l');
  plot(tc,cn,'Color',cCons,'LineStyle','--','LineWidth',LWc,'DisplayName','\theta_l^* (consigna)');
  hold off;
  title(casos{c,3},'FontSize',15);
  ylabel('\theta_l [rad]  /  \omega_l [rad/s]','FontSize',14);
  legend('Location','northeast','FontSize',11,'Box','on');
  fmtejes(FS); ax1.XTickLabel = []; ylim([-2.6 6.9]);
  ax2 = subplot(2,1,2);
  plot(ti,iq,'Color',cMed,'LineWidth',LWm);
  ylabel('i_{qs} [A]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
  fmtejes(FS); ylim([-6 6]);
  linkaxes([ax1 ax2],'x'); xlim([0 15]);
  exportgraphics(gcf, fullfile(outdir, casos{c,2}), 'Resolution', 300);
  fprintf('  generada: %s\n', casos{c,2});
end
fprintf('OK: figuras de sensores (5.2.6) generadas en\n  %s\n', outdir);

%% ===================== helpers =====================
function r = getrun(N)
  ids = Simulink.sdi.getAllRunIDs; r = [];
  for k = 1:numel(ids)
    rr = Simulink.sdi.getRun(ids(k));
    tok = regexp(rr.Name, 'Run\s+(\d+)', 'tokens', 'once');
    if ~isempty(tok) && str2double(tok{1}) == N, r = rr; return; end
  end
  error('No se encontro Run %d en el Data Inspector', N);
end

function [t, Y] = sig(r, name)
  idx = [];
  for i = 1:r.SignalCount
    if strcmp(r.getSignalByIndex(i).Name, name), idx = i; break; end
  end
  if isempty(idx), error('Senal "%s" no encontrada en %s', name, r.Name); end
  v = r.getSignalByIndex(idx).Values;
  t = v.Time(:); Y = squeeze(v.Data);
  if isrow(Y), Y = Y(:); end
  if size(Y,1) ~= numel(t) && size(Y,2) == numel(t), Y = Y.'; end
end

function fmtejes(FS)
  grid on; grid minor;
  set(gca,'FontSize',FS,'XMinorGrid','on','YMinorGrid','on', ...
          'GridAlpha',0.3,'MinorGridAlpha',0.12);
end
