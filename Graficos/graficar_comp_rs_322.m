%% graficar_comp_rs_322.m
% Figuras del punto 3.2.2: efecto de estimar R_s(T_s) en el desacople del
% modulador de torque vs. usar un valor nominal constante.
% Ensayo open-loop (sin cascada de posicion), StopTime 4 s, consigna de
% torque constante T_m* = 0,1 N.m.  La planta NL siempre tiene R_s(T_s) real.
%   Run 34 -> R_s VARIABLE (estimada en linea): T_neto se mantiene en ~0,1
%   Run 33 -> R_s CONSTANTE (nominal 25 C): T_neto cae ~15% al calentarse
% Paleta coherente (misma que graficar_lazo_cerrado.m / graficar_sensores_526.m):
%   R_s variable = azul #0072BD (caso correcto/base)
%   R_s constante = naranja #D95319 (2da condicion)
% Ambas curvas SOLIDAS, sin punteado (las rectas se separan y se ven bien).

close all;
outdir = fullfile(fileparts(mfilename('fullpath')), '..', ...
                  'Trabajo_Final___AyME___2026', 'imagenes');
assert(exist(outdir,'dir')==7, 'No existe la carpeta de imagenes:\n  %s', outdir);

% paleta
cVar = '#0072BD';   % R_s variable (estimada) -> azul
cCte = '#D95319';   % R_s constante (nominal) -> naranja
FS = 13; LW = 1.9; DPI = 300;

% nombres de senal (ASCII-safe con char())
OM = char(969);
NM_TNET = 'T_neto [N.m]';
NM_WM   = [OM 'm(t) [rad/s]'];

r_var = getrun(34);    % R_s variable
r_cte = getrun(33);    % R_s constante

%% ---------- FIG 1: comp_torque_rs ----------
[tv,tnv] = sig(r_var, NM_TNET);
[tc,tnc] = sig(r_cte, NM_TNET);
figure('Color','w','Position',[100 100 1000 520]);
plot(tv,tnv,'Color',cVar,'LineWidth',LW,'DisplayName','R_s variable (estimada)'); hold on;
plot(tc,tnc,'Color',cCte,'LineWidth',LW,'DisplayName','R_s constante (nominal)');
hold off;
title('Torque neto entregado: R_s estimada vs. R_s constante','FontSize',15);
ylabel('T_{neto} [N\cdotm]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','southwest','FontSize',12,'Box','on');
fmtejes(FS); xlim([0 4]); ylim([0 0.105]);
exportgraphics(gcf, fullfile(outdir,'comp_torque_rs.png'), 'Resolution', DPI);
fprintf('  generada: comp_torque_rs.png\n');

%% ---------- FIG 2: comp_velocidad_rs ----------
[wv,wmv] = sig(r_var, NM_WM);
[wc,wmc] = sig(r_cte, NM_WM);
figure('Color','w','Position',[100 100 1000 520]);
plot(wv,wmv,'Color',cVar,'LineWidth',LW,'DisplayName','R_s variable (estimada)'); hold on;
plot(wc,wmc,'Color',cCte,'LineWidth',LW,'DisplayName','R_s constante (nominal)');
hold off;
title('Velocidad angular \omega_m: R_s estimada vs. R_s constante','FontSize',15);
ylabel('\omega_m [rad/s]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','northwest','FontSize',12,'Box','on');
fmtejes(FS); xlim([0 4]);
exportgraphics(gcf, fullfile(outdir,'comp_velocidad_rs.png'), 'Resolution', DPI);
fprintf('  generada: comp_velocidad_rs.png\n');

fprintf('OK: figuras 3.2.2 (comparacion R_s) generadas en\n  %s\n', outdir);

%% ===================== helpers =====================
function r = getrun(N)
  ids = Simulink.sdi.getAllRunIDs; r = [];
  for k = 1:numel(ids)
    rr = Simulink.sdi.getRun(ids(k));
    tok = regexp(rr.Name,'Run\s+(\d+)','tokens','once');
    if ~isempty(tok) && str2double(tok{1})==N, r = rr; return; end
  end
  error('No se encontro Run %d en el Data Inspector', N);
end

function [t,Y] = sig(r,name)
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
