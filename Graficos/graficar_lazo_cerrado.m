%% graficar_lazo_cerrado.m
% Figuras del lazo cerrado (puntos 5.2.4 y 5.2.5) en estilo grande y con una
% PALETA DE COLORES COHERENTE en todas las figuras.
%
% Toma las corridas directamente del Simulink Data Inspector, por numero de Run:
%   Run 10 -> seguimiento trapezoidal SIN saturacion y SIN perturbacion (5.2.4)
%   Run 11 -> seguimiento CON perturbacion de carga Tld escalonada       (5.2.4)
%   Run 12 -> seguimiento CON saturacion de tension del inversor          (5.2.5)
%
% Guarda los PNG (300 dpi) en Trabajo_Final___AyME___2026/imagenes con los
% nombres que usa el informe.  Ejecutar con:  run('graficar_lazo_cerrado.m')

close all;

outdir = fullfile(fileparts(mfilename('fullpath')), '..', ...
                  'Trabajo_Final___AyME___2026', 'imagenes');
assert(exist(outdir,'dir')==7, 'No existe la carpeta de imagenes:\n  %s', outdir);

% ===================== PALETA COHERENTE =====================
cCons = '#A2142F';        % CONSIGNA / referencia   -> rojo vino ('--', siempre encima)
cMed  = '#0072BD';        % MEDIDA (respuesta base)  -> azul
cErr  = '#7f7f7f';        % ERROR                    -> gris
cPert = '#64d413';        % PERTURBACION T_ld        -> verde claro
cSat  = '#D95319';        % 2da condicion (CON sat.) -> naranja
cVel  = '#7E2F8E';        % velocidad (perfiles de trayectoria) -> violeta
cLim  = [0.24 0.24 0.24]; % limites / especificaciones -> gris oscuro
FS = 13; LWm = 1.7; LWc = 2.1; DPI = 300;
sv = @(n) exportgraphics(gcf, fullfile(outdir,n), 'Resolution', DPI);

% ===================== NOMBRES DE SENAL (ASCII-safe) =====================
TH = char(952); OM = char(969);        % theta, omega   (~ = tilde ASCII 126)
NM_THL  = [TH 'l(t) [rad]'];
NM_CONS = [TH 'l*(t) [rad]'];
NM_WL   = [OM '~_l(t) [rad/s]'];
NM_IRMS = 'I_rms';
NM_TNET = 'T_neto [N.m]';
NM_TLD  = 'Tld [N.m]';
NM_IMAX = 'I_rms_max';
NM_INOM = 'I_rms_nom';
NM_THM  = [TH 'm(t) [rad]'];
NM_WEM  = [OM '~_m(t) [rad/s]'];

% ===================== CORRIDAS =====================
r_ref  = getrun(10);       % sin saturacion, sin perturbacion
r_pert = getrun(11);       % con perturbacion
r_sat  = getrun(12);       % con saturacion de tension
r_cur  = getrun(13);       % con saturacion de la consigna de corriente
r_vel  = getrun(15);       % perfil trapezoidal de VELOCIDAD (mejora de trayectoria)
r_acc  = getrun(17);       % perfil trapezoidal de ACELERACION (curva S)
r_obs  = getrun(22);       % observador integral + perturbacion Tld (5.2.5.B)
r_term = getrun(27);       % termico: Tld=5 const, arranca en 25C=T_amb, consigna repetida (5.2.5.C)
Tmax = 15;

%% ---------- FIG 1: pos-velocidad-simulacion_cont (run 10) ----------
[t,thl]  = sig(r_ref, NM_THL);
[~,cons] = sig(r_ref, NM_CONS);
[tw,wl]  = sig(r_ref, NM_WL);
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t,thl,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l (medida)'); hold on;
plot(t,cons,'Color',cCons,'LineStyle','--','LineWidth',LWc,'DisplayName','\theta_l^* (consigna)');
hold off; title('Seguimiento de la consigna trapezoidal','FontSize',15);
ylabel('\theta_l [rad]','FontSize',14);
legend('Location','south','Orientation','horizontal','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
plot(tw,wl,'Color',cMed,'LineWidth',LWm);
ylabel('\omega_l [rad/s]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
fmtejes(FS); linkaxes([ax1 ax2],'x'); xlim([0 Tmax]);
sv('pos-velocidad-simulacion_cont.png');

%% ---------- FIG 2: overshoot_transitorio_trap_pos (run 10, zoom t=5) ----------
[t,thl]  = sig(r_ref, NM_THL);
[~,cons] = sig(r_ref, NM_CONS);
xr = [4.98 5.04]; err = thl - cons;
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t,thl,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l (medida)'); hold on;
plot(t,cons,'Color',cCons,'LineStyle','--','LineWidth',LWc,'DisplayName','\theta_l^* (consigna)');
hold off; title('Overshoot del transitorio rampa-meseta (t = 5 s)','FontSize',15);
ylabel('\theta_l [rad]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = []; setYzoom(ax1,t,[thl cons],xr);
ax2 = subplot(2,1,2);
h0 = yline(0,'LineStyle',':','Color',[.6 .6 .6],'LineWidth',1);
h0.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on;
plot(t,err,'Color',cErr,'LineWidth',LWm); hold off;
ylabel('\theta_l - \theta_l^* [rad]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
fmtejes(FS); setYzoom(ax2,t,err,xr);
linkaxes([ax1 ax2],'x'); xlim(xr);
sv('overshoot_transitorio_trap_pos.png');

%% ---------- FIG 3: torque_pos_relacion (run 10) ----------
[t,tn]   = sig(r_ref, NM_TNET);
[tp,thl] = sig(r_ref, NM_THL);
[~,cons] = sig(r_ref, NM_CONS);
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t,tn,'Color',cMed,'LineWidth',LWm);
title('Torque neto y posición de la articulación','FontSize',15);
ylabel('T_{neto} [N\cdotm]','FontSize',14);
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
plot(tp,thl,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l (medida)'); hold on;
plot(tp,cons,'Color',cCons,'LineStyle','--','LineWidth',LWc,'DisplayName','\theta_l^* (consigna)');
hold off; ylabel('\theta_l [rad]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','south','Orientation','horizontal','FontSize',11,'Box','on');
fmtejes(FS); linkaxes([ax1 ax2],'x'); xlim([0 Tmax]);
sv('torque_pos_relacion.png');

%% ---------- FIG 4: picos_corriente_trapz_pos (run 10) ----------
[t,ir]  = sig(r_ref, NM_IRMS);
[~,imx] = sig(r_ref, NM_IMAX);  imx = imx(1);
[~,inm] = sig(r_ref, NM_INOM);  inm = inm(1);
figure('Color','w','Position',[100 100 1000 500]);
plot(t,ir,'Color',cMed,'LineWidth',LWm,'DisplayName','I_{rms}'); hold on;
plot([0 Tmax],[imx imx],'--','Color',cLim,'LineWidth',1.8,'DisplayName','I_{max} = 2 A');
plot([0 Tmax],[inm inm],':','Color',[.5 .5 .5],'LineWidth',1.6,'DisplayName','I_{nom} = 0,4 A');
hold off; title('Picos de corriente rms en los transitorios','FontSize',15);
ylabel('I_{rms} [A]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); xlim([0 Tmax]);
sv('picos_corriente_trapz_pos.png');

%% ---------- FIG 5: plot_perturb_error (run 11) ----------
[t,thl]  = sig(r_pert, NM_THL);
[~,cons] = sig(r_pert, NM_CONS);
err = thl - cons;
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t,thl,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l (medida)'); hold on;
plot(t,cons,'Color',cCons,'LineStyle','--','LineWidth',LWc,'DisplayName','\theta_l^* (consigna)');
hold off; title('Seguimiento de posición bajo perturbación de carga','FontSize',15);
ylabel('\theta_l [rad]','FontSize',14);
legend('Location','south','Orientation','horizontal','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
h0 = yline(0,'LineStyle',':','Color',[.6 .6 .6],'LineWidth',1);
h0.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on;
plot(t,err,'Color',cErr,'LineWidth',LWm); hold off;
title('Error de seguimiento (crece con la perturbación)','FontSize',15);
ylabel('\theta_l - \theta_l^* [rad]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
fmtejes(FS); linkaxes([ax1 ax2],'x'); xlim([0 Tmax]);
sv('plot_perturb_error.png');

%% ---------- FIG 6: perturb_y_posicion_plot (run 11) ----------
[t,thl]  = sig(r_pert, NM_THL);
[~,cons] = sig(r_pert, NM_CONS);
[tt,tld] = sig(r_pert, NM_TLD);
figure('Color','w','Position',[100 100 1000 620]);
plot(t,thl,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l (medida)'); hold on;
plot(tt,tld,'Color',cPert,'LineWidth',LWm,'DisplayName','T_{ld} (perturbación)');
plot(t,cons,'Color',cCons,'LineStyle','--','LineWidth',LWc,'DisplayName','\theta_l^* (consigna)');
hold off; title('Perturbación, consigna y posición medida','FontSize',15);
ylabel('\theta_l [rad]  /  T_{ld} [N\cdotm]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','southeast','FontSize',12,'Box','on');
fmtejes(FS); xlim([0 Tmax]);
sv('perturb_y_posicion_plot.png');

%% ---------- FIG 7: corrientes_limitando_voltaje (run 10 sin vs run 12 con) ----------
[t0,ir0] = sig(r_ref, NM_IRMS);   % sin saturacion
[t1,ir1] = sig(r_sat, NM_IRMS);   % con saturacion de tension
[~,imx]  = sig(r_ref, NM_IMAX);  imx = imx(1);
[~,inm]  = sig(r_ref, NM_INOM);  inm = inm(1);
figure('Color','w','Position',[100 100 1000 560]);
plot(t0,ir0,'Color',cMed,'LineWidth',LWm,'DisplayName','I_{rms} sin saturación'); hold on;
plot(t1,ir1,'Color',cSat,'LineWidth',LWm,'DisplayName','I_{rms} con saturación');
plot([0 Tmax],[imx imx],'--','Color',cLim,'LineWidth',1.8,'DisplayName','I_{max} = 2 A');
plot([0 Tmax],[inm inm],':','Color',[.5 .5 .5],'LineWidth',1.6,'DisplayName','I_{nom} = 0,4 A');
hold off; title('Corriente rms del estator con y sin saturación de tensión','FontSize',15);
ylabel('I_{rms} [A]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); xlim([0 Tmax]);
sv('corrientes_limitando_voltaje.png');

%% ---------- FIG 8: pos_voltaje_comparacion (run 10 sin vs run 12 con) ----------
[t0,th0] = sig(r_ref, NM_THL);  [~,c0] = sig(r_ref, NM_CONS);
[t1,th1] = sig(r_sat, NM_THL);  [~,c1] = sig(r_sat, NM_CONS);
e0 = th0 - c0;  e1 = th1 - c1;
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t0,th0,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l sin saturación'); hold on;
plot(t1,th1,'Color',cSat,'LineWidth',LWm,'DisplayName','\theta_l con saturación');
plot(t0,c0,'Color',cCons,'LineStyle','--','LineWidth',LWc,'DisplayName','\theta_l^* (consigna)');
hold off; title('Seguimiento de posición con y sin saturación de tensión','FontSize',15);
ylabel('\theta_l [rad]','FontSize',14);
legend('Location','south','Orientation','horizontal','FontSize',11,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
h0 = yline(0,'LineStyle',':','Color',[.6 .6 .6],'LineWidth',1);
h0.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on;
plot(t0,e0,'Color',cMed,'LineWidth',LWm,'DisplayName','sin saturación');
plot(t1,e1,'Color',cSat,'LineWidth',LWm,'DisplayName','con saturación');
hold off; title('Error de seguimiento \theta_l - \theta_l^*','FontSize',15);
ylabel('\theta_l - \theta_l^* [rad]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','northwest','FontSize',11,'Box','on');
fmtejes(FS); linkaxes([ax1 ax2],'x'); xlim([0 Tmax]);
sv('pos_voltaje_comparacion.png');

%% ---------- FIG 9: corrientes_comparacion (run 13, con sat de consigna de corriente) ----------
[t,ir]  = sig(r_cur, NM_IRMS);
[~,imx] = sig(r_cur, NM_IMAX);  imx = imx(1);
[~,inm] = sig(r_cur, NM_INOM);  inm = inm(1);
figure('Color','w','Position',[100 100 1000 520]);
plot(t,ir,'Color',cSat,'LineWidth',LWm,'DisplayName','I_{rms} (consigna saturada)'); hold on;
plot([0 Tmax],[imx imx],'--','Color',cLim,'LineWidth',1.8,'DisplayName','I_{max} = 2 A');
plot([0 Tmax],[inm inm],':','Color',[.5 .5 .5],'LineWidth',1.6,'DisplayName','I_{nom} = 0,4 A');
hold off; title('Corriente rms con saturación de la consigna de corriente','FontSize',15);
ylabel('I_{rms} [A]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','northwest','FontSize',12,'Box','on');
fmtejes(FS); xlim([0 Tmax]);
sv('corrientes_comparacion.png');

%% ---------- FIG 10: perfil_trap_vel_comp (run 15, perfil trapezoidal de velocidad) ----------
[t,thl] = sig(r_vel, NM_THL);
[tw,wl] = sig(r_vel, NM_WL);
[ti,ir] = sig(r_vel, NM_IRMS);
[~,inm] = sig(r_vel, NM_INOM);  inm = inm(1);
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t,thl,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l [rad]'); hold on;
plot(tw,wl,'Color',cVel,'LineWidth',LWm,'DisplayName','\omega_l [rad/s]'); hold off;
title('Perfil trapezoidal de velocidad (t_{rampa} = 0,1 s)','FontSize',15);
ylabel('\theta_l [rad]  /  \omega_l [rad/s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
plot(ti,ir,'Color',cMed,'LineWidth',LWm,'DisplayName','I_{rms}'); hold on;
plot([0 Tmax],[inm inm],':','Color',[.5 .5 .5],'LineWidth',1.6,'DisplayName','I_{nom} = 0,4 A'); hold off;
ylabel('I_{rms} [A]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); ylim([0 0.6]); linkaxes([ax1 ax2],'x'); xlim([0 Tmax]);
sv('perfil_trap_vel_comp.png');

%% ---------- FIG 11: perfil_trap_accel_comp (run 17, perfil trapezoidal de aceleracion) ----------
[t,thl] = sig(r_acc, NM_THL);
[tw,wl] = sig(r_acc, NM_WL);
[ti,ir] = sig(r_acc, NM_IRMS);
[~,inm] = sig(r_acc, NM_INOM);  inm = inm(1);
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t,thl,'Color',cMed,'LineWidth',LWm,'DisplayName','\theta_l [rad]'); hold on;
plot(tw,wl,'Color',cVel,'LineWidth',LWm,'DisplayName','\omega_l [rad/s]'); hold off;
title('Perfil trapezoidal de aceleración (curva S)','FontSize',15);
ylabel('\theta_l [rad]  /  \omega_l [rad/s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
plot(ti,ir,'Color',cMed,'LineWidth',LWm,'DisplayName','I_{rms}'); hold on;
plot([0 Tmax],[inm inm],':','Color',[.5 .5 .5],'LineWidth',1.6,'DisplayName','I_{nom} = 0,4 A'); hold off;
ylabel('I_{rms} [A]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); ylim([0 0.6]); linkaxes([ax1 ax2],'x'); xlim([0 Tmax]);
sv('perfil_trap_accel_comp.png');

%% ---------- FIG 12: 5_2_5_B_error_observador (run 11, observador CLASICO) ----------
[t,thm] = sig(r_pert, NM_THM);
[~,wem] = sig(r_pert, NM_WEM);
wtrue = gradient(thm)./gradient(t);
eobs = wtrue - wem;
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t,wtrue,'Color',cMed,'LineWidth',LWm,'DisplayName','\omega_m verdadera'); hold on;
plot(t,wem,'Color',cSat,'LineStyle','--','LineWidth',LWc,'DisplayName','\omega_m estimada'); hold off;
title('Observador clásico (2° orden): \omega_m verdadera vs. estimada','FontSize',15);
ylabel('\omega_m [rad/s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
h0 = yline(0,'LineStyle',':','Color',[.6 .6 .6],'LineWidth',1);
h0.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on;
plot(t,eobs,'Color',cErr,'LineWidth',LWm); hold off;
title('Error de estimación: offset de régimen permanente \propto T_{ld}','FontSize',15);
ylabel('\omega_m - \omega_m^{est} [rad/s]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
fmtejes(FS); ylim([-1.6 1.6]); linkaxes([ax1 ax2],'x'); xlim([0 Tmax]);
sv('5_2_5_B_error_observador.png');

%% ---------- FIG 13: 5_2_5_B_observador_sin_offset (run 22, observador INTEGRAL) ----------
[t,thm] = sig(r_obs, NM_THM);
[~,wem] = sig(r_obs, NM_WEM);
wtrue = gradient(thm)./gradient(t);
eobs = wtrue - wem;
figure('Color','w','Position',[100 100 1000 720]);
ax1 = subplot(2,1,1);
plot(t,wtrue,'Color',cMed,'LineWidth',LWm,'DisplayName','\omega_m verdadera'); hold on;
plot(t,wem,'Color',cSat,'LineStyle','--','LineWidth',LWc,'DisplayName','\omega_m estimada'); hold off;
title('Observador con acción integral: \omega_m verdadera vs. estimada','FontSize',15);
ylabel('\omega_m [rad/s]','FontSize',14);
legend('Location','northeast','FontSize',12,'Box','on');
fmtejes(FS); ax1.XTickLabel = [];
ax2 = subplot(2,1,2);
h0 = yline(0,'LineStyle',':','Color',[.6 .6 .6],'LineWidth',1);
h0.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on;
plot(t,eobs,'Color',cErr,'LineWidth',LWm); hold off;
title('Error de estimación: offset eliminado en régimen permanente','FontSize',15);
ylabel('\omega_m - \omega_m^{est} [rad/s]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
fmtejes(FS); ylim([-1.6 1.6]); linkaxes([ax1 ax2],'x'); xlim([0 Tmax]);
sv('5_2_5_B_observador_sin_offset.png');

%% ---------- FIG 14: 5_2_5_C_comportamiento_termico (run 25, Tld=5 constante) ----------
zTs = '';
for zi = 1:r_term.SignalCount
  znm = r_term.getSignalByIndex(zi).Name;
  if contains(znm,'T_s') && contains(znm,'C]'), zTs = znm; break; end
end
[t,ts] = sig(r_term, zTs);
Tsmax = 115; ic = find(ts>=Tsmax,1);
if ~isempty(ic) && ic>1, Tcross = interp1(ts(ic-1:ic), t(ic-1:ic), Tsmax); else, Tcross = t(ic); end
fprintf('  [termico] cruce de 115 C: t = %.2f s\n', Tcross);
figure('Color','w','Position',[100 100 1000 560]);
plot(t,ts,'Color',cMed,'LineWidth',LWm,'DisplayName','T_s'); hold on;
plot([0 t(end)],[Tsmax Tsmax],'--','Color',cLim,'LineWidth',1.9,'DisplayName','T_{s,max} = 115 ^{\circ}C');
zxl = xline(Tcross,'LineStyle',':','Color',[.45 .45 .45],'LineWidth',1.3);
zxl.Annotation.LegendInformation.IconDisplayStyle = 'off';
text(Tcross, 34, strrep(sprintf(' t = %.1f s ', round(Tcross*10)/10),'.',','),'Rotation',90, ...
     'Color',[.3 .3 .3],'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom');
hold off;
title('Evolución térmica del bobinado en operación repetitiva (T_{ld} = 5 N\cdotm)','FontSize',15);
ylabel('T_s [^{\circ}C]','FontSize',14); xlabel('Tiempo [s]','FontSize',14);
legend('Location','northwest','FontSize',12,'Box','on');
fmtejes(FS); xlim([0 t(end)]); ylim([15 125]);
sv('5_2_5_C_comportamiento_termico.png');

fprintf('\nOK: 14 figuras generadas en\n  %s\n', outdir);

%% ===================== funciones auxiliares =====================
function fmtejes(FS)
    grid on; grid minor;
    set(gca,'FontSize',FS,'XMinorGrid','on','YMinorGrid','on', ...
            'GridAlpha',0.3,'MinorGridAlpha',0.12);
end

function setYzoom(ax, t, Y, xr)
    m = t>=xr(1) & t<=xr(2);
    v = Y(m,:); v = v(:); v = v(isfinite(v));
    if isempty(v), return; end
    lo = min(v); hi = max(v); d = hi-lo;
    if d==0, d = max(1e-6, abs(hi)); end
    ylim(ax, [lo-0.12*d, hi+0.12*d]);
end

function r = getrun(N)
    ids = Simulink.sdi.getAllRunIDs;
    r = [];
    for k = 1:numel(ids)
        rr = Simulink.sdi.getRun(ids(k));
        tok = regexp(rr.Name, 'Run\s+(\d+)', 'tokens', 'once');
        if ~isempty(tok) && str2double(tok{1}) == N
            r = rr; return;
        end
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
