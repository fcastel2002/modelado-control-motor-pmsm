%% graficar_516.m
% Figuras del punto 5.1.6 (simulacion NL vs LTI) desde el dataset 'datos'
% exportado del Data Inspector. Toma la ULTIMA corrida: run = datos{end}.
% Guarda PNG (300 dpi) en 'figuras_516/' con los nombres del informe.
%
% Indices (verificar con el listado del final):
%   2 LTI T_s  5 LTI i_qs  6 LTI theta_m  7 LTI omega_m  8 NL I_abc
%   9 NL I_qd0s  12 NL V_abc  13 NL V_ds  14 NL delta  15 NL theta_m
%   16 NL omega_m  18 T_l  19 V_qs*

run = datos{end};
outdir = fullfile(fileparts(mfilename('fullpath')), 'figuras_516');
if ~exist(outdir,'dir'), mkdir(outdir); end

% ===================== COLORES =====================
cNL    = '#0072BD';   % NL  -> AZUL
cLTI   = '#A2142F';   % LTI -> ROJO
cErr   = '#7f7f7f';   % error -> gris
cDelta = '#ff13a6';   % angulo de carga -> ROSA
cVqs   = '#b746ff';   % consigna v_qs* -> VIOLETA
cTl    = '#64d413';   % consigna T_l   -> VERDE
cA     = '#D95319';   % fase A -> naranja
cB     = '#4DBEEE';   % fase B -> cian
cC     = '#EDB120';   % fase C -> amarillo

FS = 13; LW = 1.6; DPI = 300;
sv = @(n) exportgraphics(gcf, fullfile(outdir,n), 'Resolution', DPI);

%% sim_senales_entrada.png  (consignas)
[tv,vq]=sig(run,19); [tt,tl]=sig(run,18);
figure('Color','w');
plot(tv,vq,'Color',cVqs,'LineWidth',LW); hold on;
plot(tt,tl,'Color',cTl ,'LineWidth',LW); hold off;
title('Señales de excitación','FontSize',15);
xlabel('Tiempo [s]','FontSize',14); ylabel('Excitación','FontSize',14);
legend('v_{qs}^* [V]','T_{ld} [N\cdotm]','FontSize',12,'Location','best');
fmtejes(FS); xlim([0 2]);
sv('sim_senales_entrada.png');

%% sim_theta_m.png
[t1,nl]=sig(run,15); [t2,lt]=sig(run,6);
fig2p('Posición angular \theta_m','\theta_m [rad]', t1,nl,t2,lt, run, cNL,cLTI,cVqs,cTl,FS,LW);
sv('sim_theta_m.png');

%% sim_theta_m_error.png
[t1,nl]=sig(run,15); [~,lt]=sig(run,6);
figErr('Error de seguimiento en \theta_m','\theta_m^{NL}-\theta_m^{LTI} [rad]', t1,nl-lt, run, cErr,cVqs,cTl,FS,LW);
sv('sim_theta_m_error.png');

%% sim_omega_m.png
[t1,nl]=sig(run,16); [t2,lt]=sig(run,7);
fig2p('Velocidad angular \omega_m','\omega_m [rad/s]', t1,nl,t2,lt, run, cNL,cLTI,cVqs,cTl,FS,LW);
sv('sim_omega_m.png');

%% sim_omega_m_error.png
[t1,nl]=sig(run,16); [~,lt]=sig(run,7);
figErr('Error de seguimiento en \omega_m','\omega_m^{NL}-\omega_m^{LTI} [rad/s]', t1,nl-lt, run, cErr,cVqs,cTl,FS,LW);
sv('sim_omega_m_error.png');

%% sim_iqs.png
[t1,Iqd]=sig(run,9); nl=Iqd(:,1); [t2,lt]=sig(run,5);
fig2p('Corriente de eje q, i_{qs}','i_{qs} [A]', t1,nl,t2,lt, run, cNL,cLTI,cVqs,cTl,FS,LW,'southeast');
sv('sim_iqs.png');

%% sim_iqs_error.png
[t1,Iqd]=sig(run,9); nl=Iqd(:,1); [~,lt]=sig(run,5);
figErr('Error de seguimiento en i_{qs}','i_{qs}^{NL}-i_{qs}^{LTI} [A]', t1,nl-lt, run, cErr,cVqs,cTl,FS,LW);
sv('sim_iqs_error.png');

%% sim_Ts.png
[t1,nl]=sig(run,11); [t2,lt]=sig(run,2);
fig2p('Temperatura del estator T_s','T_s [\circC]', t1,nl,t2,lt, run, cNL,cLTI,cVqs,cTl,FS,LW);
sv('sim_Ts.png');

%% sim_v_ds.png  (1 panel)
[t1,vds]=sig(run,13);
figure('Color','w');
plot(t1,vds,'Color',cNL,'LineWidth',LW);
title('Tensión v_{ds} (ley de control mínima)','FontSize',15);
xlabel('Tiempo [s]','FontSize',14); ylabel('v_{ds} [V]','FontSize',14);
fmtejes(FS); xlim([0 2]);
sv('sim_v_ds.png');

%% sim_corrientes_abc.png  (1 panel, colores por fase)
[t1,Iabc]=sig(run,8);
figure('Color','w');
plot(t1,Iabc(:,1),'Color',cA,'LineWidth',LW); hold on;
plot(t1,Iabc(:,2),'Color',cB,'LineWidth',LW);
plot(t1,Iabc(:,3),'Color',cC,'LineWidth',LW); hold off;
title('Corrientes de fase (Park inversa)','FontSize',15);
xlabel('Tiempo [s]','FontSize',14); ylabel('i_{abc} [A]','FontSize',14);
legend('i_a','i_b','i_c','FontSize',12,'Location','best');
fmtejes(FS); xlim([0 2]);
sv('sim_corrientes_abc.png');

%% sim_tensiones_abc.png  (NUEVO)
[t1,Vabc]=sig(run,12);
figure('Color','w');
plot(t1,Vabc(:,1),'Color',cA,'LineWidth',LW); hold on;
plot(t1,Vabc(:,2),'Color',cB,'LineWidth',LW);
plot(t1,Vabc(:,3),'Color',cC,'LineWidth',LW); hold off;
title('Tensiones de fase (Park inversa)','FontSize',15);
xlabel('Tiempo [s]','FontSize',14); ylabel('v_{abc} [V]','FontSize',14);
legend('v_a','v_b','v_c','FontSize',12,'Location','best');
fmtejes(FS); xlim([0 2]);
sv('sim_tensiones_abc.png');

%% sim_angulo_carga.png  (delta ROSA + excitacion)
[t1,d]=sig(run,14);
figure('Color','w');
ax1=subplot(2,1,1);
plot(t1,d,'Color',cDelta,'LineWidth',LW);
title('Ángulo de carga \delta','FontSize',15); ylabel('\delta [\circ]','FontSize',14);
fmtejes(FS); ax1.XTickLabel=[];
ax2=subplot(2,1,2); panel_exc(run,cVqs,cTl,FS,LW);
linkaxes([ax1 ax2],'x'); xlim([0 2]);
sv('sim_angulo_carga.png');

%% utilidad: verificar indice -> nombre (descomentar)
% for i=1:run.numElements, fprintf('%2d: %s\n', i, run{i}.Name); end

%% ===================== funciones auxiliares =====================
function fmtejes(FS)
    grid on; grid minor;
    set(gca,'FontSize',FS,'XMinorGrid','on','YMinorGrid','on', ...
            'GridAlpha',0.3,'MinorGridAlpha',0.12);
    xticks(0.1:0.2:1.9);          % etiquetas en "impares": 0.1 0.3 0.5 ... (no se giran)
    xtickangle(0);                % horizontales
end

function padY(v)
    % Da margen vertical para que picos/transitorios no queden cortados.
    lo=min(v(:)); hi=max(v(:)); d=hi-lo;
    if ~isfinite(d) || d==0, d=max(1,abs(hi)); end
    ylim([lo-0.08*d, hi+0.08*d]);
end

function [t,Y] = sig(run, idx)
    v = run{idx}.Values;  t = v.Time(:);  Y = squeeze(v.Data);
    if isrow(Y), Y = Y(:); end
    if size(Y,1) ~= numel(t) && size(Y,2) == numel(t), Y = Y.'; end
end

function panel_exc(run, cVqs, cTl, FS, LW)
    [tv,vq]=sig(run,19); [tt,tl]=sig(run,18);
    plot(tv,vq,'Color',cVqs,'LineWidth',LW); hold on;
    plot(tt,tl,'Color',cTl ,'LineWidth',LW); hold off;
    xlabel('Tiempo [s]','FontSize',14); ylabel('Excitación','FontSize',13);
    lg = legend('v_{qs}^*','T_{ld}','FontSize',9,'Location','northeast');
    lg.ItemTokenSize = [12 12];   % achica las muestras de color de la leyenda
    fmtejes(FS); xlim([0 2]);
end

function fig2p(tit, ylab, t1,nl,t2,lt, run, cNL,cLTI,cVqs,cTl,FS,LW,locLeg)
    if nargin<14 || isempty(locLeg), locLeg='best'; end
    figure('Color','w');
    ax1=subplot(2,1,1);
    plot(t1,nl,'Color',cNL,'LineWidth',LW); hold on;
    plot(t2,lt,'Color',cLTI,'LineWidth',LW,'LineStyle','--'); hold off;
    title(tit,'FontSize',15); ylabel(ylab,'FontSize',14);
    legend('NL','LTI','FontSize',12,'Location',locLeg);
    fmtejes(FS); padY([nl(:);lt(:)]); ax1.XTickLabel=[];
    ax2=subplot(2,1,2); panel_exc(run,cVqs,cTl,FS,LW);
    linkaxes([ax1 ax2],'x'); xlim([0 2]);
end

function figErr(tit, ylab, t1,err, ~, cErr,~,~,FS,LW)
    figure('Color','w');                 % un solo panel, sin excitacion debajo
    plot(t1,err,'Color',cErr,'LineWidth',LW);
    title(tit,'FontSize',15); ylabel(ylab,'FontSize',14);
    xlabel('Tiempo [s]','FontSize',14);
    fmtejes(FS); padY(err); xlim([0 2]);
end
