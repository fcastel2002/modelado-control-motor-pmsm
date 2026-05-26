%% =====================================================================
%  metricas_5_1_6b.m
%  Métricas de transitorio CORRECTAS para el Sub-ítem 5.1.6.b
%
%  Las métricas se miden sobre la VELOCIDAD omega_m del modelo LTI
%  aumentado (caso nominal i_ds(0)=0). omega_m está gobernada por el par
%  de polos complejos SUBAMORTIGUADOS hallado en 5.1.3
%  (lambda_{2,3} = -91.7 +/- j148, zeta = 0.527, w_n = 174 rad/s),
%  por lo tanto SU RESPUESTA PRESENTA SOBREPICO.
%
%  (La tabla anterior trataba a i_qs como un primer orden que se
%   establecía en v_qs*/R_s = 18.53 A con M_p = 0%: eso correspondía a un
%   modelo DESACOPLADO que NO es el que se simula. En el modelo real i_qs
%   es un PICO de torque que decae a ~0 al aparecer la fuerza
%   contraelectromotriz, y omega_m alcanza la velocidad de vacío con
%   sobrepico.)
%
%  USO:  correr primero  simulacion_dinamica_516.m  (deja en el workspace
%        res_lti, res_nl y t_sim) y luego este script.
% =====================================================================

if ~exist('res_lti','var') || ~exist('t_sim','var')
    error(['Faltan variables en el workspace. Corré primero ', ...
           'simulacion_dinamica_516.m para generar res_lti y t_sim.']);
end

wm  = res_lti{1}(:,2);     % omega_m [rad/s]  (caso i_ds(0)=0)
iqs = res_lti{1}(:,3);     % i_qs    [A]

% Instante de inicio y etiqueta de cada uno de los 10 transitorios
% (coinciden con vqs_func y tld_func de simulacion_dinamica_516.m)
t_ev = [0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9];
evt  = {'vqs* -> +19.6V','Tld -> +6.28','Tld -> -6.28', ...
        'vqs* -> 0','Tld -> 0','vqs* -> -19.6V', ...
        'Tld -> +6.28','Tld -> -6.28','vqs* -> 0','Tld -> 0'};
seg  = 0.2;     % duración de cada segmento [s]

fprintf('\n=============== MÉTRICAS sobre omega_m  (LTI, i_ds(0)=0) ===============\n');
fprintf('%-3s %-16s %10s %8s %8s %8s %8s\n', ...
        'Tr','Evento','wm_ss','Mp[%]','tp[ms]','tr[ms]','ts[ms]');
for k = 1:numel(t_ev)
    idx = (t_sim >= t_ev(k)) & (t_sim < t_ev(k)+seg);
    tt  = t_sim(idx) - t_ev(k);
    yy  = wm(idx);
    y0  = yy(1);                       % valor inicial del segmento
    yf  = mean(yy(end-19:end));        % valor de régimen del segmento
    if abs(yf - y0) < 1e-3
        fprintf('%-3d %-16s %10.1f %8s %8s %8s %8s\n', ...
                k, evt{k}, yf, '--','--','--','--');
        continue;
    end
    % Se normaliza como un escalón de 0 a (yf-y0) -> robusto entre versiones
    S = stepinfo(yy - y0, tt, yf - y0, 'SettlingTimeThreshold', 0.01);
    fprintf('%-3d %-16s %10.1f %8.1f %8.1f %8.1f %8.1f\n', ...
            k, evt{k}, yf, S.Overshoot, S.PeakTime*1e3, ...
            S.RiseTime*1e3, S.SettlingTime*1e3);
end

fprintf('\n=============== Pico y régimen de i_qs  (LTI, i_ds(0)=0) ==============\n');
fprintf('%-3s %-16s %12s %12s\n','Tr','Evento','iqs_pico[A]','iqs_reg[A]');
for k = 1:numel(t_ev)
    idx = (t_sim >= t_ev(k)) & (t_sim < t_ev(k)+seg);
    s   = iqs(idx);
    [~,im] = max(abs(s));
    fprintf('%-3d %-16s %12.2f %12.3f\n', ...
            k, evt{k}, s(im), mean(s(end-19:end)));
end

% (Opcional) lo mismo sobre el modelo NL, para verificar coincidencia
if exist('res_nl','var')
    wm_nl = res_nl{1}(:,2);
    fprintf('\n=============== Verificación: sobrepico de omega_m (NL) ===============\n');
    for k = [1 4 6 9]     % solo escalones de v_qs* (donde el sobrepico es relevante)
        idx = (t_sim >= t_ev(k)) & (t_sim < t_ev(k)+seg);
        tt  = t_sim(idx) - t_ev(k);
        yy  = wm_nl(idx);  y0 = yy(1);  yf = mean(yy(end-19:end));
        S = stepinfo(yy - y0, tt, yf - y0, 'SettlingTimeThreshold', 0.01);
        fprintf('Tr %2d (%s): wm_ss=%7.1f  Mp=%5.1f%%  ts=%5.1f ms\n', ...
                k, evt{k}, yf, S.Overshoot, S.SettlingTime*1e3);
    end
end
fprintf('\n======================================================================\n');
