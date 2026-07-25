%% =====================================================================
%  mapa_polos_5_2_2.m
%  Punto 5.2.2: Mapa de polos en el plano s del controlador de movimiento
%  y migracion de los polos de lazo cerrado ante variacion de carga (J_eq).
%
%  Compara y grafica:
%   - Polos de la PLANTA a lazo abierto (LTI aumentado): 0, par complejo
%     electromecanico (-91.7 +/- j148), polo residual i_ds (-160.3).
%   - Polo de los LAZOS DE CORRIENTE (proporcional): -5000 rad/s.
%   - Polos de LAZO CERRADO del PID externo (diseno nominal):
%     p(s) = J_eq*s^3 + b_a*s^2 + K_sa*s + K_sia = J_eq*(s+wn)(s^2+2*zeta*wn*s+wn^2)
%     con wn = 800 rad/s, zeta = 0.75  ->  -800 y -600 +/- j529.
%   - Migracion de esos polos al variar m_l (=> J_eq) con ganancias FIJAS
%     (nominales): la inercia NO esta agendada (gain-scheduled), por lo que
%     un cambio de carga corre los polos de lazo cerrado.
%
%  Requiere parametros_sistema_completo.m en el path.
% =====================================================================
clear; clc;
parametros_sistema_completo;     % carga J_m, J_cm, m, l_cm, l_l, r, b_m, b_l,
                                 % L_q, L_d, R_sREF, lambda_r, P_p, ...

R_s0 = R_sREF*(1 + 3.9e-3*(29.5-20));         % 1.058 Ohm (T media de operacion)
kt   = 1.5*P_p*lambda_r;                       % constante de torque

% Inercia y friccion equivalentes NOMINALES (m_l = 0)
Jl_nom  = (m*l_cm^2 + J_cm) + 0*l_l^2;
Jeq_nom = J_m + Jl_nom/r^2;
beq_nom = b_m + b_l/r^2;

% ---------- 1) Planta LTI aumentado a lazo abierto ----------
A_aug = [0,                1,          0,          0;
         0,     -beq_nom/Jeq_nom,  kt/Jeq_nom,     0;
         0, -lambda_r*P_p/L_q,  -R_s0/L_q,         0;
         0,                0,          0,   -R_s0/L_d];
p_planta = eig(A_aug);

% ---------- 2) Lazo de corriente ----------
p_corriente = -5000;             % rad/s (los 3 ejes qd0)

% ---------- 3) PID externo (diseno nominal) ----------
zeta = 0.75;  wn = 800;  n = 2*zeta + 1;       % n = 2.5
ba_nom   = Jeq_nom*n*wn;
Ksa_nom  = Jeq_nom*n*wn^2;
Ksia_nom = Jeq_nom*wn^3;
p_pid_nom = roots([Jeq_nom, ba_nom, Ksa_nom, Ksia_nom]);

fprintf('\n================ POLOS (diseno nominal) ================\n');
fprintf('Planta lazo abierto (LTI aumentado), polos [rad/s]:\n');
disp(p_planta);
fprintf('Lazo de corriente: %g rad/s\n', p_corriente);
fprintf('PID lazo cerrado (nominal), polos [rad/s]:\n');
disp(p_pid_nom);
fprintf('Ganancias nominales: b_a=%.4f, K_sa=%.3f, K_sia=%.1f\n', ba_nom, Ksa_nom, Ksia_nom);

% ---------- 4) Migracion ante variacion de m_l (J_eq), ganancias fijas ----------
ml_vec  = linspace(0, 1.5, 16);
P_mig   = zeros(3, numel(ml_vec));
Jeq_vec = zeros(size(ml_vec));
for k = 1:numel(ml_vec)
    Jl_k   = (m*l_cm^2 + J_cm) + ml_vec(k)*l_l^2;
    Jeq_k  = J_m + Jl_k/r^2;
    Jeq_vec(k) = Jeq_k;
    P_mig(:,k) = roots([Jeq_k, ba_nom, Ksa_nom, Ksia_nom]);
end

fprintf('\n================ MIGRACION (m_l: 0 -> 1.5 kg, ganancias fijas) ================\n');
fprintf('%6s %12s   %-22s %-12s %-8s %-8s\n','m_l','J_eq','par complejo','polo real','wn','zeta');
for k = 1:numel(ml_vec)
    pk = P_mig(:,k);
    pr = pk(abs(imag(pk))<1e-3);              % polo(s) real(es)
    pc = pk(imag(pk)>1e-3);                    % complejo (parte imaginaria +)
    if isempty(pc)
        fprintf('%6.2f %12.3e   (3 reales: %7.1f %7.1f %7.1f)\n', ...
            ml_vec(k), Jeq_vec(k), sort(real(pk),'descend'));
    else
        wn_k = abs(pc(1));  zeta_k = -real(pc(1))/wn_k;
        fprintf('%6.2f %12.3e   %7.1f +/- j%6.1f   %8.1f   %6.1f   %5.3f\n', ...
            ml_vec(k), Jeq_vec(k), real(pc(1)), imag(pc(1)), pr(1), wn_k, zeta_k);
    end
end

% ===================== FIGURA 1: mapa de polos =====================
fig1 = figure('Color','w','Position',[100 100 1150 460]);

% (a) vista completa
subplot(1,2,1); hold on; grid on; box on;
plot(real(p_planta), imag(p_planta), 's','MarkerSize',10,'MarkerFaceColor',[0 0.447 0.741],'MarkerEdgeColor','k','LineWidth',0.8);
plot(p_corriente, 0, 's','MarkerSize',10,'MarkerFaceColor',[0.466 0.674 0.188],'MarkerEdgeColor','k','LineWidth',0.8);
plot(real(p_pid_nom), imag(p_pid_nom), 's','MarkerSize',10,'MarkerFaceColor',[0.635 0.078 0.184],'MarkerEdgeColor','k','LineWidth',0.8);
plot([0 0],[-700 700],'k-'); plot([-5300 200],[0 0],'k-');
xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');
title('Vista completa (separacion de escalas)');
legend({'Planta (lazo abierto)','Lazo de corriente (-5000)','PID lazo cerrado'}, ...
       'Location','southwest','FontSize',8);
xlim([-5300 200]); ylim([-700 700]);

% (b) zoom cerca del origen
subplot(1,2,2); hold on; grid on; box on;
plot(real(p_planta), imag(p_planta), 's','MarkerSize',11,'MarkerFaceColor',[0 0.447 0.741],'MarkerEdgeColor','k','LineWidth',0.8);
plot(real(p_pid_nom), imag(p_pid_nom), 's','MarkerSize',11,'MarkerFaceColor',[0.635 0.078 0.184],'MarkerEdgeColor','k','LineWidth',0.8);
plot([0 0],[-650 650],'k-'); plot([-950 80],[0 0],'k-');
xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');
title('Zoom: planta vs polos del PID');
legend({'Planta (lazo abierto)','PID lazo cerrado'},'Location','southwest','FontSize',8);
xlim([-950 80]); ylim([-650 650]);

sgtitle('Polos en el plano s: planta original, lazo de corriente y controlador PID de movimiento');
exportgraphics(fig1, 'imagenes/mapa_polos_5_2_2.png', 'Resolution', 150);
fprintf('\n-> Exportada: imagenes/mapa_polos_5_2_2.png\n');

% ===================== FIGURA 2: migracion =====================
fig2 = figure('Color','w','Position',[100 100 720 560]); hold on; grid on; box on;
cmap = parula(numel(ml_vec));
for k = 1:numel(ml_vec)
    plot(real(P_mig(:,k)), imag(P_mig(:,k)), 's', ...
        'MarkerSize',7,'MarkerFaceColor',cmap(k,:),'MarkerEdgeColor',cmap(k,:),'HandleVisibility','off');
end
plot(real(P_mig(:,1)),   imag(P_mig(:,1)),   'd','MarkerSize',13,'MarkerFaceColor',[0 0.447 0.741],'MarkerEdgeColor','k','LineWidth',1.2);
plot(real(P_mig(:,end)), imag(P_mig(:,end)), 's','MarkerSize',13,'MarkerFaceColor',[0.635 0.078 0.184],'MarkerEdgeColor','k','LineWidth',1.2);
plot([0 0],ylim,'k-','HandleVisibility','off');
colormap(parula); cb = colorbar; cb.Label.String = 'm_l [kg]'; caxis([0 1.5]);
xlabel('Eje Real [rad/s]'); ylabel('Eje Imaginario [rad/s]');
title('Migracion de los polos de lazo cerrado (PID) al variar m_l (J_{eq})');
legend({'m_l = 0 (nominal)','m_l = 1.5 (peor caso)'},'Location','best','FontSize',9);
exportgraphics(fig2, 'imagenes/migracion_polos_pid.png', 'Resolution', 150);
fprintf('-> Exportada: imagenes/migracion_polos_pid.png\n');
fprintf('\n================ LISTO ================\n');
