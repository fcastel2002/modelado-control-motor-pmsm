%% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
%  Trabajo Final - Automática y Máquinas Eléctricas
%  Autores: Joaquín Calderón - Francisco Castel

%% Subsistema Mecánico
r = 120;            % Relación de transmisión
g = 9.80665;        % [m/s2] Aceleración de la gravedad

b_l = 0.1;          % +/- 0.03  [N.m/(rad/s)] Coef. fricción viscosa articulación
b_m = 15.0e-6;      % +/- 1%    [N.m/(rad/s)] Coef. fricción viscosa (motor y caja)

J_m = 14e-6;        % +/- 1%    [kg.m2] Momento de inercia motor y caja
J_cm = 0.0208;      %           [kg.m2] Inercia equiv. centro de masa
m = 1;              %           [kg] Masa del brazo manipulador
l_cm = 0.25;        %           [m] Longitud equiv. centro de masa
m_l = 0;            % [0, 1.5]  [kg] Masa de carga útil (nominal = 0)
l_l = 0.5;          %           [m] Longitud total extremo

% Variables derivadas
b_eq = b_m + b_l/(r^2);                 % [N.m/(rad/s)] Fricción viscosa equivalente
J_l  = (m*l_cm^2 + J_cm) + m_l*l_l^2;  % [kg.m2] Momento de inercia al eje de rotación
J_eq = J_m + J_l/r^2;                   % [kg.m2] Momento de inercia equivalente
k_l  = m*l_cm + m_l*l_l;               % [kg.m] Coeficiente k_l torque de carga

T_ld = 0;           % +/- 5.0   [N.m] Torque de perturbación por contacto

%% Subsistema Térmico
a_cu = 3.9e-3;      %           [1/°C] Coef. de aumento de Rs con Ts
C_ts = 0.818;       % +/- 1%    [W/(°C/s)] Capacitancia térmica del estator
R_ts_amb = 146.7;   % +/- 1%    [°C/W] Resistencia térmica estator-ambiente
R_sREF = 1.02;      % +/- 1%    [Ohm] Resistencia de referencia estator (T=20°C)
T_sREF = 20;        %           [°C] Temperatura de referencia
T_amb = 25;         %           [°C] Temperatura ambiente (operación)

%% Subsistema Electromagnético
P_p = 3;            % Pares de polos
lambda_r = 0.016;   % +/- 1%    [Wb] o [V/(rad/s)] Flujo magnético equiv.
L_q = 5.8e-3;       % +/- 1%    [H] Inductancia estator (eje en cuadratura)
L_d = 6.6e-3;       % +/- 1%    [H] Inductancia estator (eje directo)
L_ls = 0.8e-3;      % +/- 1%    [H] Inductancia de dispersión

%% Especificaciones de Operación
n_l_nom = 55;       % [rpm] Velocidad nominal de carga (= n_m,nom/r = 6600/120)
T_q_nom = 17;       % [N.m] Torque nominal (salida, régimen continuo)
T_q_max = 45;       % [N.m] Torque pico (salida, corta duración)

%% Lógica de Escenarios y Control
% Sobreescribe variables si ESCENARIO no es 'NOMINAL'
run('calculo_incertidumbre_y_control.m');

%% Struct para Modelo LPV (MATLAB Function en Simulink)
params = struct();

% Subsistema electromagnético
params.Ld        = L_d;
params.Lq        = L_q;
params.Lls       = L_ls;
params.Rs_ref    = R_sREF;
params.alpha_cu  = a_cu;
params.lambda_m  = lambda_r;
params.Pp        = P_p;

% Subsistema mecánico
params.J_total   = J_eq;
params.b_total   = b_eq;
params.r_red     = r;
params.kl        = k_l;
params.g         = g;

% Subsistema térmico
params.Cts       = C_ts;
params.Rts_amb   = R_ts_amb;
params.Ts_ref    = T_sREF;
