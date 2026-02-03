function AO = calcular_matriz_AO(x_op, params)
    % Entradas:
    % x_op: Vector del punto de operación [theta_mO; w_mO; i_qsO; i_dsO; i_0sO; TsO]
    % params: Estructura con constantes (Ld, Lq, Rs_ref, etc.)
    
    % 1. Desempaquetar el Punto de Operación (Operating Point)
    theta_mO = x_op(1);
    w_mO     = x_op(2);
    i_qsO    = x_op(3);
    i_dsO    = x_op(4);
    i_0sO    = x_op(5);
    % Nota: TsO (Temperatura) suele ser una entrada exógena o estado lento.
    % Si la linealización incluye la dinámica térmica, se usa x_op(6).
    % Si no, TsO es un parámetro que varía lento. Asumimos estado x_op(6) 
    % por tu matriz de 6x6.
    TsO      = x_op(6); 

    % 2. Desempaquetar parámetros (o usar params.Variable)
    % Esto asume que tienes una estructura 'params' definida en el workspace
    % o creada via Mask.
    Ld      = params.Ld;
    Lq      = params.Lq;
    Lls     = params.Lls;
    Rs_ref  = params.Rs_ref;
    alpha_cu= params.alpha_cu;
    lambda_m= params.lambda_m;
    Pp      = params.Pp;
    J_eq    = params.J_total; % Inercia total
    b_eq    = params.b_total; % Fricción total
    r       = params.r_red;
    kl      = params.kl;      % Constante de carga (m*g*l_cm...)
    g       = params.g;
    Cts     = params.Cts;
    Rts_amb = params.Rts_amb;
    
    % 3. Cálculo auxiliar de Rs en el punto de operación
    % (Necesario para los términos que no son derivadas parciales de R)
    % Ec. [cite_start]3.9 de la guía [cite: 103]
    % Asumiendo T_ref = 20 grados (o el valor que tengas en params)
    delta_T = TsO - params.Ts_ref;
    Rs_op   = Rs_ref * (1 + alpha_cu * delta_T);

    % 4. Cálculo de Gamma_th (Elemento A(6,6))
    % Definición dada en tu ecuación
    term_currents = 2*(i_0sO^2) + i_dsO^2 + i_qsO^2;
    Gamma_th = -((1/Rts_amb) - (3/2)*Rs_ref*alpha_cu*term_currents) / Cts;

    % 5. Construcción de la Matriz AO (6x6)
    % Inicializamos en ceros para eficiencia
    AO = zeros(6,6);

    % Fila 1: Derivada de Theta (Velocidad)
    AO(1,2) = 1;

    % Fila 2: Derivada de Omega (Mecánica)
    % Termino d(Tm)/d(theta) (Gravedad linealizada)
    AO(2,1) = -(g * kl * cos(theta_mO / r)) / (J_eq * r^2); 
    % Fricción viscosa
    AO(2,2) = -b_eq / J_eq;
    % Acoplamiento corriente q (Torque)
    AO(2,3) = (1.5 * Pp * (lambda_m + i_dsO*(Ld - Lq))) / J_eq;
    % Acoplamiento corriente d (Torque de reluctancia)
    AO(2,4) = (1.5 * Pp * i_qsO * (Ld - Lq)) / J_eq;

    % Fila 3: Derivada de i_qs (Eje q)
    % Acoplamiento w (f.c.e.m.)
    AO(3,2) = -(Pp * (lambda_m + Ld * i_dsO)) / Lq;
    % Resistencia / Inductancia
    AO(3,3) = -Rs_op / Lq;
    % Acoplamiento eje d
    AO(3,4) = -(Ld * Pp * w_mO) / Lq;
    % Sensibilidad térmica de R
    AO(3,6) = -(Rs_ref * alpha_cu * i_qsO) / Lq;

    % Fila 4: Derivada de i_ds (Eje d)
    % Acoplamiento eje q
    AO(4,2) = (Lq * Pp * i_qsO) / Ld; % OJO: En tu imagen dice Lq*Pp*i_qs0 / Ld en pos(4,2)?
    % Nota: Revisa si este termino depende de w o de i_qs. 
    % La derivada de (w*Lq*iq) respecto a w es Lq*iq. Correcto.
    
    AO(4,3) = (Lq * Pp * w_mO) / Ld;
    AO(4,4) = -Rs_op / Ld;
    AO(4,6) = -(Rs_ref * alpha_cu * i_dsO) / Ld;

    % Fila 5: Derivada de i_0s (Secuencia cero)
    % Si el neutro flota, i_0s suele ser 0 y desacoplada, pero incluimos tu ec.
    AO(5,5) = -Rs_op / Lls;
    AO(5,6) = -(Rs_ref * alpha_cu * i_0sO) / Lls;

    % Fila 6: Derivada de Temperatura
    % Joule losses derivatives
    AO(6,3) = (3 * Rs_op * i_qsO) / Cts;
    AO(6,4) = (3 * Rs_op * i_dsO) / Cts;
    AO(6,5) = (6 * Rs_op * i_0sO) / Cts; % Factor 6 según tu fórmula (3 fases * 2?)
    AO(6,6) = Gamma_th;

end