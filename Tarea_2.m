clc; clear; close all;
%% Batimetría
fname = 'REU2004bathy.txt';  
A = readmatrix(fname, 'CommentStyle','%'); 
x = A(:,1);  z = A(:,2);

g   = 9.81;            
h = -z; 

%% === Constantes para rotura (Baldock 1998) ===
rho = 1000;           % kg/m^3
Bcoef = 1.0;          % B ~ 1
h_floor = 0.02;       % piso numérico de profundidad

%% === Con rotura - Primer caso (R39) ===
S = load('R39.mat'); 
y1 = S.R39.LWF.H;
x1 = S.R39.xreal;

fp    = 1/8;                 % Hz
Hrms0 = 0.37;                % m

omega = 2*pi*fp;
k  = arrayfun(@(hh) k_dispersion(omega, hh, g), h);        
C  = omega ./ k;                                           
Cg = 0.5*C .* (1 + (2*k.*h) ./ sinh(2*k.*h));

% ---------- Baldock 98 ----------
Tp = 1/fp; L0 = 1.56*Tp^2; S0 = Hrms0 / L0;
gamma = 0.39 + 0.56*tanh(33*S0);
Aeps  = 0.25 * rho * g * fp * Bcoef;

N = numel(x);
Hrms = nan(N,1); F = nan(N,1);
Hrms(1) = Hrms0;
E0 = 0.125*rho*g*Hrms0^2; 
F(1) = E0*Cg(1);

for i = 1:N-1
    hi = max(h(i), h_floor);
    H  = max(Hrms(i), 1e-6);
    Hb = gamma * hi;
    eps_b = Aeps * exp( - (Hb/H)^2 ) * ( Hb^2 + H^2 );
    dx = x(i+1)-x(i);
    F(i+1) = max(F(i) - eps_b*dx, 0);
    Hrms(i+1) = sqrt( max( 8*F(i+1)/(rho*g*Cg(i+1)), 0) );
end
Hrms(h<=0) = NaN;

% === FIGURA 1: Baldock — R39 ARRIBA ===
figure(1)
subplot(3,1,1);
plot(x, Hrms, 'b','LineWidth',2); hold on
plot(x1, y1, 'ko', 'MarkerSize',6, 'MarkerFaceColor','r','DisplayName','Puntos'); 
grid on; box on;
xlabel('x [m]'); ylabel('H_{rms} [m]');
title('H_{rms}(x) con rotura — Baldock (1998) V/S R39');
xlim([min(x) max(x)]); ylim([0 1]);
legend('H_{rms}','R39','Location','best');
text(0.02, 0.98, 'a)', 'Units','normalized','FontWeight','bold', ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');

%% === Con rotura - Segundo caso (R40) ===
S = load('R40.mat'); 
y1 = S.R40.LWF.H;
x1 = S.R40.xreal;

fp    = 1/4;                 % Hz
Hrms0 = 0.40;                % m

omega = 2*pi*fp;
k  = arrayfun(@(hh) k_dispersion(omega, hh, g), h);        
C  = omega ./ k;                                           
Cg = 0.5*C .* (1 + (2*k.*h) ./ sinh(2*k.*h));

Tp = 1/fp; L0 = 1.56*Tp^2; S0 = Hrms0 / L0;
gamma = 0.39 + 0.56*tanh(33*S0);
Aeps  = 0.25 * rho * g * fp * Bcoef;

N = numel(x);
Hrms = nan(N,1); F = nan(N,1);
Hrms(1) = Hrms0;
E0 = 0.125*rho*g*Hrms0^2; 
F(1) = E0*Cg(1);

for i = 1:N-1
    hi = max(h(i), h_floor);
    H  = max(Hrms(i), 1e-6);
    Hb = gamma * hi;
    eps_b = Aeps * exp( - (Hb/H)^2 ) * ( Hb^2 + H^2 );
    dx = x(i+1)-x(i);
    F(i+1) = max(F(i) - eps_b*dx, 0);
    Hrms(i+1) = sqrt( max( 8*F(i+1)/(rho*g*Cg(i+1)), 0) );
end
Hrms(h<=0) = NaN;

% === FIGURA 1: Baldock — R40 AL MEDIO ===
figure(1)
subplot(3,1,2);
plot(x, Hrms, 'b','LineWidth',2); hold on
plot(x1, y1, 'ko', 'MarkerSize',6, 'MarkerFaceColor','r','DisplayName','Puntos');
grid on; box on;
xlabel('x [m]'); ylabel('H_{rms} [m]');
title('H_{rms}(x) con rotura — Baldock (1998) V/S R40');
xlim([min(x) max(x)]); ylim([0 1]);
legend('H_{rms}','R40','Location','best');
text(0.02, 0.98, 'b)', 'Units','normalized','FontWeight','bold', ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');

% === FIGURA 1: BATIMETRÍA ABAJO ===
figure(1)
subplot(3,1,3);
plot(x, z, 'k','LineWidth',1.4);
grid on; box on;
xlabel('x [m]'); ylabel('z [m]'); yline(0,'--k');
title('Batimetría a lo largo del perfil');
text(0.02, 0.95, 'c)', 'Units','normalized','FontWeight','bold', ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');
xlim([min(x) max(x)]);

%% ================= FIGURA 2: Lippmann (1996) =================
alpha_deg = 10;               % ángulo del frente del roller (5–15° típico)
tanA = tand(alpha_deg); 
cosA = cosd(alpha_deg);

% ---------- R39 (arriba) ----------
S = load('R39.mat'); 
y1 = S.R39.LWF.H;   x1 = S.R39.xreal;
fp    = 1/8;        Hrms0 = 0.37;

omega = 2*pi*fp;
k  = arrayfun(@(hh) k_dispersion(omega, hh, g), h);
C  = omega ./ k;                                           
Cg = 0.5*C .* (1 + (2*k.*h) ./ sinh(2*k.*h));

Tp = 1/fp; L0 = 1.56*Tp^2; S0 = Hrms0 / L0;
gamma = 0.39 + 0.56*tanh(33*S0);

N = numel(x);
HrmsL = nan(N,1); Ft = nan(N,1);
HrmsL(1) = Hrms0;
Ew0 = 0.125*rho*g*Hrms0^2;

% roller offshore ~ 0 (Qb tiny)
Hb = gamma * max(h, h_floor);
Lw = 2*pi ./ max(k,1e-8);
Qb = exp( -(Hb./HrmsL(1)).^2 ); Qb(:) = Qb(1);   % solo para inicializar
Aroll = Qb .* (Hb.^2) ./ (4*max(h,h_floor)*tanA);
Er = 0.5 * rho .* (C.^2) .* (Aroll ./ Lw);

Ft(1) = Ew0*Cg(1) + Er(1)*C(1);

for i = 1:N-1
    hi = max(h(i), h_floor);
    Hb_i = gamma * hi;
    % fracción rompiendo en i (Rayleigh): Qb = exp(-(Hb/H)^2)
    Qb_i = exp( -(Hb_i / max(HrmsL(i),1e-6))^2 );
    % disipación por cizalle (ponderada por Qb)
    epsL = 0.25 * rho * g * fp * C(i) * cosA * (Hb_i^2 / hi) * Qb_i;

    dx = x(i+1)-x(i);
    Ft(i+1) = max(Ft(i) - epsL*dx, 0);

    % resolver H en i+1 con iteración fija (Er depende de H via Qb)
    hi1 = max(h(i+1), h_floor); Hb1 = gamma * hi1;
    Hguess = HrmsL(i);
    for it = 1:6
        Qb1 = exp( -(Hb1 / max(Hguess,1e-6))^2 );
        A1  = Qb1 * (Hb1^2) / (4*hi1*tanA);
        Er1 = 0.5 * rho * (C(i+1)^2) * (A1 / Lw(i+1));
        Ew1 = max(Ft(i+1) - Er1*C(i+1), 0) / max(Cg(i+1),1e-8);
        Hnew = sqrt( max(8*Ew1/(rho*g),0) );
        if abs(Hnew - Hguess) < 1e-6, break; end
        Hguess = 0.6*Hguess + 0.4*Hnew;    
    end
    HrmsL(i+1) = Hnew;
end
HrmsL(h<=0) = NaN;

figure(2)
subplot(3,1,1);
plot(x, HrmsL, 'g','LineWidth',2, 'DisplayName','Lippmann''96 (roller)'); hold on
plot(x1, y1, 'ko', 'MarkerSize',6, 'MarkerFaceColor','r','DisplayName','R39'); 
grid on; box on;
xlabel('x [m]'); ylabel('H_{rms} [m]');
title('H_{rms}(x) con rotura — Lippmann (1996) vs R39');
xlim([min(x) max(x)]); ylim([0 1]);
legend('Location','best');
text(0.02, 0.98, 'a)', 'Units','normalized','FontWeight','bold', ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');

% ---------- R40 (medio) ----------
S = load('R40.mat'); 
y1 = S.R40.LWF.H;   x1 = S.R40.xreal;
fp    = 1/4;        Hrms0 = 0.40;

omega = 2*pi*fp;
k  = arrayfun(@(hh) k_dispersion(omega, hh, g), h);
C  = omega ./ k;                                           
Cg = 0.5*C .* (1 + (2*k.*h) ./ sinh(2*k.*h));

Tp = 1/fp; L0 = 1.56*Tp^2; S0 = Hrms0 / L0;
gamma = 0.39 + 0.56*tanh(33*S0);

N = numel(x);
HrmsL = nan(N,1); Ft = nan(N,1);
HrmsL(1) = Hrms0;
Ew0 = 0.125*rho*g*Hrms0^2;

Hb = gamma * max(h, h_floor);
Lw = 2*pi ./ max(k,1e-8);

% init roller ~ 0 offshore
Qb = exp( -(Hb./HrmsL(1)).^2 ); Qb(:) = Qb(1);
Aroll = Qb .* (Hb.^2) ./ (4*max(h,h_floor)*tanA);
Er = 0.5 * rho .* (C.^2) .* (Aroll ./ Lw);

Ft(1) = Ew0*Cg(1) + Er(1)*C(1);

for i = 1:N-1
    hi = max(h(i), h_floor);
    Hb_i = gamma * hi;
    Qb_i = exp( -(Hb_i / max(HrmsL(i),1e-6))^2 );
    epsL = 0.25 * rho * g * fp * C(i) * cosA * (Hb_i^2 / hi) * Qb_i;

    dx = x(i+1)-x(i);
    Ft(i+1) = max(Ft(i) - epsL*dx, 0);

    hi1 = max(h(i+1), h_floor); Hb1 = gamma * hi1;
    Hguess = HrmsL(i);
    for it = 1:6
        Qb1 = exp( -(Hb1 / max(Hguess,1e-6))^2 );
        A1  = Qb1 * (Hb1^2) / (4*hi1*tanA);
        Er1 = 0.5 * rho * (C(i+1)^2) * (A1 / Lw(i+1));
        Ew1 = max(Ft(i+1) - Er1*C(i+1), 0) / max(Cg(i+1),1e-8);
        Hnew = sqrt( max(8*Ew1/(rho*g),0) );
        if abs(Hnew - Hguess) < 1e-6, break; end
        Hguess = 0.6*Hguess + 0.4*Hnew;
    end
    HrmsL(i+1) = Hnew;
end
HrmsL(h<=0) = NaN;

figure(2)
subplot(3,1,2);
plot(x, HrmsL, 'g','LineWidth',2, 'DisplayName','Lippmann''96 (roller)'); hold on
plot(x1, y1, 'ko', 'MarkerSize',6, 'MarkerFaceColor','r','DisplayName','R40');
grid on; box on;
xlabel('x [m]'); ylabel('H_{rms} [m]');
title('H_{rms}(x) con rotura — Lippmann (1996) vs R40');
xlim([min(x) max(x)]); ylim([0 1]);
legend('Location','best');
text(0.02, 0.98, 'b)', 'Units','normalized','FontWeight','bold', ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');

% ---------- Batimetría abajo ----------
figure(2)
subplot(3,1,3);
plot(x, z, 'k','LineWidth',1.4);
grid on; box on;
xlabel('x [m]'); ylabel('z [m]'); yline(0,'--k');
title('Batimetría a lo largo del perfil');
text(0.02, 0.95, 'c)', 'Units','normalized','FontWeight','bold', ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');
xlim([min(x) max(x)]);

%% ====== Función auxiliar: dispersión lineal (Newton-Raphson) ======
function k = k_dispersion(omega, h, g)
    if h <= 0, k = eps; return; end
    k = max(omega^2/g, 1e-6);                  
    for it = 1:30
        th = tanh(k*h);
        f  = g*k*th - omega^2;
        df = g*th + g*k*h*(1 - th^2);
        step = f/df;
        k = k - step;
        if abs(step) < 1e-12, break; end
    end
    k = max(k, 1e-8);
end
