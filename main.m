% Blade Element Theory


% Physical Parameters

rho = 1.225; % Air density at sea level
c = 0.05; % chord length
R = 1; % Blade Radius
theta = deg2rad(10); % twist angle
B = 2; % Number of blades
V_inf = 10; % freestream velocity
omega = 300 * 2*pi/60; % Angular velocity
N = 20;


% Airfoil data
a0 = 5.7;
alpha0 = 0;
Cd0 = 0.01;
k = 0.01;

% Discretization

r = linspace(0.1*R,R,N);
dr = r(2)-r(1);

T = 0; Q = 0;


for i = 1:N
    V_a = V_inf;
    V_t = omega*r(i);
    V_res = sqrt(V_a^2 + V_t^2);

    phi = atan2(V_a,V_t);
    alpha = phi - theta;

    Cl = a0 * alpha;
    Cd = Cd0 + (k*(Cl^2));

    dL = 0.5 * rho * (V_res^2)* c * Cl * dr;
    dD = 0.5 * rho * (V_res^2) * c * Cd * dr;

    dT = B*((dL*cos(phi)) - (dD*sin(phi)));
    dQ = B * r(i)*((dL*sin(phi)) + (dD*cos(phi)));

    T = T + dT;
    Q = Q + dQ;
end

P = omega * Q;
eta = (T*V_inf)/P;


fprintf('Total Thrust:     %.2f N\n', T);
fprintf('Total Torque:     %.2f Nm\n', Q);
fprintf('Power Required:   %.2f W\n', P);
fprintf('Prop Efficiency:  %.2f\n', eta);








