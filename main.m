% Blade Element Theory


% Automating xfoil

function runXfoil(airfoil_name, alpha_seq, Re, output_file)
    
    input_data = [
        airfoil_name,...
        "",...
        "PANE",...
        "OPER",...
        "VISC "+num2str(Re),...
        "PACC",...
        output_file,...
        "",...
        "aseq " + alpha_seq,...
        "1",...
        "QUIT"
        ];
    
    input_filename = 'input_commands.inp';
    fid = fopen(input_filename,"w");
    for i = 1:length(input_data)
        fprintf(fid,"%s\n",input_data(i));
    end
    fclose(fid);

    system(['xfoil.exe <',input_filename]);


    data = readmatrix(output_file,'NumHeaderLines',12);

    trimmed_data = data(:,1:3);

    writematrix(trimmed_data, output_file);


end


% Physical Parameters

rho = 1.225; % Air density at sea level
c = input('Enter the choed length (c):'); % chord length
R = input('Enter the Blade Radius (R) :'); % Blade Radius
theta_inp = input('Enter the twist angle (in degrees) :');% twist angle
theta = deg2rad(theta_inp);
B = input('Enter the number of blades (B) :'); % Number of blades
V_inf = input('Enter the free stream velocity (in m/s)'); % freestream velocity
omega = input('Enter the angular velocity (in rad/s)'); % Angular velocity
N = input('Enter the number of elements, N :');


%Input

airfoil = input('Enter the airfoil name :','s');
alpha_range = input('Enter the range of alpha values :','s');
Re = input('Enter the Reynolds number value :','s');

runXfoil(airfoil,alpha_range,Re,'alpha_cl_cd.txt');

% load airfoil data
data = load('alpha_cl_cd.txt');

alpha_table = data(:,1);
Cl_table = data(:,2);
Cd_table = data(:,3);

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
    alpha_deg = rad2deg(alpha);

    Cl = interp1(alpha_table, Cl_table, alpha_deg, 'linear', 'extrap');
    Cd = interp1(alpha_table, Cd_table, alpha_deg, 'linear', 'extrap');

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








