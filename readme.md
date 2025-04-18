# Blade Element Theory

## Theory
Blade Element Theory(BET) is a mathematical approach used to analyze the aerodynamic performance of rotating blades, sucha as those found in propellers, wind turbines, helicopter rotors, and aircraft wings.

It breaks down the blade into small segments or "elements" and calculates the forces ( lift and drag ) acting on each segment based on local airflow conditions, blade geometry and angle of attack. By integrating these forces along the length of the blade, BET provides an estimate of the overall thrust, torque, and power generated or absorbed by the system.

BET treats a rotating blade (like a propeller or turbine blade) as a series of independent, infinitesimally small segments or 'elements' along its length. Each element is analyzed as a two dimensional airfoil, subject to local airflow conditions. By calculating the aerodynamic forces ( lift and drag ) on each element and summing ( or integrating ) these contributions across the entire blade, BET determines the total thrust, torque, and power produced or absorbed by the system.

The theory assumes that the flow over each blade element is primarily two-dimensional and unaffected by neighbouring elements, simplifying the complex three-dimensional flow around a rotating blade into a more tractacble problem.

[Complete Derivation of BET](https://www.aerodynamics4students.com/propulsion/blade-element-propeller-theory.php)

## Parameters involved in BET

### 1. Geometric Parameters

* Number of Blades (B)
    * The total number of blades in the rotor system.
    * Affects the total forces and solidity og the rotor. More blades increase thrust/torque but also drag and structural complexity.

* Blade Radius (R)
    * The total length of the blade from the center of rotation to the tip.
    * Determines the swept area ($\pi R^2$), which dictates how much fluid (air/water) the rotor interacts with, impacting power output.

* Hub Radius (r<sub>hub</sub>)
    * The radius from the center of the rotation to the start of the aerodynamic portion of the b lade (root or hub)
    * Defibnes the non-lifting portion of the rotor near the hub, where airflow is less effective due to structural constraints
    * Typically 5-20 % of R

* Radial Position (r) and Element width (dr)
    * r : Array of radial positions along the blade from, r<sub>hub</sub> to R.
    * dr : Width of each element (dr = (R - r<sub>hub</sub>)/N)

* Chord Length (c)
    * The width of the blade at a given radial position, perpendicular to the span.
    * Influences the lift and drag forces (proportional to c in L = 0.5 $\rho$ V<sup>2</sup>cC<sub>l</sub>)
    * Varies along the blade to optimize performance

* Twist angle ($\theta$)
    * The angle between the chord line of the blade element and the plane of rotation, varying with radius.

* Number of Elements (N)
    * The number of discrete segments the blade is divided into for analysis

### 2.Operating Condition Parameters

* Rotational Speed ($\omega$)
    * The angular velocity of the rotor
    * Determines the tangential velocity (V<sub>t</sub> = $\omega$ r), a key component of the resultant velocity affecting lift and drag

* Free stream velocity ($V_\infty$)
    * Velocity of the undisturbed fluid approaching the rotor
    * The axial component of velocity, driving the flow through the rotor and affecting power extraction or thrust

* Air density ($\rho$)
    * The density of the fluid medium

### 3.Aerodynamic Parameters

* Resultant Velocity ($V_{res}$)
    * The total velocity seen by each blade element, combining axial (V<sub>a</sub>) and tangential ( V<sub>t</sub>) components.
    * Drives lift and drag calculations
    * V<sub>res</sub> = $\sqrt{V_a^2 + V_t^2}$
    * In BET $V_a = V_\infty , V_t = \omega r$

* Inflow Angle ($\phi$)
    * The angle between the resultant velocity and the plane of rotation.
    * $\phi = tan^{-1}\frac{V_a}{V_t}$

* Angle of Attack($\alpha$)
    * The angle between the chord line and the resultant velocity vector
    * Determines $C_l and C_d$ from airfoil data, which are crucial for calculations

* Lift Coefficient ($C_l$)
    * A dimensionless factor representing the lift generated by an airfoil at a given $\alpha$

* Drag Coefficient ($C_d$)
    * A dimensionless factor representing drag at a given $\alpha$

### 4.Output Parameters

* Thrust(T)
    * The axial force produced (propeller) or extracted (turbine)

* Torque(Q)
    * The rotational force around the rotor axis
    * Drives power calculation

* Power(P)
    * The ratio of extracted power to available power in the wind
    * Measures turbine efficiency

## Working of Simulation

### Getting the input parameters

The geometric and aerodynamic parameters like Chord Length (c), Blade Radius (R), Twist angle ($\theta$), Number of Blades (B), Freestream Velocity (V $\infty$), Angular velocity($\omega$) are defined.
```
c = input('Enter the chord length (c):'); % chord length
R = input('Enter the Blade Radius (R) :'); % Blade Radius
r_hub = input('Enter the hub radius :');
theta_inp = input('Enter the twist angle (in degrees) :');% twist angle
theta = deg2rad(theta_inp);
B = input('Enter the number of blades (B) :'); % Number of blades
V_inf = input('Enter the free stream velocity (in m/s)'); % freestream velocity
omega = input('Enter the angular velocity (in rad/s)'); % Angular velocity

```

### Using X-foil to get $C_l$ and $C_d$ values

We use X-foil software to get the values of $C_l and C_d$ for different angle of attacks

```
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
```

In the above runXfoil function, the input data is written into a file called 'input_commands.inp' and that file is given to the xfoil software using the system function in matlab.

The output file from xfoil has a lot of data, but we only need $\alpha$,$C_l$ and $C_d$ values, so for getting rid of the unwanted data, the last three lines of the runXfoil function.
```
data = readmatrix(output_file,'NumHeaderLines',12);
trimmed_data = data(:,1:3);
writematrix(trimmed_data, output_file);

```

The first 12 lines in the output file will be some text which says about the airfoil, so 12 lines are left out by using readmatrix function.
Only the first three columns of the data is what is required ($\alpha$,$C_l$ and $C_d$), it is done by using the index. Finally the trimmed data is wriiten back to the output file.

### Getting airfoil data

```
airfoil = input('Enter the airfoil name :','s');
alpha_range = input('Enter the range of alpha values :','s');
Re = input('Enter the Reynolds number value :','s');
```

### Loading airfoil data
```
% Load airfoil data
runXfoil(airfoil,alpha_range,Re,'alpha_cl_cd.txt');
data = load('alpha_cl_cd.txt');
alpha_table = data(:,1);
Cl_table = data(:,2);
Cd_table = data(:,3);

```

We run the runXfoil function and get a text file named 'alpha_cl_cd.txt'. The file is loaded into MATLAB using the load function.
The $\alpha$,$C_l$ and $C_d$ are derived from the loaded data, which will be later used for interpolation.

### Discretization of the blade
The blade is split into discrete sections from the hub to the tip. Each section has a specific radial position (r), chord length (c), and twist angle ($\theta$)

```
% Discretization
r = linspace(R-r_hub,R,N);
dr = r(2)-r(1);
```
The blade is divided in to N elements and the positions of each element is defined by the matix 'r'. The elemental width dr is defined. It is the difference between two consecutive positions in the r matrix.

### Initializing Thrust and Torque
```
T = 0;  % Thrust
Q = 0;  % Torque
```

### Calculating Local velocities

```
V_a = V_inf; % Axial Velocity
V_t = omega*r(i); % Tangential Velocity
V_res = sqrt(V_a^2 + V_t^2); % Resultant velocity
```

The axial velocity is equal to the freestream velocity here. Tangential velocity is calculated here as $\omega r$.
The resultant velocity is calculated as above.\

### Calculating $C_l$ and $C_d$ values

```
phi = atan2(V_a,V_t);
alpha = phi - theta;
alpha_deg = rad2deg(alpha);
Cl = interp1(alpha_table, Cl_table, alpha_deg, 'linear', 'extrap');
Cd = interp1(alpha_table, Cd_table, alpha_deg, 'linear', 'extrap');
```
The inflow angle($\phi$), angle of attack($\alpha$) are calculated here.
$C_l$ and $C_d$ are calculated using interpolation involving the previously loaded data for the specific airfoil.

### Calculating Lift, Drag, Thrust and Torque

```
dL = 0.5 * rho * (V_res^2)* c * Cl * dr;
dD = 0.5 * rho * (V_res^2) * c * Cd * dr;

dT = B*((dL*cos(phi)) - (dD*sin(phi)));
dQ = B * r(i)*((dL*sin(phi)) + (dD*cos(phi)));
```
Lift(L), Drag(D), Torque(Q), Thrust(T) are all calculated for each element.


### Work it in Your system

Clone the repositary using the below command

```
git clone 'https://github.com/Deepak7781/Blade_Element_Theory.git'
```

And you are ready to go. Run the main.m file in MATLAB.