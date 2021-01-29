set(0, 'DefaultFigureWindowStyle', 'docked');
close all

% Physical Constants
global C;
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

m_n = 0.26 * C.m_0;
T_i = 300; % initial temperature (K)
tau_mn = 0.2e-12;
t = 0;
n = 1;
dt = 5e-15;
t_max = 1000*dt;

v_th = sqrt((2 * C.kb * T_i) / m_n); % thermal velocity of electrons
l_mn = v_th / tau_mn;

num_e = 10;
x_max = 200e-9; % maximum x position (nm)
y_max = 100e-9; % maximum y position (nm)

% randomly distribute electrons within the boundaries
x = x_max * rand(1, num_e);
xp = x;
y = y_max * rand(1, num_e);
yp = y;
% assign constant speed and random direction
theta = 2 * pi * rand(1, num_e);
vx = v_th * cos(theta);
vy = v_th * sin(theta);
% temp as function of electron velocities
T = @(vx, vy) (m_n .* mean((vx.^2) + (vy.^2))) / (2 * C.kb);

hold on

while t < t_max
    % compute new position
    x = x + vx * dt;
    y = y + vy * dt;
    % detect collisions and apply boundary conditions
    x_collision = (x < 0) | (x > x_max);
    x(x_collision) = mod(x(x_collision), x_max);
    xp(x_collision) = x(x_collision);
    y_collision = (y < 0) | (y > y_max);
    vy(y_collision) = -vy(y_collision);
    Temp = T(vx, vy);
    % advance clock
    t = t + dt;
    fprintf("Time: %3.3E s / %3.3E s\n", t, t_max);
    title("Temperature = " + Temp + " K");
    set(gca, 'ColorOrderIndex',1);
    plot([xp; x], [yp; y]);
    xlim([0 x_max]);
    ylim([0 y_max]);
    pause(0.042);
    
    xp = x;
    yp = y;
end


