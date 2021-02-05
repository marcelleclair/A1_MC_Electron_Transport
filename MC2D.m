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
tau_mn = 0.2e-12; % given mean time between collisions
t = 0; % current time
n = 0; % current number of steps
dt = 5e-15; % time step duration
t_max = 1000*dt; % run for 1000 cycles
P_sca = 1 - exp(-dt / tau_mn); % scattering prob. per e- per dt
% P_sca = 0;

sig = sqrt((C.kb * T_i) / m_n); %std for maxwell-boltzmann
v_th = sqrt((2 * C.kb * T_i) / m_n); % mean thermal velocity of electrons
l_mn = v_th / tau_mn;

num_e = 10000; % number of electrons
num_disp = 10; % number of electrons displayed
x_max = 200e-9; % maximum x position (nm)
y_max = 100e-9; % maximum y position (nm)

t_slc = zeros(1, num_e); % Times since last colision
col_total = 0; % cumulative number of electron collisions and scatters
MFP = 0;
tau_calc = 0;
t_free = []; % list of all times between collisions
l_free = []; % list of all free path lengths

% randomly distribute electrons within the boundaries
x = x_max * rand(1, num_e);
xp = x;
y = y_max * rand(1, num_e);
yp = y;

% assign random speed by maxwell-boltzmann
vx = normrnd(0,sig,1,num_e);
vy = normrnd(0,sig,1,num_e);

% % assign constant speed and random direction
% theta = 2 * pi * rand(1, num_e);
% vx = v_th * cos(theta);
% vy = v_th * sin(theta);

% temp as function of electron velocities
temp = @(vx, vy) (m_n .* mean((vx.^2) + (vy.^2))) / (2 * C.kb);
T = temp(vx, vy);
Tp = T; % temperature at last step
T_avg = T;

% initialize figures
fig_traj = figure("Name", "Trajectories");
ax_traj = gca;
pbaspect([2 1 1]);
xlim(ax_traj, [0 x_max]);
ylim(ax_traj, [0 y_max]);
hold(ax_traj, 'on');
fig_temp = figure("Name", "Temperature");
ax_temp = gca;
title(ax_temp, "Semiconductor Temperature with Maxwell-Boltzmann Scattering");
ylabel(ax_temp, "Temperature (K)");
xlabel(ax_temp, "Time (s)");
hold(ax_temp, 'on');
grid(ax_temp, 'on');
fig_hist = figure("Name", "Velocity Distribution");
ax_hist = gca;

v = sqrt(vx.^2 + vy.^2);

while t < t_max
    % compute new position
    x = x + vx * dt;
    y = y + vy * dt;
    % scatter electrons
    scatter = rand(1,num_e) < P_sca;
    vx(scatter) = normrnd(0,sig,1,nnz(scatter));
    vy(scatter) = normrnd(0,sig,1,nnz(scatter));
    % detect collisions and apply boundary conditions
    x_collision = (x < 0) | (x > x_max);
    x(x_collision) = mod(x(x_collision), x_max);
    xp(x_collision) = x(x_collision);
    y_collision = (y < 0) | (y > y_max);
    vy(y_collision) = -vy(y_collision);
    
    t_slc = t_slc + dt;
    col_curr = nnz(scatter | y_collision);
    % update rolling average MFP and tau
    tau_calc = ((col_total * tau_calc) + (col_curr * mean(t_slc(scatter | y_collision))))...
        / (col_total + col_curr);
    MFP = ((col_total * MFP) + (col_curr * mean(v(scatter | y_collision) .* t_slc(scatter | y_collision))))...
        / (col_total + col_curr);
    col_total = col_total + col_curr;
    % append times for scattered or reflected electrons before resetting
    % t_free = [t_free t_slc(scatter | y_collision)];
    % l_free = [l_free (v(scatter | y_collision) .* t_slc(scatter | y_collision))];
    t_slc(scatter | y_collision) = 0;
    
    % update temp and velocity magnitude, v is left to now
    % to allow computing mean-free-path
    T = temp(vx, vy);
    v = sqrt(vx.^2 + vy.^2);
    % advance clock
    t = t + dt;
    n = n + 1;
    T_avg = T_avg + (T - T_avg)/n;
    fprintf("Time: %3.3E s / %3.3E s\n", t, t_max);
    title(ax_traj, "Temperature = " + T + " K");
    set(ax_traj, 'ColorOrderIndex',1);
    plot(ax_traj, [xp(1:num_disp); x(1:num_disp)], [yp(1:num_disp); y(1:num_disp)]);
    plot(ax_temp, [(t-dt) t], [Tp T], 'r');
    plot(ax_temp, t, T_avg, '.g');
    histogram(ax_hist, v, 'Normalization', 'probability');
    title(ax_hist, "Electron Velocity Distribution");
    xlabel(ax_hist, "Velocity (m/s)");
    ylabel(ax_hist, "Relative Frequency");
    pause(0.00001); 
    % present becomes past
    Tp = T;
    xp = x;
    yp = y;
end

%MFP = mean(l_free);
%tau_calc = mean(t_free);
fprintf("SIMULATION END AFTER %d STEPS\n", n);
fprintf("Mean Free Path = %3.3E m\n", MFP);
fprintf("Mean time between collisions = %3.3E\n", tau_calc);
