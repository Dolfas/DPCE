
% Initialization
clear
close all
clc

% Exercise 2: Implement closed-loop MPC controller (P4)

% Load model
MODEL = load('singleheater_model.mat');
load('singleheater_model.mat','A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n = size(A,1);
e_std = sqrt(e_var); % input disturbance standard deviation

% Build the functions for applying the control and reading the temperature,
% mimicking the TCLab interface
x_ss = [eye(n)-A; C]\[B*u_ss; y_ss];
c1 = ((eye(n)-A)*x_ss - B*u_ss);
c2 = (y_ss - C*x_ss);
h1 = @(x,u) A*x + B*u + Ke*e_std*randn + c1; % apply control
T1C = @(x) C*x + e_std*randn + c2; % read temperature

% Simulation parameters
T = 4000; % Experiment duration [s]
N = T/Ts; % Number of samples to collect

% Initial conditions (start at ambient temperature, i.e. equilibrium for u = 0)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n);

% Define parameters
H = 100;   % Prediction horizon
R = 0.1;    % Control weight

% Initial condition
x0 = Dx0 + x_ss;

% Initialize arrays to store data
u_mpc = zeros(1, N - 1);
x_mpc = zeros(n, N - 1);
x_mpc(:, 1) = x0;

% Initialize arrays to store output y
y_mpc = zeros(1, N - 1);
y_mpc(:, 1) = T1C(x0);

t = nan(1, N);
Dy = nan(1, N);
Du = nan(1, N);
Dx = nan(n, N);

Dx(:,1) = x0 - x_ss;
Dy(:,1) = y_mpc(:, 1) - y_ss;

% Main loop for MPC control
for k = 1:N - 1 
    % Solve MPC problem
    u = mpc_solve(x0, H, R, A, B, C, y_ss);
    u_mpc(k) = u;
    Du(:, k) = u - u_ss;
    
    % Apply control action to the plant
    Dx(:, k+1) = h1(Dx(:, k), Du(:, k));
    x_mpc(:, k+1) = Dx(:, k+1) + x_ss;
    
    % Compute output y based on updated state x_mpc
    Dy(:, k+1) = T1C(Dx(:, k+1));
    y_mpc(k) = Dy(:, k) + y_ss;

    % Update initial condition for next iteration
    x0 = x_mpc(:, k+1);
    t(k) = (k-1)*Ts;
end

% Plot results
figure('Units','normalized','Position',[0.2 0.5 0.3 0.4])
subplot(2,1,1), hold on, grid on   
title('Absolute input/output (MPC)')
plot(t(1:end-1), y_mpc,'.','MarkerSize',5)
yl = yline(y_ss,'k--');
xlabel('Time [s]')
ylabel('y [°C]')
legend(yl,'$\bar{y}$','Interpreter','latex','Location','best')
subplot(2,1,2), hold on, grid on   
stairs(t(1:end-1), u_mpc,'LineWidth',2)
yl = yline(u_ss,'k--');
yline(0,'r--')
yline(100,'r--')
xlabel('Time [s]')
ylabel('u [%]')
legend(yl,'$\bar{u}$','Interpreter','latex','Location','best');

% Plot incremental variables
figure('Units','normalized','Position',[0.5 0.5 0.3 0.4])
subplot(2,1,1), hold on, grid on   
title('Incremental input/output (MPC)')
plot(t,Dy,'.','MarkerSize',5)
xlabel('Time [s]')
ylabel('\Delta{y} [°C]')
subplot(2,1,2), hold on, grid on   
stairs(t,Du,'LineWidth',2)
yline(-u_ss,'r--')
yline(100-u_ss,'r--')
xlabel('Time [s]')
ylabel('\Delta{u} [%]')

function u0 = mpc_solve(x0, H, R, A, B, C, y_ss)
    % Compute weight matrices
    % Q = transpose(C) * C;
    
    % Initialize matrices for MPC
    W = zeros(H, H);
    Pi = zeros(H, size(A, 1));
    for j = 1:H
        for k = 1:j
            W(j, k) = C * A^(j-k) * B;
        end
        Pi(j, :) = (C * A^j);
    end
    y_ref = ones(1, H)*y_ss;

    % Compute optimal RH gain using quadprog
    F = 2 * (W' * W + R * eye(H));
    % f = 2 * x0' * Pi' * W; Makes y tend to zero
    f =  2 * (x0' * Pi' * W - y_ref * W); % Makes y tend to y_ss
    Aineq = []; bineq = [];
    Aeq = []; beq = [];
    lb = ones(H, 1) * (0); ub = ones(H, 1)*(100);

    % Set options to suppress quadprog output
    options = optimoptions('quadprog', 'Display', 'off');

    K_RH = quadprog(F, f, Aineq, bineq, Aeq, beq, lb, ub, x0, options);

    % Extract first control action
    u0 = K_RH(1);
    fprintf('control input: %f ', u0);
end
