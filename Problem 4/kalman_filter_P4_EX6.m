
% Initialization
clear
close all
clc

% Exercise 2: Implement closed-loop MPC controller (P4)

% Load model
MODEL = load('singleheater_model.mat');
load('singleheater_model.mat','A','B','C','Ke','e_var','y_ss','u_ss','Ts');
load('openloop_data_1.mat','y','u','t2');
n = size(A,1);
e_std = sqrt(e_var); % input disturbance standard deviation

y = y(1, 1:800)-y_ss;

% Build the functions for applying the control and reading the temperature,
% mimicking the TCLab interface
x_ss = [eye(n)-A; C]\[B*u_ss; y_ss];
c1 = ((eye(n)-A)*x_ss - B*u_ss);
c1 = c1 + 0.1*c1;
c2 = (y_ss - C*x_ss);
h1 = @(x,u) A*x + B*u + Ke*e_std*randn + c1; % apply control
T1C = @(x) C*x + e_std*randn + c2; % read temperature


% Compute covariances
Qe = (Ke * e_std) * (Ke * e_std)'; % state covariance matrix
Re = e_var;

% State augmentation
Ad = [A, B; zeros(1, size(A,2)), 1];
Bd = [B; 0];
Cd = [C, 0];
hd = @(x,u) Ad*x + Bd*u;
TdC = @(x) Cd*x;

de = 10; % tuning parameter
Qed = blkdiag(Qe, de);

% Simulation parameters
T = 4000; % Experiment duration [s]
N = T/Ts; % Number of samples to collect

% Initial conditions (start at ambient temperature, i.e. equilibrium for u = 0)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n);

% Define the reference
Dref = 0;
ref = Dref + y_ss;

Dx_ref = pinv(C) * Dref;
Du_ref = pinv(B) * (Dx_ref - A * Dx_ref);

% Define parameters
H = 10;   % Prediction horizon
R = 0.01;    % Control weight
alpha = 1e6;
eta = Dref;

% Initial condition
x0 = Dx0 + x_ss + randn(n, 1);

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

dy = nan(1, N);
du = nan(1, N);
dx = nan(n+1, N);

pred_dx = nan(n, N);

Dx(:,1) = x0 - x_ss;
dx(:,1) = [Dx(:,1) - Dx_ref; e_std];

Dy(:,1) = y_mpc(:, 1) - y_ss;
dy(:,1) = Dy(:,1) - Dref;

% Compute optimal observer gain
L = dlqe(Ad, eye(size(Ad, 1)), Cd, Qed, Re);

% Main loop for MPC control
for k = 1:N - 1 
    % Solve MPC problem
    % [u, eta] = mpc_solve(x0, H, R, alpha, A, B, C, ref, eta, k);
    u = u_ss;
    u_mpc(k) = u;
    Du(:, k) = u - u_ss;
    du(:, k) = Du(:, k) - Du_ref;
    
    % State prediction
    dx(:, k+1) = hd(dx(:, k), du(:, k));
    Dx(:, k+1) = dx(1:end-1, k+1) + Dx_ref;
    x_mpc(:, k+1) = Dx(:, k+1) + x_ss;
    
    % Compute output y based on prediction
    dy(:, k+1) = TdC(dx(:, k+1));
    Dy(:, k+1) = dy(:, k+1) + Dref;
    y_mpc(k+1) = Dy(:, k+1) + y_ss;

    % State correction
    dx(:, k+1) = dx(:, k+1) + L * (y(:, k+1) - dy(:, k+1));
    
    % Update initial condition for next iteration
    x0 = x_mpc(:, k+1);
    t(k) = (k-1)*Ts;
end

diferenca = y - dy;

% Plot results
figure('Units','normalized','Position',[0.2 0.5 0.3 0.4])
subplot(2,1,1), hold on, grid on   
title('Absolute input/output (MPC)')
plot(t(1:end), y_mpc,'.','MarkerSize',5)
yl = yline(ref,'k--');
xlabel('Time [s]')
ylabel('y [°C]')
legend(yl,'$ref$','Interpreter','latex','Location','best')
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
plot(t,dy,'.','MarkerSize',5)
xlabel('Time [s]')
ylabel('\delta{y} [°C]')
subplot(2,1,2), hold on, grid on   
stairs(t,du,'LineWidth',2)
yline(-u_ss,'r--')
yline(100-u_ss,'r--')
xlabel('Time [s]')
ylabel('\delta{u} [%]')

figure();
plot(diferenca);

function [u0, eta] = mpc_solve(x0, H, R, alpha, A, B, C, y_ss, eta, iter)
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
    ymax = 55 + eta; % soft constraint for max temp

    % Compute optimal RH gain using quadprog
    F = 2 * (W' * W + R * eye(H)); F = blkdiag(F, alpha * eye(H));
    f =  2 * (x0' * Pi' * W - y_ref * W); f = [f, zeros(1, H)];
    Aineq = eye(H) * W; bineq = ones(H, 1) * ymax - eye(H) * Pi * x0;
%     Aineq = triu(ones(H)) * W; bineq = ones(H, 1) * ymax - triu(ones(H)) * Pi * x0;
%     matrix = zeros(H, 2*H-1);
%     for i = 1:H
%         % Assign the shifted row to the matrix
%         matrix(i, i:i+H-1) = Aineq;
%     end
%     Aineq = matrix;
%     
%     disp(Aineq);
%     fprintf('Size of matrix: %d x %d\n', size(matrix, 1), size(matrix, 2));

    Aineq = blkdiag(Aineq, zeros(H, H)); bineq = [bineq; zeros(H, 1)];
    Aeq = []; beq = [];
    lb = [ones(H, 1) * (0); ones(H, 1) * (0)]; ub = [ones(H, 1)*(100); ones(H, 1) * 50];
    fprintf('Iteration: %d\n', iter);
    fprintf('MAX temp: %f\n', ymax);
%     fprintf('Size of F: %d x %d\n', size(F, 1), size(F, 2));
%     fprintf('Size of f: %d x %d\n', size(f, 1), size(f, 2));
%     fprintf('Size of Aineq: %d x %d\n', size(Aineq, 1), size(Aineq, 2));
%     fprintf('Size of bineq: %d x %d\n', size(bineq, 1), size(bineq, 2));
    
    % Set options to suppress quadprog output
    options = optimoptions('quadprog', 'Display', 'off');
    [z, ~, exitflag] = quadprog(F, f, Aineq, bineq, Aeq, beq, lb, ub, x0, options);

    % Check exitflag
    if exitflag == 1
        disp('Constrained optimization successful.');
    else
        disp('Constrained optimization failed.');
        [z, ~, ~] = quadprog(F, f, [], [], Aeq, beq, lb, ub, x0, options);
    end

    % Extract optimal control action and slack variables
    u0 = z(1);
    eta = z(H+1);
    
    fprintf('control input: %f\n', u0);
    fprintf('temp constraint margin: %f\n', eta);

end
