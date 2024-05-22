% Initialization
clear
close all
clc

%% Exercise 2: Implement closed-loop MPC controller (P4)

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
T = 2000; % Experiment duration [s]
N = T/Ts; % Number of samples to collect

% Initial conditions (start at ambient temperature, i.e. equilibrium for u = 0)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n);

% Define parameters
H = 20;   % Prediction horizon
R = 0.01;    % Control weight

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

Dx(:,1) = x_mpc(:, 1) - x_ss;
Dy(:,1) = y_mpc(:, 1) - y_ss;

z0 = zeros(209,1);
% Main loop for MPC control
for k = 1:N - 1 
    % Solve MPC problem
    z = mpc_solve(x0, H, R, A, B, C, y_ss);
    Du(:,k) = z(190);
    u_mpc(k) = Du(:,k) + u_ss;

    % u_mpc(k) = u;
    % Du(:, k) = u - u_ss;
    
    % Apply control action to the plant
    x_mpc(:, k+1)  = h1(x_mpc(:, k), u_mpc(k));
    Dx(:, k+1) = x_mpc(:, k+1) - x_ss;
    
    % Compute output y based on updated state x_mpc
    y_mpc(k) = T1C(x_mpc(:, k));
    Dy(:, k) = y_mpc(k) - y_ss ;

    % Update initial condition for next iteration
    x0 = Dx(:, k+1);
    
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

function z = mpc_solve(x0, H, R, A, B, C, ref)
    % Compute the augmented matrices
    Q = C' * C;
    Q_aug = blkdiag(Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q);  % Important not to touch H 
    Q_aug = [zeros(180,9),Q_aug];
    Q_aug = [zeros(9,189);Q_aug;];
    
    R_aug = R*eye(H,H); 
    
    F = 2 .* blkdiag(Q_aug,R_aug);
    f = zeros(209,1);
    
    A_aug = blkdiag(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A);  
    zero_column = [eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9)]'.*0;
    zero_row = [eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9)].*0;
    A_aug = [A_aug, zero_column];
    A_aug = [zero_row; A_aug];
    
    B_aug = blkdiag(B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B);
    B_aug_aux = zeros(9,20);
    B_aug = [B_aug_aux;B_aug];
    
    E = [eye(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),eye(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9)]';
   
    Aeq = [A_aug - eye(189,189)  B_aug];
    beq = -E*x0;
    Aineq = []; bineq = [];

    lb = zeros(209,1); ub = zeros(209,1);
    lb(1:189,1) = -inf; ub(1:189,1)= inf;
    lb(190:209,1)= 0; ub(190:209,1) = 70;  

    % Set options to suppress quadprog output
    options = optimoptions('quadprog', 'Display', 'off');
    [z, ~, ~] = quadprog(F, f, Aineq, bineq, Aeq, beq, lb, ub, [], options);

    % Extract optimal control action and slack variables
    
    
end


