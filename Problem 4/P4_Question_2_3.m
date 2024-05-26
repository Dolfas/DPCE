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


% Main loop for MPC control
for k = 1:N - 1 
    % Solve MPC problem
    z = mpc_solve(x0, H, R, A, B, C, y_ss);
    Du(:,k) = z((H+1)*n + 1);
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
    n = 9;
    Q = C' * C;

    for i = 1:H
        if(i==1)
            Q_aug=Q;
        else
            Q_aug=blkdiag(Q_aug,Q);
        end
    end
    Q_aug = [zeros(size(Q_aug,1),n),Q_aug];
    Q_aug = [zeros(n,size(Q_aug,2));Q_aug];
    
    R_aug = R*eye(H,H); 
    
    F = 2 .* blkdiag(Q_aug,R_aug);
    f = zeros(size(F,2),1);
    
    for i = 1:H
        if(i==1)
            A_aug=A;
        else
            A_aug=blkdiag(A_aug,A);
        end
    end
    A_aug = [A_aug, zeros(size(A_aug,1),n)];
    A_aug = [zeros(n,size(A_aug,2)); A_aug];

    for i = 1:H
        if(i==1)
            B_aug=B;
        else
            B_aug=blkdiag(B_aug,B);
        end
    end
    B_aug = [zeros(n,size(B_aug,2));B_aug];
    
    E = [eye(n,n);zeros(H*n,n)];

    Aeq = [A_aug - eye(size(A_aug,1),size(A_aug,2)),  B_aug];
    beq = -E*x0;
    Aineq = []; bineq = [];

    lb = zeros(size(F,1),1); ub = zeros(size(F,1),1);
    lb(1:((H+1)*n),1) = -inf; ub(1:((H+1)*n),1)= inf;
    lb((((H+1)*n)+1):size(F,1),1)= -30; ub((((H+1)*n)+1):size(F,1),1) = 70;  

    % Set options to suppress quadprog output
    options = optimoptions('quadprog', 'Display', 'off');
    [z, ~, ~] = quadprog(F, f, Aineq, bineq, Aeq, beq, lb, ub, [], options);

    % Extract optimal control action and slack variables
    
    
end


