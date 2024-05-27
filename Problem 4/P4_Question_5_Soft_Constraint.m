%% Initialization
clear
close all
clc

%% Exercise 5: 

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
T = 1000; % Experiment duration [s]
N = T/Ts; % Number of samples to collect

% Initial conditions (start at ambient temperature, i.e. equilibrium for u = 0)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n);

% Define the reference
Dref = 20;
ref = Dref + y_ss;

Dx_ref = pinv(C) * Dref;
Du_ref = pinv(B) * (Dx_ref - A * Dx_ref);

% Define parameters
H = 20;   % Prediction horizon **DO NOT TOUCH**
R = 0.01;    % Control weight
alpha = 10;  %Alpha value for eta 

% Initial condition
x0 = Dx0 + x_ss;
%x0 = zeros(9,1);

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
dx = nan(n, N);


Dx(:,1) = x_mpc(:,1) - x_ss;
dx(:,1) = Dx(:,1) - Dx_ref;

Dy(:,1) = y_mpc(:, 1) - y_ss;
dy(:,1) = Dy(:,1) - Dref;

z0 = zeros(209,1);
z0(1:9,1) = dx(:,1);


% Main loop for MPC control
for k = 1:N - 1 
    % Solve MPC problem
    du0 = mpc_solve(x0, H, R, A, B, C, ref,alpha);
    du(:,k) = du0;                                              
    u_mpc(k) = du(:,k) + u_ss +  Du_ref;

    % u_mpc(k) = u;
    % Du(:, k) = u - u_ss;
    % du(:, k) = Du(:, k) - Du_ref;
    
    % Apply control action to the plant

    x_mpc(:, k+1) = h1(x_mpc(:, k), u_mpc(k));
    Dx(:, k+1) = x_mpc(:, k+1)  - x_ss;
    dx(:, k+1) = Dx(:, k+1) - Dx_ref;

    % Compute output y based on updated state x_mpc
    y_mpc(k) = T1C(x_mpc(:, k));
    Dy(:, k) = y_mpc(k) -  y_ss;
    dy(:, k) = Dy(:, k) - Dref;
    
    % Update initial condition for next iteration
    x0 = dx(:, k+1);
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
legend(yl,'$y_ss$','Interpreter','latex','Location','best')
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

function u0 = mpc_solve(x0, H, R, A, B, C, ref,alpha)

    maxTemp = 55;

    n=size(A,1);
    % Compute the augmented matrices

    Q = C' * C;
    for k = 1:H
        if(k==1)
            Q_aug=Q;
        else
            Q_aug=blkdiag(Q_aug,Q);
        end
    end
    Q_aug = [zeros(size(Q_aug,1),n),Q_aug];
    Q_aug = [zeros(n,size(Q_aug,2));Q_aug];

    R_aug = R*eye(H,H); 

    F = 2 .* blkdiag(Q_aug,R_aug,alpha*eye(H,H));
    f = zeros(size(F,2),1);

    for k = 1:H
        if(k==1)
            A_aug=A;
        else
            A_aug=blkdiag(A_aug,A);
        end
    end
    zero_column = [eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9)]'.*0;
    zero_row = [eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9),eye(9,9)].*0;
    A_aug = [A_aug, zeros(size(A_aug,1),n)];
    A_aug = [zeros(n,size(A_aug,2)); A_aug];

    for k = 1:H
        if(k==1)
            B_aug=B;
        else
            B_aug=blkdiag(B_aug,B);
        end
    end
    B_aug = [zeros(n,size(B_aug,2));B_aug];

    E = [eye(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),eye(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9),zeros(9,9)]';
    E = [eye(n,n);zeros(H*n,n)];

    Aeq = [A_aug - eye(size(A_aug,1),size(A_aug,2)), B_aug];
    Aeq = [Aeq, zeros(size(A_aug,1),H)]; % Change to Aeq
    beq = -E*x0;

    g_aug = repmat(min([maxTemp-ref,0]),H,1);
    G_aug = eye(H,H);

    for k = 1:H
        if(k==1)
            C_aug=C;
        else
            C_aug=blkdiag(C_aug,C);
        end
    end
    C_aug = [C_aug, zeros(size(C_aug,1),n)];



    Aineq = [G_aug*C_aug, zeros(H,H), -eye(H,H)];
    bineq = g_aug; 

    lb = -inf*ones(H*(n+2)+n,1); ub = inf*ones(H*(n+2)+n,1);   
    lb(end-2*H+1:end-H,1)= -30; ub(end-2*H+1:end-H,1) = 68; 
    lb(end-H+1:end,1)= 0; ub(end-H+1:end,1) = inf; % Constraint for eta

    % Set options to suppress quadprog output
    options = optimoptions('quadprog', 'Display', 'off');
    [z, ~, exitflag] = quadprog(F, f, Aineq, bineq, Aeq, beq, lb, ub, [], options);

    if exitflag == 1
        disp('Constrained optimization successful.');
    else
        disp('Constrained optimization failed.');
    end
    
    % Extract optimal control action and slack variables
    u0=z(n*(H+1)+1);
end    