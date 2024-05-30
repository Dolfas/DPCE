
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
T = 500; % Experiment duration [s]
N = T/Ts; % Number of samples to collect

% Initial conditions (start at ambient temperature, i.e. equilibrium for u = 0)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n);

% Define the reference
ref_vec = [50*ones(1,ceil(N/4)),40*ones(1,ceil(N/4)),60*ones(1,ceil(N/4)),45*ones(1,ceil(N/4))];
ref = ref_vec(1);
Dref = ref - y_ss;
Dx_ref = pinv(C) * Dref;
Du_ref = pinv(B) * (Dx_ref - A * Dx_ref);

% Define parameters
H = 20;   % Prediction horizon
R = 0.01;    % Control weight
alpha = 10;

% Initial condition
x0 = Dx0 + x_ss + randn(n, 1);

x_sim = zeros(n, N - 1);
x_sim(:,1) = [0.238477245433900 0.174633394882822 -0.0552960226120982 -0.00362697099962213 -0.0248462216691613 0.0295908101449648 -0.00123115615559610 0.000991305650849106 -0.0474868449989831];
% Initialize arrays to store data
u_mpc = zeros(1, N - 1);
x_mpc = zeros(n, N - 1);
x_mpc(:, 1) = x0;

% Initialize arrays to store output y
y_mpc = zeros(1, N - 1);
y_mpc(:, 1) = T1C(x0);

y_sim = zeros(1, N - 1);

t = nan(1, N);
Dy = nan(1, N);
Du = nan(1, N);
Dx = nan(n, N);

dy = nan(1, N);
du = nan(1, N);
dx = nan(n+1, N);

dy_sim = nan(1, N);

pred_dx = nan(n, N);

Dx(:,1) = x0 - x_ss;
dx(:,1) = [Dx(:,1) - Dx_ref; e_std];

Dy(:,1) = y_mpc(:, 1) - y_ss;
dy(:,1) = Dy(:,1) - Dref;

y_sim(:,1) = 18.9149560117302;
dy_sim(:, 1)= y_sim(:, 1) - ref + 8.5;

% Compute optimal observer gain
L = dlqe(Ad, eye(size(Ad, 1)), Cd, Qed, Re);

% Main loop for MPC control
for k = 1:N - 1 
    % Reference changes
    ref = ref_vec(k);
    Dref = ref - y_ss;
    Dx_ref = pinv(C) * Dref;
    Du_ref = pinv(B) * (Dx_ref - A * Dx_ref);
     
    % Solve MPC problem
    du0 = mpc_solve(x0, H, R, A, B, C, ref,alpha);
    du(:,k) = du0;                                             
    u_mpc(k) = du(:,k) + u_ss + Du_ref;

    % Simulate the real system with input from MPC (100%)
    x_sim(:,k+1) = h1(x_sim(:,k), u_mpc(k));
    y_sim(:, k+1) = T1C(x_sim(:,k+1));
    dy_sim(:,k+1) = (y_sim(:, k+1)-ref);

    % Kalman Filter (100% works)
        %prediction
    dx1e = hd(dx(:, k), du(:,k));                  %dx1e = prediction of dx(k+1) 
    dy1e = TdC(dx1e);                              %dy1e = prediction of dy(k+1) 
        %correction
    dx(:, k+1) = dx1e + L * (dy_sim(k+1) - dy1e);  %dx1c = correction of dx(k+1) 
    
    % Extra Info
    Dx(:, k+1) = dx(1:end-1, k+1) + Dx_ref;
    x_mpc(:, k+1) = Dx(:, k+1) + x_ss;

    dy(:, k+1) = TdC(dx(:, k+1));
    Dy(:, k+1) = dy(:, k+1) + Dref;
    y_mpc(k+1) = Dy(:, k+1) + y_ss;
    
    % Update initial condition for next iteration
    x0 = dx(1:n, k+1);
    t(k) = (k-1)*Ts;
end

diferenca = y_sim - y_mpc;

% Plot results
figure('Units','normalized','Position',[0.2 0.5 0.3 0.4])
subplot(2,1,1), hold on, grid on   
title('Absolute input/output (MPC)')
plot(t(1:end), y_sim,'.','MarkerSize',5)
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
plot(t,dy_sim,'.','MarkerSize',5)
xlabel('Time [s]')
ylabel('\delta{y} [°C]')
subplot(2,1,2), hold on, grid on   
stairs(t,du,'LineWidth',2)
yline(-u_ss,'r--')
yline(100-u_ss,'r--')
xlabel('Time [s]')
ylabel('\delta{u} [%]')

figure();
plot(t,diferenca);
xlabel('Time [s]');
ylabel('error');

function u0 = mpc_solve(x0, H, R, A, B, C, ref,alpha)

    maxTemp = 55;

    n=size(A,1); % Dimension of the state
    
    % Compute the Augmented Cost matrices
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

    
    % Sparse Formulation functions 
    F = 2 .* blkdiag(Q_aug,R_aug,alpha*eye(H,H)); 
    f = zeros(size(F,2),1);

    % Augmented Dynamics matrices
    for k = 1:H
        if(k==1)
            A_aug=A;
        else
            A_aug=blkdiag(A_aug,A);
        end
    end

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
    
    for k = 1:H
        if(k==1)
            C_aug=C;
        else
            C_aug=blkdiag(C_aug,C);
        end
    end
    C_aug = [C_aug, zeros(size(C_aug,1),n)];

    E = [eye(n,n);zeros(H*n,n)];
    
    % Equality constraint matrices 
    Aeq = [A_aug - eye(size(A_aug,1),size(A_aug,2)), B_aug];
    Aeq = [Aeq, zeros(size(A_aug,1),H)]; 
    beq = -E*x0;
    
    % Inequality constraint matrices 
    g_aug = repmat(min([maxTemp-ref,0]),H,1);
    G_aug = eye(H,H);

    Aineq = [G_aug*C_aug, zeros(H,H), -eye(H,H)];
    bineq = g_aug; 

    % Lower and Higher bounds for the different variables 
    lb = -inf*ones(H*(n+2)+n,1); ub = inf*ones(H*(n+2)+n,1); % Bounds for the state
    lb(end-2*H+1:end-H,1)= -30; ub(end-2*H+1:end-H,1) = 68;  % Bounds for the control action
    lb(end-H+1:end,1)= 0; ub(end-H+1:end,1) = inf;           % Bounds for the eta

    % Set options to suppress quadprog output
    options = optimoptions('quadprog', 'Display', 'off');
    [z, ~, exitflag] = quadprog(F, f, Aineq, bineq, Aeq, beq, lb, ub, [], options);

    if exitflag == 1
        disp('Constrained optimization successful.');
    else
        disp('Constrained optimization failed.');
    end
    
    % Extract the optimal control action 
    u0=z(n*(H+1)+1);
end    