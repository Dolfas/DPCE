% Closed-loop experiment for data collection in TCLab
%
% Copy of TCLab_openloop with some placeholder suggestions for where to
% place controller and state observer.
%
% If you see the warning 'Computation time exceeded sampling time by x
% seconds at sample k', it is because the computation time in a given
% loop took more than the sampling period Ts. Try to disable the rt_plot
% flag to fix it.
%
% Functions called: tclab.
%
% J. Miranda Lemos and Afonso Botelho, IST, May 2023
%__________________________________________________________________________

% Initialization
clear all %#ok<CLALL> 
close all
clc
tclab;

% Load model
load('singleheater_model.mat','A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n = size(A,1);
e_std = sqrt(e_var); % input disturbance standard deviation
x_ss = [eye(n)-A; C]\[B*u_ss; y_ss];

% tuning parameter
de = 1;
H = 10;
R = 0.01;
alpha = 100;

% Experiment parameters
T = 4000; % experiment duration [s]
N = T/Ts; % number of samples to collect

% Kalman filter design
% TO DO: Compute augmented matrices Ad, Bd, Cd
Ad = [A, B; zeros(1, size(A,2)), 1];
Bd = [B; 0];
Cd = [C, 0];
Qe = (Ke * e_std) * (Ke * e_std)'; % state covariance matrix
Re = e_var;
Qed = blkdiag(Qe, de);

% TO DO: L = dlqe(...);
L = dlqe(Ad, eye(size(Ad, 1)), Cd, Qed, Re);

% Initial conditions (start at ambient temperature, i.e. equilibrium for u = 0)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n); % to initialize filter

% Reference
%r = [50*ones(1,ceil(N/4)),40*ones(1,ceil(N/4)),60*ones(1,ceil(N/4)),45*ones(1,ceil(N/4))];
r = y_ss.*ones(1,N);

% Real-time plot flag. If true, plots the input and measured temperature in
% real time. If false, only plots at the end of the experiment and instead
% prints the results in the command window.
rt_plot = true;

% Initialize figure and signals
if rt_plot
    figure
    drawnow;
end
t = nan(1,N);
u = zeros(1,N);
y = zeros(1,N);
dy = nan(1,N);
du = nan(1,N);
xd_est = nan(n+1,N);
dy_error = nan(1,N);
xd_est(:,1) = [Dx0; 0.1]; % Kalman filter initialization

% String with date for saving results
timestr = char(datetime('now','Format','yyMMdd_HHmmSS'));

% Signals the start of the experiment by lighting the LED
led(1)
disp('Temperature test started.')

for k=1:N
    tic;
    % Reference changes
    ref = r(k);

    % Computes analog time
    t(k) = (k-1)*Ts;

    % Reads the sensor temperatures
    y(1,k) = T1C();

    % Compute incremental variables
    dy(:,k) = y(:,k) - ref;

    % Kalman filter correction step
    % TO DO: xd_est(:,k) = xd_est(:,k) + L*(dy(:,k) - Cd*xd_est(:,k));
    xd_est(:,k) = xd_est(:,k) + L*(dy(:,k) - Cd*xd_est(:,k));
    dy_error(:,k) = dy(:,k) - Cd*xd_est(:,k);

    % Computes the control variable to apply
    % TO DO: [...] = mpc_solve(...)
    %du(:,k)=mpc_solve(xd_est(1:n,k), H, R, A, B, C, ref, alpha);
    du(:,k)=0;
    u(:,k) = du(:,k) + u_ss;

    % Kalman filter prediction step
    % TO DO: xd_est(:,k+1) = Ad*xd_est(:,k) + Bd*du(:,k);
    xd_est(:,k+1) = Ad*xd_est(:,k) + Bd*du(:,k);

    % Applies the control variable to the plant
    h1(u(1,k));

    if rt_plot
        close all;
        f1 = figure();
        % Plots results
        subplot(2,1,1), hold on, grid on   
        plot(t(1:k),y(1,1:k),'.','MarkerSize',10)
        stairs(t,r,'r--')
        xlabel('Time [s]')
        ylabel('y [°C]')
        subplot(2,1,2), hold on, grid on   
        stairs(t(1:k),u(1,1:k),'LineWidth',2)
        xlabel('Time [s]')
        ylabel('u [%]')
        ylim([0 100]);
        drawnow;
        
        f2 = figure();
        hold on, grid on
        title('Estimation error over time');
        plot(t(1:k),dy_error(1:k));
        stairs(t,zeros(size(t,2),1),'r--');
        xlabel('Time [s]')
        ylabel('y error [ºC]')
        drawnow;

        f3 = figure();
        hold on, grid on
        title('Disturbance over time');
        plot(t(1:k),xd_est(end,1:k));
        stairs(t,zeros(size(t,2),1),'r--');
        xlabel('Time [s]')
        ylabel('disturbance')
        drawnow;

    else
        fprintf('t = %d, y1 = %.1f C, y2 = %.1f C, u1 = %.1f, u2 = %.1f\n',t(k),y(1,k),y(2,k),u(1,k),u(2,k)) %#ok<UNRCH> 
    end

    % Check if computation time did not exceed sampling time
    if toc > Ts
        warning('Computation time exceeded sampling time by %.2f s at sample %d.',toc-Ts,k)
    end
    % Waits for the begining of the new sampling interval
    pause(max(0,Ts-toc));
end

% Turns off both heaters at the end of the experiment
h1(0);
h2(0);

% Signals the end of the experiment by shutting off the LED
led(0)

disp('Temperature test complete.')

if ~rt_plot
    close all;
        f1 = figure();
        % Plots results
        subplot(2,1,1), hold on, grid on   
        plot(t(1:k),y(1,1:k),'.','MarkerSize',10)
        stairs(t,r,'r--')
        xlabel('Time [s]')
        ylabel('y [°C]')
        subplot(2,1,2), hold on, grid on   
        stairs(t(1:k),u(1,1:k),'LineWidth',2)
        xlabel('Time [s]')
        ylabel('u [%]')
        ylim([0 100]);
        drawnow;
        
        f2 = figure();
        hold on, grid on
        title('Estimation error over time');
        plot(t(1:k),dy_error(1:k));
        stairs(t,zeros(size(t,2),1),'r--');
        xlabel('Time [s]')
        ylabel('y error [ºC]')
        drawnow;

        f3 = figure();
        hold on, grid on
        title('Disturbance over time');
        plot(t(1:k),xd_est(end,1:k));
        stairs(t,zeros(size(t,2),1),'r--');
        xlabel('Time [s]')
        ylabel('disturbance')
        drawnow;

end
x_disturbance = xd_est(end,1:k);
%--------------------------------------------------------------------------

% Save figure and experiment data to file
exportgraphics(f1,['openloop_plot_',timestr,'.png'],'Resolution',300)
exportgraphics(f2,['openloop_error_plot_',timestr,'.png'],'Resolution',300)
exportgraphics(f3,['openloop_disturbance_plot_',timestr,'.png'],'Resolution',300)
save(['openloop_data_',timestr,'.mat'],'y','u','t',"dy_error","x_disturbance");

%--------------------------------------------------------------------------
% End of File


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
    u0 = z(n*(H+1)+1); 
end    