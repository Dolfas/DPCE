clear , clc, close all;


%% Exercise 1 - (2nd Plant) - R = 1

% Define system matrices
A = 0.8;
B = 1;
C = 1;

% Define weight matrices
Q = transpose(C) * C;
R = 1;

% Compute optimal LQ gain using dlqr
[K_LQ, S, ~] = dlqr(A, B, Q, R);

eigenvalues = eig(A - B*K_LQ);

%% Exercise 2 - (2nd Plant) - R = 1

% Initialize horizon H values
H_values = 1:30;

% Initialize arrays to store optimal RH gains
K_RH_values = zeros(size(H_values));
RH_eigenvalues = zeros(size(H_values));

% Compute optimal RH gain for each horizon value
for i = 1:length(H_values)
    H = H_values(i);
    
    % Compute matrices W and Pi
    W = zeros(H);
    Pi = zeros(1, H);
    for j = 1:H
        for k = 1:j
            W(j, k) = C * A^(j-k) * B;
        end
        Pi(j) = C * A^(j);
    end

    % Compute optimal RH gain
    M = W' * W + R * eye(H);
    K_RH = inv(M) * W' * Pi';
    K_RH_values(i) = K_RH(1);
    RH_eigenvalues(i) = eig(A - B*K_RH(1));
end

% Plot optimal RH gain as a function of H
figure;
plot(H_values, K_RH_values, 'b-o', 'LineWidth', 2);
hold on;
grid on;
optimal_LQ_gain = repmat(K_LQ, size(H_values));
plot(H_values, optimal_LQ_gain, 'r--', 'LineWidth', 2);
xlabel('Horizon (H)');
ylabel('Optimal Gain');
legend('RH Gain', 'Optimal LQ Gain');
title('Optimal LQ Gain vs RH 1st Gain (R = 1)');

%% Exercise 3 - (2nd Plant) - R = 1


figure;
plot(H_values, RH_eigenvalues, 'b-o', 'LineWidth', 2);
hold on;
grid on;
yline(1, 'r--', 'LineWidth', 2);
xlabel('Horizon (H)');
ylabel('Eigenvalues');
title('Closed-Loop eigenvalues (R = 1)');

figure;
zplane([], RH_eigenvalues'), grid on

%% Exercise 1 - (2nd Plant) - R = 10

% Define system matrices
A = 0.8;
B = 1;
C = 1;

% Define weight matrices
Q = transpose(C) * C;
R = 10;

% Compute optimal LQ gain using dlqr
[K_LQ, S, ~] = dlqr(A, B, Q, R);

eigenvalues = eig(A - B*K_LQ);

%% Exercise 2 - (2nd Plant) - R = 10

% Initialize horizon H values
H_values = 1:30;

% Initialize arrays to store optimal RH gains
K_RH_values = zeros(size(H_values));
RH_eigenvalues = zeros(size(H_values));

% Compute optimal RH gain for each horizon value
for i = 1:length(H_values)
    H = H_values(i);
    
    % Compute matrices W and Pi
    W = zeros(H);
    Pi = zeros(1, H);
    for j = 1:H
        for k = 1:j
            W(j, k) = C * A^(j-k) * B;
        end
        Pi(j) = C * A^(j);
    end

    % Compute optimal RH gain
    M = W' * W + R * eye(H);
    K_RH = inv(M) * W' * Pi';
    K_RH_values(i) = K_RH(1);
    RH_eigenvalues(i) = eig(A - B*K_RH(1));
end

% Plot optimal RH gain as a function of H
figure;
plot(H_values, K_RH_values, 'b-o', 'LineWidth', 2);
hold on;
grid on;
optimal_LQ_gain = repmat(K_LQ, size(H_values));
plot(H_values, optimal_LQ_gain, 'r--', 'LineWidth', 2);
xlabel('Horizon (H)');
ylabel('Optimal Gain');
legend('RH Gain', 'Optimal LQ Gain');
title('Optimal LQ Gain vs RH 1st Gain (R = 10)');

%% Exercise 3 - (2nd Plant) - R = 10

figure;
plot(H_values, RH_eigenvalues, 'b-o', 'LineWidth', 2);
hold on;
grid on;
yline(1, 'r--', 'LineWidth', 2);
xlabel('Horizon (H)');
ylabel('Eigenvalues');
title('Closed-Loop eigenvalues (R = 10)');

figure;
zplane([], RH_eigenvalues'), grid on

%% Exercise 1 - (2nd Plant) - R = 100

% Define system matrices
A = 0.8;
B = 1;
C = 1;

% Define weight matrices
Q = transpose(C) * C;
R = 100;

% Compute optimal LQ gain using dlqr
[K_LQ, S, ~] = dlqr(A, B, Q, R);

eigenvalues = eig(A - B*K_LQ);

%% Exercise 2 - (2nd Plant) - R = 100

% Initialize horizon H values
H_values = 1:30;

% Initialize arrays to store optimal RH gains
K_RH_values = zeros(size(H_values));
RH_eigenvalues = zeros(size(H_values));

% Compute optimal RH gain for each horizon value
for i = 1:length(H_values)
    H = H_values(i);
    
    % Compute matrices W and Pi
    W = zeros(H);
    Pi = zeros(1, H);
    for j = 1:H
        for k = 1:j
            W(j, k) = C * A^(j-k) * B;
        end
        Pi(j) = C * A^(j);
    end

    % Compute optimal RH gain
    M = W' * W + R * eye(H);
    K_RH = inv(M) * W' * Pi';
    K_RH_values(i) = K_RH(1);
    RH_eigenvalues(i) = eig(A - B*K_RH(1));
end

% Plot optimal RH gain as a function of H
figure;
plot(H_values, K_RH_values, 'b-o', 'LineWidth', 2);
hold on;
grid on;
optimal_LQ_gain = repmat(K_LQ, size(H_values));
plot(H_values, optimal_LQ_gain, 'r--', 'LineWidth', 2);
xlabel('Horizon (H)');
ylabel('Optimal Gain');
legend('RH Gain', 'Optimal LQ Gain');
title('Optimal LQ Gain vs RH 1st Gain (R = 100)');

%% Exercise 3 - (2nd Plant) - R = 100

figure;
plot(H_values, RH_eigenvalues, 'b-o', 'LineWidth', 2);
hold on;
grid on;
yline(1, 'r--', 'LineWidth', 2);
xlabel('Horizon (H)');
ylabel('Eigenvalues');
title('Closed-Loop eigenvalues (R = 100)');

figure;
zplane([], RH_eigenvalues'), grid on
