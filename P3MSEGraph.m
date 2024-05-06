n=[1:19];
MSE1 = [0.1231 0.0978 0.1008 0.0971 0.0974 0.1148 0.0971 0.1363 0.0896 0.0886 0.0924 0.0873 0.0884 0.09 0.0986 0.1518 0.0906 0.0913 0.1057];
MSE2 = [3.4428 0.1841 0.2279 0.1721 0.1857 0.2622 0.1793 0.2973 0.1084 0.141 0.1076 0.155 0.1401 0.1199 0.1378 0.1453 0.1142 0.0886 0.1571];
MSE_total = [MSE1;MSE2];

% Plot results
bar(n,MSE_total)
title('MSE resaults varying State Dimension')
xlabel('State Dimension')
ylabel('MSE')
legend('Dataset1','Dataset2')
grid on;
