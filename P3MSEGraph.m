n=[1:19];
MSE1 = [0.1231 0.0978 0.1008 0.0971 0.0974 0.1148 0.0971 0.1363 0.0896 0.0886 0.0924 0.0873 0.0884 0.09 0.0986 0.1518 0.0906 0.0913 0.1057];
MSE2 = [3.4428 0.1841 0.2279 0.1721 0.1857 0.2622 0.1793 0.2973 0.1084 0.141 0.1076 0.155 0.1401 0.1199 0.1378 0.1453 0.1142 0.0886 0.1571];
MSE_total = [MSE1;MSE2];


% Plot results
b = bar(n,MSE_total)

% Log scale
%set(gca,'YScale','log');

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,2));
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontWeight','bold')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(round(b(2).YData,2));
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontWeight', 'bold')


title('MSE resaults varying State Dimension', 'FontWeight', 'bold')
xlabel('State Dimension', 'FontWeight', 'bold')
ylabel('MSE', 'FontWeight', 'bold')
lgd = legend('Dataset1','Dataset2')
set(lgd,'FontSize',18);
set(lgd,'FontWeight','bold');
grid on;
