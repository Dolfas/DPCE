clc, clear all
%% Compute the unconstrained minimum 

x0=[-1; 1];  %Intial estimate
options = optimoptions('fminunc','Algorithm','quasi-newton');
xopt=fminunc(@Rosenbrock_function,x0,options)

%% Compute the constrained minimum

x0=[-1; 1];  %Intial estimate
A = [1, 0];  %Constrain
B = 0.5; 
xoptconstr=fmincon(@Rosenbrock_function,x0,A,B)

%% Plots 
% Range of independent variables to consider 
x1min=-2;
x1max=2;
x2min=-2;
x2max=2;

% Number of intervals in the mesh grid
N1=100;
N2=100;

xv1 = linspace(x1min,x1max,N1);
xv2 = linspace(x2min,x2max,N2);
[xx1,xx2] = meshgrid(xv1,xv2);

% Computes the function at the different points of the mesh grid
for ii=1:N1
    for jj=1:N2
        x=[xx1(ii,jj); xx2(ii,jj)];
        ff(ii,jj)=Rosenbrock_function(x);
    end
end

% Plots the level curves using the Matlab function contoutr
Nlevel=10;  % Number of level curves in the contour plot
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
figure(1), contour(xv1,xv2,ff,Nlevel,LW,1.2), colorbar
axis([x1min x1max x2min x2max]), axis square
hold on
% Plots the initial point as a red circle
gg=plot(x0(1),x0(2),'or');
set(gg,'Linewidth',1.5);

% Plots the final estimate of the unconstrained minimum as a red cross
gg=plot(xopt(1),xopt(2),'xr');
set(gg,'Linewidth',1.5);

% Plots the final estimate of the constrained minimum as a red star
gg=plot(xoptconstr(1),xoptconstr(2),'*r');
set(gg,'Linewidth',1.5);

%plots the constraint boundary

% z1c=[-2:0.1:0.5];
% z2c=[2:-0.1:0.5];
gg = xline(0.5,'k');
set(gg,'Linewidth',1.5);

% Identifies axis
gg=xlabel('x_1');
set(gg,'FontSize',14);

gg=ylabel('x_2');
set(gg,'FontSize',14);

gg = title('Level lines of the Rosenbrock function');
set(gg,'FontSize',14);

gg = legend('Level lines','Initial Estimate', 'Unconstrained Minimum', 'Constrained Minimum');


hold off

% Plot the function in 3D 
figure(2)
surf(xx1,xx2,ff);

% Identifies axis
gg=xlabel('x_1');
set(gg,'FontSize',14);

gg=ylabel('x_2');
set(gg,'FontSize',14);

gg=zlabel('f(x)');
set(gg,'FontSize',14);

gg = title('Rosenbrock Function');
set(gg,'FontSize',14);