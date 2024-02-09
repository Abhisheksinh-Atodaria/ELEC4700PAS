set(0,'DefaultFigureWindowStyle','docked')
% Parameters
nx = 100; % Number of points in x-direction
ny = 100; % Number of points in y-direction
ni = 1000; % Maximum number of iterations

% Initializing variables
V = zeros(nx, ny); % Initial guess (could be anything)
V(:,nx) = 1; % Left side Dirichlet BC
V(:,ny) = 0; % Right side Dirichlet BC

% Iterative solution
for k = 1:ni
    V_new = V;
    for i = 2:nx-1
        for j = 2:ny-1
            V_new(i,j) = 0.25 * (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1));
        end
    end
    V = V_new;
end

% Plotting the solution
if mod(k,50) == 0
    surf(V')
    pause(0.05)
end
% [X,Y] = meshgrid(1:ny, 1:nx);
% surf(X, Y, V);
% xlabel('y');
% ylabel('x');
% zlabel('Voltage (V)');
% title('2D Laplace Equation Solution');

% Boundary conditions reset for the second part
V = zeros(nx, ny);
V(:,[1,ny]) = 1; % Left and right side Dirichlet BC
V([1,nx],:) = 0; % Top and bottom Dirichlet BC

% Simulation using imboxfilt
V_smoothed = imboxfilt(V, 3);

% Plotting smoothed solution
figure;
surf(X, Y, V_smoothed);
xlabel('y');
ylabel('x');
zlabel('Voltage (V)');
title('Smoothed Solution using imboxfilt');