set(0,'DefaultFigureWindowStyle','docked')
% Parameters
nx = 100; % Number of points in x-direction
ny = 100; % Number of points in y-direction
ni = 10000; % Maximum number of iterations

% Initializing variables
V = zeros(nx, ny); % Initial guess
% V(nx,:) = 1; 
% V(1,:) = 0;
% V(:,1) = 0;
% V(:,ny) = 0;
V_new2 = V;

% Iterative solution
for k = 1:ni
    V_new = V;
    for i = 1:nx
        for j = 1:ny

            if i == 1 
                V_new(i,j) = 1;
                 V_new2(i,j) = 1;
            elseif i == nx
                V_new(i,j) = 0;
                 V_new2(i,j) = 1;

            elseif j == 1 
                %V_new(i,j) = 0;
                V_new(i,j) = V(i,j+1);
            elseif j == ny
                %V_new(i,j) = 0;
                V_new(i,j) = V(i,j-1);
            
            else
                 V_new(i,j) = 0.25 * (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1));
             
            end
            %V(:,ny) = V(:,1);
            % V(nx,:) = 1;
            % V(1,:) = 0;
        end
       
    end
     if mod(k,50) == 0
        
        subplot(2,1,1)
        surf(V_new')
        pause(0.05)
        
        V_new2 = imboxfilt(V, 3);
       subplot(2,1,2)
       surf(V_new2')
        
       
       % V(:,1) = 0;
       % V(:,ny) = 0;
     end
        
    V = V_new;
end


figure

[Ex,Ey] = gradient(V);

figure
quiver(-Ex',-Ey',0.8)













% [X,Y] = meshgrid(1:ny, 1:nx);
% surf(X, Y, V);
% xlabel('x');
% ylabel('y');
% zlabel('Voltage (V)');
% title('2D Laplace Equation Solution');

% Boundary conditions reset for the second part
% V = ones(nx, ny);
% V(:,[1,ny]) = 1; % Left and right side Dirichlet BC
% V([1,nx],:) = 0; % Top and bottom Dirichlet BC
% 
% % Simulation using imboxfilt
% V_smoothed = imboxfilt(V, 3);

% % Plotting smoothed solution
% figure;
% surf(X, Y, V_smoothed);
% xlabel('y');
% ylabel('x');
% zlabel('Voltage (V)');
% title('Smoothed Solution using imboxfilt');