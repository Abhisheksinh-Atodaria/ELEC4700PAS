% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2 = 3.44;          % Core
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 10;         % number of modes to compute

[x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3],[h1,h2,h3], ...
                                            rh,rw,side,dx,dy); 

% First consider the fundamental TE mode:

% Loop through modes
for mode_index = 1:nmodes
    % Compute TE mode
    [Hx, Hy, neff] = wgmodes(lambda, n2, mode_index, dx, dy, eps, '000A');
    fprintf(1, 'TE Mode %d - neff = %.6f\n', mode_index, neff);
    
    % Plot TE mode
    figure(mode_index);
    subplot(2, 2, 1);
    contourmode(x, y, Hx(:,:,mode_index));
    title(['TE Mode Hx', num2str(mode_index)]); xlabel('x'); ylabel('y');
    for v = edges, line(v{:}); end
    
    subplot(2, 2, 2);
    contourmode(x, y, Hy(:,:,mode_index));
    title(['Hy TE Mode', num2str(mode_index)]); xlabel('x'); ylabel('y');
    for v = edges, line(v{:}); end

    % Compute TM mode
    [Hx, Hy, neff] = wgmodes(lambda, n2, mode_index, dx, dy, eps, '000S');
    fprintf(1, 'TM Mode %d - neff = %.6f\n', mode_index, neff);

    % Plot TM mode
    figure(mode_index);
    subplot(2, 2, 3);
    contourmode(x, y, Hx(:,:,mode_index));
    title(['TM Mode Hx ', num2str(mode_index)]); xlabel('x'); ylabel('y');
    for v = edges, line(v{:}); end

    subplot(2, 2, 4);
    contourmode(x, y, Hy(:,:,mode_index));
    title(['TM Mode Hy ', num2str(mode_index)]); xlabel('x'); ylabel('y');
    for v = edges, line(v{:}); end
end