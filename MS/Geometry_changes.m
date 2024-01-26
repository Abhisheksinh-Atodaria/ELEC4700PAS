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
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 1;        % number of modes to compute

% Loop through Ridge half-width values
rhw = linspace(0.325, 1.0, 10);

n_e=[];
    
   
    % Loop through modes
    i=1;
    r= (1-0.325)/10
    for index = 0.325:r:1

         % Generate the waveguide mesh
    [x, y, xc, yc, nx, ny, eps, edges] = waveguidemesh([n1, n2, n3], [h1, h2, h3], rh, index, side, dx, dy);
    
        % Compute TE mode
        [Hx, Hy, neff] = wgmodes(lambda, n2,  nmodes, dx, dy, eps, '000A');
        fprintf(1, ' neff = %.6f\n',  index, neff);
        
        % Plot TE mode
        figure(i);
        subplot(2, 1, 1);
        contourmode(x, y, Hx);
        title(['Ridge width TE Mode Hx', num2str(i)]); xlabel('x'); ylabel('y');
        for v = edges, line(v{:}); end

        subplot(2, 1, 2);
        contourmode(x, y, Hy);
        title(['Ridge width Hy TE Mode' num2str(i)]); xlabel('x'); ylabel('y');
        for v = edges, line(v{:}); end
        i=i+1;
        n_e =[n_e,neff];
    end
        figure(11);
        subplot(1, 1, 1);
        plot(n_e)
        contourmode(x, y);
        title(['N-effective' ]); xlabel('x'); ylabel('y');
        for v = edges, line(v{:}); end
        % % Compute TM mode
        % [Hx, Hy, neff] = wgmodes(lambda, n2, i, dx, dy, eps, '000S');
        % fprintf(1, 'Ridge width %.3f - TM Mode %d - neff = %.6f\n', rw, i, neff);

        % % Plot TM mode
        % figure(i);
        % subplot(2, 2, 3);
        % contourmode(x, y, Hx(:,:,i));
        % title(['Ridge width ', num2str(rw), ' TM Mode Hx ', num2str(i)]); xlabel('x'); ylabel('y');
        % for v = edges, line(v{:}); end
        % 
        % subplot(2, 2, 4);
        % contourmode(x, y, Hy(:,:,i));
        % title(['Ridge width ', num2str(rw), ' TM Mode Hy ', num2str(i)]); xlabel('x'); ylabel('y');
        % for v = edges, line(v{:}); end