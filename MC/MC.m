set(0, 'DefaultFigureWindowStyle', 'docked')
C.m_0 = 9.10938215e-31;            
n = 200;
dt = 1;
F = 1e-20;
t = 1 : n;
Ps = 0.05;
numParticles = 10;

x = zeros(numParticles, n);
v = zeros(numParticles, n);
v(1) = 0;
x(1) = 0;

for i = 2 : n
    for j = 1 : numParticles
        P = rand();
        if P <= Ps
            v(j, i) = 0;
        else  
            v(j, i) = v(j, i - 1) + F * dt / C.m_0;
        end

        x(j, i) = x(j, i - 1) + v(j, i) * dt;

        driftVelocity = mean(v);
    
        figure(1) 
        subplot(2, 1, 1)
        plot(t, x)
        xlabel('Time (s)')
        ylabel('Postion')
    
        subplot(2, 1, 2)
        plot(t, v)
        xlabel('Time (s)')
        ylabel('Velocity')
        %title(['Drift Velocity: ', num2str(driftVelocity,'%0.5g')])
    end
end

% figure(2)
% plot(t, x)
% 
% figure(3)
% plot(t, v)