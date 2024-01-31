set(0,'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultLineLineWidth',2);
set(0,'Defaultaxeslinewidth',2)

set(0,'DefaultFigureWindowStyle','docked')

%defining constants
c_c = 299792458;
c_eps_0 = 808542149e-12;
c_eps_0_cm = c_eps_0/100;
c_mu_0 = 1/c_eps_0/c_c^2;
c_q = 1.60217653e-19;
c_hb = 1.05457266913e-34;
c_h = c_hb*2*pi;

RL = 0.9i;
RR = 0.9i;

%Input parameters 
InputParasL.E0 = 1e12;
InputParasL.we = 4e12;
InputParasL.t0 = 2e-12;
InputParasL.wg = 100;
InputParasL.phi = 0;
InputParasR = 10;
%InputParasR.E0 = 1e5;
%InputParasR.we = 10e10;
%InputParasR.t0 = 2e-12;
%InputParasR.wg = 10e-10;
%InputParasR.phi = 0;

%defining the refractive index and wavelength
n_g = 3.5;
vg = c_c/n_g*1e2;
Lambda = 1550e-9;

%simulation setup
plotN = 10;

%boundry functions
L = 1000e-6*1e2;
XL = [0,L];
YL = [-InputParasL.E0,InputParasL.E0];


Nz = 100;
dz = L/(Nz-1);
dt = dz/vg;
fsync = dt*vg/dz;

Nt = floor(5*Nz);
tmax = Nt*dt;
t_L = dt*Nz;

z = linspace(0,L,Nz).';
time = nan(1,Nt);
InputL = nan(1,Nt);
InputR = nan(1,Nt);
OutputL = nan(1,Nt);
OutputR = nan(1,Nt);

%Initializing electric field vectors
Ef =zeros(size(z));
Er = zeros(size(z));

%defining the electric field vectors
Ef1 = @SourceFct;
ErN = @SourceFct;

%Input functions
t = 0;
time(1) = t;
InputL(1) = Ef1(t,InputParasL);
InputR(1) = ErN(t,InputParasR);

%output functions

OutputR(1) = Er(1);
OutputL(1) = Ef(Nz);


%plotting the vectors
Ef(1) = InputL(1);
Er(Nz) = InputR(1);

figure('name','Fields')
subplot(3,1,1)
plot(z*10000,real(Ef),'r');
hold off
xlabel('z(\mum)')
ylabel('E_f')
subplot(3,1,2)
plot(z*10000,real(Er),'b');
xlabel('z(\mum)')
ylabel('E_r')
hold off
subplot(3,1,3)
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR),'r--');
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputL),'b--');
xlabel('time(ps)')
ylabel('E')

hold off

for i = 2:Nt
    t = dt*(i-1);
    time(i) = t;

    InputL(i) = Ef1(t,InputParasL);%setting input equal to source fucntion
    InputR(i) = ErN(t,InputParasR);%setting input equal to source fucntion

    
    Ef(1) = InputL(i) + RL*Er(1);%setting The electric field vectors equal to the Inputs with the reflectio added
    Er(Nz) = InputR(i) + RR*Ef(Nz);%setting The electric field vectors equal to the Inputs with the reflectio added
    
    

    Ef(2:Nz) = fsync*Ef(1:Nz-1);
    Er(1:Nz-1) = fsync*Er(2:Nz);

    OutputR(i) = Ef(Nz)*(1-RR);
    OutputL(i) = Er(1)*(1-RL);

    %plotting Ef
    if mod(i,plotN) == 0
        subplot(3,1,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off
        subplot(3,1,2)
        plot(z*10000,real(Er),'b'); hold on
        plot(z*10000,imag(Er),'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re','\Im')
        
        hold off
        %plotting Er
        subplot(3,1,3)
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12,real(OutputR),'g');
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12,real(OutputL),'m');
        xlim([0,Nt*dt*1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output' ...
            ,'Location','east')
        hold off
        pause(0.01)
    end
end