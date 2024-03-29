clear all;
Is=0.01e-12;
Ib=0.1e-12;
Vb = 1.3;
Gp = 0.1;

V = linspace(-1.95,0.7,200);
I = zeros(size(V));

I = (Is .* (exp((1.2 / 0.025) .* V) - 1)) + (Gp .* V) - (Ib .* (exp((-1.2 / 0.025) .* (V + Vb)) - 1));
%I = (Is(exp(1.2*V/0.025))-1)+ Gp*V - (Ib(exp(-1.2*(V+Vb)/0.025)-1));

I_Noise = I + 0.2 .* I .* randn(1, 200);


polyFit_4th = polyfit(V, I, 4);
polyFit_8th = polyfit(V, I, 8);
polyFit_4thNoise = polyfit(V, I, 4);
polyFit_8thNoise = polyfit(V, I, 8);

figure(1);
subplot(2,1,1)
plot(V, I)
title('I-V Curve of a PN Diode With 4th Order Polynomial Fitting')
xlabel('V (V)')
ylabel('I (A)')
hold on
plot(V, polyval(polyFit_4th, V))
hold off

subplot(2,1,2)
semilogy(V, abs(I))
title('Logarithmic I-V Curve of a PN Diode With 4th Order Polynomial Fitting')
xlabel('V (V)')
ylabel('abs I (A)')
hold on
semilogy(V, abs(polyval(polyFit_4th, V)))
hold off

figure(2)
subplot(2,1,1)
plot(V, I_Noise)
title('I-V Curve of a PN Diode with Experimental Noise and With 4th Order Polynomial Fitting')
xlabel('V (V)')
ylabel('I (A)')
hold on
plot(V, polyval(polyFit_4thNoise, V))
hold off

subplot(2,1,2)
semilogy(V, abs(I_Noise))
title('Logarithmic I-V Curve of a PN Diode with Experimental Noise and With 4th Order Polynomial Fitting')
xlabel('V (V)')
ylabel('abs I (A)')
hold on
semilogy(V, abs(polyval(polyFit_4thNoise, V)))
hold off

figure(3)
subplot(2,1,1)
plot(V, I)
title('I-V Curve of a PN Diode with 8th Order Polynomial Fitting')
xlabel('V (V)')
ylabel('I  (A)')
hold on
plot(V, polyval(polyFit_8thNoise, V))
hold off

subplot(2,1,2)
semilogy(V, abs(I))
title('Logarithmic I-V Curve of a PN Diode with 8th Order Polynomial Fitting')
xlabel('V (V)')
ylabel('abs I (A)')
hold on
semilogy(V, abs(polyval(polyFit_8thNoise, V)))
hold off

figure(4)
subplot(2,1,1)
plot(V, I_Noise)
title('I-V Curve of a PN Diode with Experimental Noise and With 8th Order Polynomial Fitting')
xlabel('V (V)')
ylabel('I (A)')
hold on
plot(V, polyval(polyFit_8thNoise, V))
hold off

subplot(2,1,2)
semilogy(V, abs(I_Noise))
title('Logarithmic I-V Curve of a PN Diode with Experimental Noise and With 8th Order Polynomial Fitting')
xlabel('V(V)')
ylabel('abs I (A)')
hold on
semilogy(V, abs(polyval(polyFit_8thNoise, V)))
hold off

% Task #4: Nonlinear Curve Fitting

% Part A (Is and Ib)
foA = fittype('A .* (exp((1.2/0.025) .* x) - 1) + (0.1 .* x) - (C .* (exp((-1.2/0.025) .* (x + 1.3)) - 1))');
ffA = fit(V', I', foA);
IfA = ffA(V);
figure(5);
plot(V, IfA)
title('I-V Curve of a PN Diode with Fitted Parameters I_s and I_b')
xlabel('Voltage (V)')
ylabel('Current (A)')
hold on

% Part B (Is, Gb and Ib)
foB = fittype('A .* (exp((1.2/0.025) .* x) - 1) + (B .* x) - (C .* (exp((-1.2/0.025) .* (x + 1.3)) - 1))');
ffB = fit(V', I', foB);
IfB = ffB(V);
plot(V, IfB)
title('I-V Curve of a PN Diode with Fitted Parameters I_s, G_b and I_b')
xlabel('Voltage (V)')
ylabel('Current (A)')

% Part C (Is, Gb, Ib and Vb)
foC = fittype('A .* (exp((1.2/0.025) .* x) - 1) + (B .* x) - (C .* (exp((-1.2/0.025) .* (x + D)) - 1))');
ffC = fit(V', I', foC);
IfC = ffC(V);
plot(V, IfC)
title('I-V Curve of a PN Diode with Fitted Parameters I_s, G_b, I_b and V_b')
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('I_s and I_b', 'I_s, G_b and I_b', 'I_s, G_b, I_b and V_b')
hold off

% Task #5: Fitting using the Neural Net model

inputs = V.';
targets = I.';
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net, tr] = train(net, inputs, targets);
outputs = net(inputs);
errors = gsubtract(outputs, targets);
performance = perform(net, targets, outputs);
view(net)
Inn = outputs;

figure(6)
plot(V, I)
hold on
plot(V, Inn, 'r--')
hold off
legend("Raw Data", "Neural Net Fit")