G = sparse(7,7);
G1 = 1;
G2 = 1/2;
G3 = 1/10;
G4 = 1/0.1;
G5 = 1/1000;

C1 = 0.25;
L = 0.2;
a = 100;
% X = [V1,V2,Is,Il,V3,V4,V5];
w= 1;

Vin = 1;

G = [1,0,0,0,0,0,0;
     G1,-G1,1,0,0,0,0;
     -G1,G1+G2,0,1,0,0,0;
     0,0,0,-1,G3,0,0;
     0,0,0,0,-a*G3,1,0;
     0,1,0,0,-1,0,0;
     0,0,0,0,0,-G4,G4+G5];
C =[0,0,0,0,0,0,0;
    C1,-C1,0,0,0,0,0;
    -C1,C1,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,-L,0,0,0;
    0,0,0,0,0,0,0];
F = [Vin;
    0;
    0;
    0;
    0;
    0;
    0];

% C*diff(X) + G*X = F;

VIN = linspace(-10,10,100);

for i = 1:100
    F(1) = VIN(i);
    Vop = G\F;
    % w = W(i);
    % Vop = (G + 1i*w*C)\F;
    V5(i) = Vop(7);
    V3(i) = Vop(5);

end
figure
    plot(VIN,V5)
    plot(VIN,V3)
% for i = 1:100
%     F(1) = VIN(i);
%     Vop = G\F;
%     w = w +1;
%     Vop = (G + 1i*w*C)\F;
%     V5(i) = Vop(7);
%     V3(i) = Vop(5);
%      imatrix(i) = w;
% end
    %figure
    %plot(imatrix,real(V5),'r'); hold on
    %plot(imatrix,V3,'b');
    %hold off
