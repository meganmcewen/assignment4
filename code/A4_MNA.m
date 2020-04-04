%this is the code used in the MNA-PA 
%it sets up the equations for the circuit and sweeps voltage

%set values
C1 = 0.25;
R1 = 1;
R2 = 2;
L = 0.2;
R3 = 10;
a = 100;
R4 = 0.1;
Ro = 1000;
omega = pi;
Vin = 10;

%set up matrices
C = zeros(7,7);
G = zeros(7,7);
V = zeros(7,1);
F = zeros(7,1);

%add values to matrices
C(1,1) = 0.25;
C(1,2) = -0.25;
C(2,1) = -0.25;
C(2,2) = 0.25;
C(6,7) = -L;

G(1,1) = 1/R1;
G(1,2) = -1/R1;
G(1,6) = 1;
G(2,1) = -1/R1;
G(2,2) = 1/R2 + 1/R1;
G(2,7) = 1;
G(3,3) = 1/R3;
G(3,7) = -1;
G(4,4) = 1;
G(4,7) = -a;
G(5,4) = -1/R4;
G(5,5) = 1/R4 + 1/Ro;
G(6,2) = 1;
G(6,3) = -1;
G(7,1) = 1;

F(7,1) = 1;

%do DC sweep
Vo = [];
V3 = [];
for i = -10:1:10
    F(7,1) = i;
    V = G\F;
    Vo = [Vo V(5)];
    V3 = [V3 V(3)];
end
   
figure(1)
plot(-10:1:10,V3);
hold on
plot(-10:1:10,Vo);
title('DC Sweep');
legend('V3','Vo');

f = logspace(0,5,5000);
VGain = [];

for i = 1:length(f)
    w = 2*pi*f(i);
    s = 1i*f(i);
    A = G + (s*C);
    V = A\F;
    Vo(i) = abs(V(5));
    VGain(i) = 20*log(abs(V(5)));
end

figure(2)
semilogx(f,Vo);
title('Output Voltage vs Frequency');


figure(3)
semilogx(f,VGain);
title('Gain vs Frequency');

%histogram with changing C
VCap = [];
F(7) = 1;
std = 0.05;
range = 10000;
rand_c = std.*randn(range,1) + 0.25;
for n=1:range
    j = randn()*0.05;
    C(1,1)= rand_c(n);
    C(2,2)= rand_c(n);
    C(1,2)= -rand_c(n);
    C(2,1)= -rand_c(n);
    V = (G+(1j*pi*pi*2*C))\F;
    VCap = [VCap 20*log10(abs(V(5)/F(7)))];
end

figure(4)
histogram(VCap);
xlabel('Gain')
ylabel('Count')
title('Histogram of Gain Across Cap')