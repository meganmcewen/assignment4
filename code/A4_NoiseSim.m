%this code modifies the circuit to include noise
%it includes two models for the added current source: random numbers and
%Gaussian

%begin with MNA setup for circuit
%set values
R1 = 1;
R2 = 2;
L = 0.2;
R3 = 10;
a = 100;
R4 = 0.1;
Ro = 1000;
omega = pi;
Vin = 10;
C1 = 0.25;
C_new = 0.0001; %change this value to see how changing Cn affects bandwidth
In = 0.001;

%set up matrices
C = zeros(7,7);
G = zeros(7,7);
V = zeros(7,1);
F = zeros(7,1);

%add values to matrices
C(1,1) = C1;
C(1,2) = -C1;
C(2,1) = -C1;
C(2,2) = C1;
C(3,3) = C_new;
C(5,7) = -L;


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

start = 1; %starting point for timesteps
step = 1000; %total number of time steps
del = start/step; %time difference between start and stop

%input 3: Gaussian pulse
%set up constants for Gaussian
std = 0.03;
delay = 0.06;
mag = 1;

%set up inputs
F_gauss = zeros(7,1,step);
V_gauss = zeros(7,1,step);

%set up F with Gaussian and random noise
for i = 1:1:step
    F_gauss(3,1,i) = In*randn;
    F_gauss(7,1,i) = exp(-1*((del*i - delay)/std)^2); %sine population
end

%step through FD solution
for i = 2:1:step %need to start at 2 to get indices for V_step to work
    S = C/del + G;
    A = C*V_gauss(:,:,i-1)/del+F_gauss(:,:,i);
    V_gauss(:,:,i) = S\A; %set up V for each step
end

%extract Vin and Vout
Vo = V_gauss(5,1,:);
Vin = V_gauss(1,1,:);

%plot Gauss function and Fourier
figure(1)
plot((1:step).*del, Vin(1,:))
hold on
plot((1:step).*del, Vo(1,:))
title(['Gauss Function (Cn = ',num2str(C_new),', timestep = ',num2str(del),'.)']);
legend('Vin','Vo')

figure(2)
ff = abs(fftshift(fft(Vin(1,:))));
FF = abs(fftshift(fft(Vo(1,:))));
plot(((1:length(ff))/step)-0.5,ff)
hold on
plot(((1:length(FF))/step)-0.5,FF);
title(['Gauss Function with Fourier and Noise (Cn = ',num2str(C_new),', timestep = ',num2str(del),'.)']);
legend('Vin','Vo');
xlim([-0.03 0.03]); %zoom in to area

