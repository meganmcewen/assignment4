%this code works through the transient simulation questions for the
%assignment. It models three different input signals: step, sine, Gaussian.

%begin with MNA setup for circuit
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
C(1,1) = C1;
C(1,2) = -C1;
C(2,1) = -C1;
C(2,2) = C1;
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

%now set up 3 input signals and plot results

%input 1: step function

start = 1; %starting point for timesteps
step = 1000; %total number of time steps
del = start/step; %time difference between start and stop

%set up stepping of input for V
V_step = zeros(7,1,step); %all inputs are initially 0, set up 3x3 to add in step

%set up stepping of input for F
F_step = zeros(7,1,step); 
%we want to simulate 1000 steps for 1s and change F_step to 1 at 0.03s, so
%.03*1000 = 30 (= first step we need to set to 1).
F_step(7,1,30:step) = 1;

%step through FD solution
for i = 2:1:step %need to start at 2 to get indices for V_step to work
    S = C/del + G;
    A = C*V_step(:,:,i-1)/del+F_step(:,:,i);
    V_step(:,:,i) = S\A; %set up V for each step
end

%extract Vin and Vout
Vo = V_step(5,1,:);
Vin = V_step(1,1,:);

%plot step function and Fourier
figure(1)
plot((1:step).*del, Vin(1,:))
hold on
plot((1:step).*del, Vo(1,:))
title('Step Function with transition to 1 starting at 0.03')
legend('Vin','Vo')

figure(2)
ff = abs(fftshift(fft(Vin(1,:))));
FF = abs(fftshift(fft(Vo(1,:))));
plot(((1:length(ff))/step)-0.5,ff)
hold on
plot(((1:length(FF))/step)-0.5,FF);
title('Step Function Fourier Transform')
legend('Vin','Vo');
xlim([-0.05 0.05]); %zoom in to area

%input 2: sine function
f = 1/0.03; %frequency in Hz; change this to see how freq. affects solution

%set up F and V matrices
F_sine = zeros(7,1,step);
V_sine = zeros(7,1,step);

%set up F with the sine function
for i = 1:1:step
    F_sine(7,1,i)= sin(2*pi*i*f*del);
end

%step through FD solution
for i = 2:1:step %need to start at 2 to get indices for V_step to work
    S = C/del + G;
    A = C*V_sine(:,:,i-1)/del+F_sine(:,:,i);
    V_sine(:,:,i) = S\A; %set up V for each step
end

%extract Vin and Vout
Vo = V_sine(5,1,:);
Vin = V_sine(1,1,:);

%plot sine function and Fourier
figure(3)
plot((1:step).*del, Vin(1,:))
hold on
plot((1:step).*del, Vo(1,:))
title('Sine function (f = 1/0.03 Hz)')
legend('Vin','Vo')

figure(4)
ff = abs(fftshift(fft(Vin(1,:))));
FF = abs(fftshift(fft(Vo(1,:))));
plot(((1:length(ff))/step)-0.5,ff)
hold on
plot(((1:length(FF))/step)-0.5,FF);
title('Sine Function Fourier Transform')
legend('Vin','Vo');
xlim([-0.15 0.15]); %zoom in to area

%input 3: Gaussian pulse
%set up constants for Gaussian
std = 0.03;
delay = 0.06;
mag = 1;

%set up inputs
F_gauss = zeros(7,1,step);
V_gauss = zeros(7,1,step);

%set up F with Gaussian
for i = 1:1:step
    F_gauss(7,1,i) = exp(-1*((del*i - delay)/std)^2);
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
figure(5)
plot((1:step).*del, Vin(1,:))
hold on
plot((1:step).*del, Vo(1,:))
title('Gauss Function (std = 0.05, mag = 1, delay = 0.06')
legend('Vin','Vo')

figure(6)
ff = abs(fftshift(fft(Vin(1,:))));
FF = abs(fftshift(fft(Vo(1,:))));
plot(((1:length(ff))/step)-0.5,ff)
hold on
plot(((1:length(FF))/step)-0.5,FF);
title('Gauss Function Fourier Transform')
legend('Vin','Vo');
xlim([-0.05 0.05]); %zoom in to area

