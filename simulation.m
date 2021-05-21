clc;
clear;
close all;

% the parameters I want to estimate
a = 2; b = 1;

% input
u = @(t) 5*sin(3*t);

% output
start = 0; step = 0.01; 
t = start:step:100;
dx = @(t,x) -a*x + b*u(t);
[t, x] = ode45(dx,t,0);

% convert vector to function of time
x = @(t) x(round((t - start)/step + 1));

% plot
%figure,hold on
%plot(t,x(t)),plot(t,u(t)),legend('output','input')

% gradient descent
[aHat,bHat] = gradientDescent(x,u,t);

% plot
figure
hold on
plot(t,aHat),plot(t,bHat),legend('a','b')
title('Gradient descent')

% lyapunov method Parallel
[xHat,aHat,bHat] = lyapunovParallel(x,u,t);

% plot
figure
subplot(211)
hold on
plot(t,aHat),plot(t,bHat),legend('a','b')
title('Lyapunov Parallel')

% lyapunov method Mix
[xHat,aHat,bHat] = lyapunovMix(x,u,t);

% plot
subplot(212)
hold on
plot(t,aHat),plot(t,bHat),legend('a','b')
title('Lyapunov Mix')

% NOISE
amplitude = 1;
noise = @(t) 0.15*sin(2*pi*20*t);
x = @(t) x(t) + noise(t);

% lyapunov method Parallel
[xHat,aHat,bHat] = lyapunovParallel(x,u,t);

% plot
figure
subplot(211)
hold on
plot(t,aHat),plot(t,bHat),legend('a','b')
title('Lyapunov Parallel with noise')

% lyapunov method Mix
[xHat,aHat,bHat] = lyapunovMix(x,u,t);

% plot
subplot(212)
hold on
plot(t,aHat),plot(t,bHat),legend('a','b')
title('Lyapunov Mix with noise')

% MIMO


