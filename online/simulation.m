%% SISO

clc;
clear;
close all;

% the parameters I want to estimate
a = 2; b = 1;

% input
u = @(t) 5*sin(3*t);

% output
start = 0; step = 0.001; finish = 10; 
t = start:step:finish;
dx = @(t,x) -a*x + b*u(t);
[t, x] = ode45(dx,t,0);

% convert vector to function of time
x = @(t) x(round((t - start)/step + 1));

% plot
figure
hold on
tPlot = 0:0.01:min(10,finish);
plot(tPlot,x(tPlot))
plot(tPlot,u(tPlot))
xlabel('Time (s)'),ylabel('Amplitude'),legend('output','input')
print('in_out','-dpng','-r300')

%% GRADIENT DESCENT

% the root of the stable filter
eig = -1;

gammaRange = [1 2:4:24];
xHat = zeros(numel(t),numel(gammaRange));
aHat = zeros(numel(t),numel(gammaRange));
bHat = zeros(numel(t),numel(gammaRange));
counter = 1;
for i=gammaRange
    % learning rate
    gamma = i;

    [phi,aHatPerGamma,bHatPerGamma] = gradientDescent(x,u,t,eig,gamma);
    phi1 = phi(:,1); phi2 = phi(:,2);

    % convert vectors to function of time
    theta1 = @(t) -eig -aHatPerGamma(round((t - start)/step + 1));
    theta2 = @(t) bHatPerGamma(round((t - start)/step + 1));
    phi1 = @(t) phi1(round((t - start)/step + 1));
    phi2 = @(t) phi2(round((t - start)/step + 1));

    % estimate output
    xHatPerGamma = @(t) theta1(t).*phi1(t) + theta2(t).*phi2(t);
    
    xHat(:,counter) = xHatPerGamma(t);
    aHat(:,counter) = aHatPerGamma;
    bHat(:,counter) = bHatPerGamma;
    
    legendInfo{counter} = ['$\gamma$ = ' num2str(i)]; 
    counter = counter + 1;
end

% gamma = 14, index = 5, best choice
figure
plot(t,xHat(:,5),t,x(t)),legend('$\hat{x}$','x','interpreter','latex')
xlabel('Time (s)'),ylabel('Amplitude')
print('out_est','-dpng','-r300')

% plot
figure,plot(t,aHat),legend(legendInfo,'interpreter','latex'),xlabel('Time (s)')
title('Analysis $\hat{\alpha}$','interpreter','latex')
print('a_param','-dpng','-r300')
figure,plot(t,bHat),legend(legendInfo,'interpreter','latex'),xlabel('Time (s)')
title('Analysis $\hat{b}$','interpreter','latex')
print('b_param','-dpng','-r300')

%% LYAPUNOV PARALLEL

gammaRange = [1 10:10:20];

[xHat,aHat,bHat,legendInfo] = gammaAnalysis(@(x,u,t,gamma,thetaM)lyapunovParallel(x,u,t,gamma,thetaM),...
    gammaRange,x,u,t);

% plot
figure,plot(t,aHat),legend(legendInfo,'interpreter','latex')
title('Lyapunov parallel $\hat{\alpha}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_paral_a','-dpng','-r300')
figure,plot(t,bHat),legend(legendInfo,'interpreter','latex')
title('Lyapunov parallel $\hat{b}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_paral_b','-dpng','-r300')

% estimate output
figure,plot(t,xHat(:,1),t,x(t))
legend('$\hat{x}$', 'x', 'interpreter', 'latex')
xlabel('Time (s)')
print('lyap_parallel_x','-dpng','-r300')

%% LYAPUNOV MIX

[xHat,aHat,bHat,legendInfo] = gammaAnalysis(@(x,u,t,gamma,thetaM)lyapunovMix(x,u,t,gamma,thetaM),...
    gammaRange,x,u,t);

% plot
figure,plot(t,aHat),legend(legendInfo,'interpreter','latex')
title('Lyapunov mix $\hat{\alpha}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_mix_a','-dpng','-r300')
figure,plot(t,bHat),legend(legendInfo,'interpreter','latex')
title('Lyapunov mix $\hat{b}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_mix_b','-dpng','-r300')

% thetaM analysis
thetamRange = [1 2:2:10];
[xHat,aHat,bHat,legendInfo] = thetamAnalysis(thetamRange,x,u,t);

% plot
figure,plot(t,aHat),legend(legendInfo,'interpreter','latex')
title('Lyapunov mix $\hat{\alpha}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_mix_a_thetam','-dpng','-r300')
figure,plot(t,bHat),legend(legendInfo,'interpreter','latex')
title('Lyapunov mix $\hat{b}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_mix_b_thetam','-dpng','-r300')

% estimate output
figure,plot(t,xHat(:,1),t,x(t))
legend('$\hat{x}$', 'x', 'interpreter', 'latex')
xlabel('Time (s)')
print('lyap_mix_x','-dpng','-r300')

%% AMPLITUDE NOISE

amplitudeRange = 0.2:0.2:0.8;
% LYAPUNOV PARALLEL NOISE
[xHat,aHat,bHat,legendInfo] = amplitudeNoiseAnalysis(@(x,u,t,gamma,thetaM)lyapunovParallel(x,u,t,gamma,thetaM)...
    ,amplitudeRange,x,u,t);

% plot
figure,plot(t,aHat),legend(legendInfo,'interpreter','latex')
title('Noise Lyapunov parallel $\hat{\alpha}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_parallel_a_noise_ampl','-dpng','-r300')
figure,plot(t,bHat),legend(legendInfo,'interpreter','latex')
title('Noise Lyapunov parallel $\hat{b}$','interpreter','latex')
print('lyap_parallel_b_noise_ampl','-dpng','-r300')

% LYAPUNOV MIX NOISE
[xHat,aHat,bHat,legendInfo] = amplitudeNoiseAnalysis(@(x,u,t,gamma,thetaM)lyapunovMix(x,u,t,gamma,thetaM)...
    ,amplitudeRange,x,u,t);

% plot
figure,plot(t,aHat),legend(legendInfo,'interpreter','latex')
title('Noise Lyapunov mix $\hat{\alpha}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_mix_a_noise_ampl','-dpng','-r300')
figure,plot(t,bHat),legend(legendInfo,'interpreter','latex')
title('Noise Lyapunov mix $\hat{b}$','interpreter','latex')
print('lyap_mix_b_noise_ampl','-dpng','-r300')
xlabel('Time (s)')

%% FREQUENCY NOISE

freqRange = [10:10:50 1000];
% LYAPUNOV PARALLEL NOISE
[xHat,aHat,bHat,legendInfo] = freqNoiseAnalysis(@(x,u,t,gamma,thetaM)lyapunovParallel(x,u,t,gamma,thetaM)...
    ,freqRange,x,u,t);

% plot
figure,plot(t,aHat),legend(legendInfo,'interpreter','latex')
title('Noise Lyapunov parallel $\hat{\alpha}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_parallel_a_noise_freq','-dpng','-r300')
figure,plot(t,bHat),legend(legendInfo,'interpreter','latex')
title('Noise Lyapunov parallel $\hat{b}$','interpreter','latex')
print('lyap_parallel_b_noise_freq','-dpng','-r300')

% LYAPUNOV MIX NOISE
[xHat,aHat,bHat,legendInfo] = freqNoiseAnalysis(@(x,u,t,gamma,thetaM)lyapunovMix(x,u,t,gamma,thetaM)...
    ,freqRange,x,u,t);

% plot
figure,plot(t,aHat),legend(legendInfo,'interpreter','latex')
title('Noise Lyapunov mix $\hat{\alpha}$','interpreter','latex')
xlabel('Time (s)')
print('lyap_mix_a_noise_freq','-dpng','-r300')
figure,plot(t,bHat),legend(legendInfo,'interpreter','latex')
title('Noise Lyapunov mix $\hat{b}$','interpreter','latex')
print('lyap_mix_b_noise_freq','-dpng','-r300')
xlabel('Time (s)')

%% MIMO

% the parameters I want to estimate
a = [-0.25 3; -5 -1];
b = [1;2];

% input
u = @(t) 10*sin(2*t) + 5*sin(7.5*t);

% time
t = 0:0.001:30;

% learning rate
%gamma = [1 1];
gammaRange = [1 10];

% compare to the previous implementations, now we compute everything 
% (output,estimations) within one ode function in lyapunovParallelMIMO()
[x,xHat,aHat,bHat,legendInfo] = gammaAnalysisMIMO(gammaRange,u,t,a,b);

% plot
a = a';
for j=1:4
    figure
    hold on
    yline(a(j),'r--','HandleVisibility','off');
    for i=1:numel(gammaRange)
        plot(t,aHat(:,j,i))
    end
    legend(legendInfo,'interpreter','latex')
    xlabel('Time (s)')
    print(['lyap_parallel_a_mimo_' num2str(j)],'-dpng','-r300')
end

for j=1:2
    figure
    hold on
    yline(b(j),'r--','HandleVisibility','off');
    for i=1:numel(gammaRange)
        plot(t,bHat(:,j,i))
    end
    legend(legendInfo,'interpreter','latex')
    xlabel('Time (s)')
    print(['lyap_parallel_b_mimo_' num2str(j)],'-dpng','-r300')
end

figure,plot(t,xHat(:,1,1),t,x(:,1))
legend('$\hat{x_1}$', '$x_1$', 'interpreter', 'latex')
xlabel('Time (s)')
print('mimo_x1','-dpng','-r300')

figure,plot(t,xHat(:,2,1),t,x(:,2))
legend('$\hat{x_2}$', '$x_2$', 'interpreter', 'latex')
xlabel('Time (s)')
print('mimo_x2','-dpng','-r300')

%% FUNCTIONS

function [xHat,aHat,bHat,legendInfo] = gammaAnalysis(func,gammaRange,x,u,t) 
    xHat = zeros(numel(t),numel(gammaRange));
    aHat = zeros(numel(t),numel(gammaRange));
    bHat = zeros(numel(t),numel(gammaRange));
    counter = 1;
    for i=gammaRange
        % learning rate
        gamma = [i i];
        thetaM = 1;
        [xHatPerGamma,aHatPerGamma,bHatPerGamma] = func(x,u,t,gamma,thetaM);

        xHat(:,counter) = xHatPerGamma;
        aHat(:,counter) = aHatPerGamma;
        bHat(:,counter) = bHatPerGamma;

        legendInfo{counter} = ['$\gamma$ = ' num2str(i)]; 
        counter = counter + 1;
    end
end

function [xHat,aHat,bHat,legendInfo] = thetamAnalysis(thetamRange,x,u,t)
    xHat = zeros(numel(t),numel(thetamRange));
    aHat = zeros(numel(t),numel(thetamRange));
    bHat = zeros(numel(t),numel(thetamRange));
    counter = 1;
    for i=thetamRange
        % learning rate
        gamma = [1 1];
        [xHatPerGamma,aHatPerGamma,bHatPerGamma] = lyapunovMix(x,u,t,gamma,i);

        xHat(:,counter) = xHatPerGamma;
        aHat(:,counter) = aHatPerGamma;
        bHat(:,counter) = bHatPerGamma;

        legendInfo{counter} = ['$\hat{\theta_m}$ = ' num2str(i)]; 
        counter = counter + 1;
    end
end

function [xHat,aHat,bHat,legendInfo] = amplitudeNoiseAnalysis(func,amplitudeRange,x,u,t)
    xHat = zeros(numel(t),numel(amplitudeRange));
    aHat = zeros(numel(t),numel(amplitudeRange));
    bHat = zeros(numel(t),numel(amplitudeRange));
    counter = 1;
    for i=amplitudeRange
        amplitude = i;
        f = 20;
        noise = @(t) amplitude*sin(2*pi*f*t);
        x = @(t) x(t) + noise(t);

        gamma = [1 1];
        thetaM = 1;
        [xHatPerAmpl,aHatPerAmpl,bHatPerAmpl] = func(x,u,t,gamma,thetaM);
        
        xHat(:,counter) = xHatPerAmpl;
        aHat(:,counter) = aHatPerAmpl;
        bHat(:,counter) = bHatPerAmpl;

        legendInfo{counter} = ['$\eta_0$ = ' num2str(i)]; 
        counter = counter + 1;
    end
end

function [xHat,aHat,bHat,legendInfo] = freqNoiseAnalysis(func,amplitudeRange,x,u,t)
    xHat = zeros(numel(t),numel(amplitudeRange));
    aHat = zeros(numel(t),numel(amplitudeRange));
    bHat = zeros(numel(t),numel(amplitudeRange));
    counter = 1;
    for i=amplitudeRange
        amplitude = 0.15;
        f = i;
        noise = @(t) amplitude*sin(2*pi*f*t);
        x = @(t) x(t) + noise(t);

        gamma = [1 1];
        thetaM = 1;
        [xHatPerFreq,aHatPerFreq,bHatPerFreq] = func(x,u,t,gamma,thetaM);
        
        xHat(:,counter) = xHatPerFreq;
        aHat(:,counter) = aHatPerFreq;
        bHat(:,counter) = bHatPerFreq;

        legendInfo{counter} = ['f = ' num2str(i)]; 
        counter = counter + 1;
    end
end

function [x,xHat,aHat,bHat,legendInfo] = gammaAnalysisMIMO(gammaRange,u,t,a,b)
    xHat = zeros(numel(t),2,numel(gammaRange));
    aHat = zeros(numel(t),4,numel(gammaRange));
    bHat = zeros(numel(t),2,numel(gammaRange));
    counter = 1;
    for i=gammaRange
        % learning rate
        gamma = [i i];
        [x,xHatPerGamma,aHatPerGamma,bHatPerGamma] = lyapunovParallelMIMO(u,t,a,b,gamma);
    
        xHat(:,:,counter) = xHatPerGamma;
        aHat(:,:,counter) = aHatPerGamma;
        bHat(:,:,counter) = bHatPerGamma;

        legendInfo{counter} = ['$\gamma$ = ' num2str(i)]; 
        counter = counter + 1;
    end
end