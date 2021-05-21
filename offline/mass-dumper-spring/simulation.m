%%% Theodoros Katzalis AEM:9282

clear;
clc;

m = 15;
b = 0.2;
k = 2;

u = @(t)5*sin(2*t) + 10.5;
timeRange = 0:0.1:10;

y0(1) = 0;
y0(2) = 0;
[time, y] = ode45(@(time, y)firstOrder(time, y, u, m, b, k), timeRange, y0);

polesRange = 0.1:0.1:10; % best pole around 0.5 based on plots
count = 0;
for i = polesRange
    count = count + 1;
    eig1 = -i;
    eig2 = -i;
    filterCoeff(:,count) = [-(eig1+eig2) eig1*eig2];
    thetaEstimate(:,count) = estimator(eig1, eig2, y, u, timeRange);

    theta(1,count) = thetaEstimate(1,count) + filterCoeff(1,count);
    theta(2,count) = thetaEstimate(2,count) + filterCoeff(2,count);
    theta(3,count) = thetaEstimate(3,count);

    mEst(count) = 1/theta(3,count);
    kEst(count) = theta(2,count) * mEst(count);
    bEst(count) = theta(1,count) * mEst(count);
end

figure();
plot(polesRange, mEst, polesRange, kEst, polesRange, bEst);
grid on;
legend('m', 'k', 'b');

figure();
title("Estimated parameters");
subplot(3,1,1);
plot(polesRange, mEst);
grid on;
legend("mEst");
subplot(3,1,2);
plot(polesRange, kEst, 'color', 'magenta');
grid on;
legend("kEst");
subplot(3,1,3);
plot(polesRange, bEst, 'color', 'red');
grid on;
legend("bEst");
xlabel('poles $\Lambda(s)$', 'interpreter','latex');
print('estimations', '-dpng', '-r300');

figure();
title("Estimation error");
subplot(3,1,1);
plot(polesRange, m - mEst);
grid on;
legend("Error m");
subplot(3,1,2);
plot(polesRange, k - kEst, 'color', 'magenta');
grid on;
legend("Error k");
subplot(3,1,3);
plot(polesRange, b - bEst, 'color', 'red');
grid on;
legend("Error b");
xlabel('poles $\Lambda(s)$', 'interpreter','latex');
print('errors', '-dpng', '-r300');

fprintf("Estimations:\nm = %f\nk = %f\nb = %f\n", mEst(1), kEst(1), bEst(1));

function theta = estimator(eig1, eig2, y, u, timeRange)

filter = [1 -(eig1+eig2) eig1*eig2];

phi = zeros(length(timeRange), 3);
phi(:,1) = lsim(tf(-[1 0], filter), y(:,1), timeRange);
phi(:,2) = lsim(tf(-[0 1], filter), y(:,1), timeRange);
phi(:,3) = lsim(tf([0 1], filter), u(timeRange), timeRange);

theta = y(:,1)' * phi / (phi' * phi);

end

function dy = firstOrder(t, y, u, m, b, k) 

dy(1) = y(2);
dy(2) = -(b/m)*y(2) -(k/m)*y(1) + u(t)/m; 

dy = dy';

end