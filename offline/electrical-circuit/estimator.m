function theta = estimator(eig1, eig2, N, VC, input1, input2, timeRange)

filter = [1 -(eig1+eig2) eig1*eig2];

phi = zeros(N, 6);
phi(:,1) = lsim(tf(-[1 0], filter), VC, timeRange);
phi(:,2) = lsim(tf(-[0 1], filter), VC, timeRange);
phi(:,3) = lsim(tf([1 0],filter), input1(timeRange), timeRange);
phi(:,4) = lsim(tf([0 1],filter), input1(timeRange), timeRange);
phi(:,5) = lsim(tf([1 0],filter), input2, timeRange);
phi(:,6) = lsim(tf([0 1],filter), input2, timeRange); 

theta = VC * phi / (phi.' * phi);

end