function [phi,aHat,bHat] = gradientDescent(x,u,t,eig,gamma)
[t, y] = ode45(@(t,y)diffSystem(t,y,u,x,gamma,eig),t,[0 0 0 0]);
phi = [y(:,1) y(:,2)];
aHat = -eig -y(:,3);
bHat = y(:,4);
end

function dy = diffSystem(t,y,u,x,gamma,eig)
phi1 = y(1); phi2 = y(2); theta(1) = y(3); theta(2) = y(4);

% thetaDot = gamma*error*phi
dPhi1 = eig*phi1 + x(t);
dPhi2 = eig*phi2 + u(t);
error = x(t) - theta(1)*phi1 - theta(2)*phi2;
dtheta(1) = gamma*error*phi1;
dtheta(2) = gamma*error*phi2;

dy = [dPhi1; dPhi2; dtheta(1); dtheta(2)];
end