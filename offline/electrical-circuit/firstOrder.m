function dx = firstOrder(t,x,invRC,invLC)

% diff(Vc,t,2) + diff(Vc,t,1)*1/RC + Vc/LC = diff(u1,t,1)/RC + diff(u2,t,1)+ u2/LC
% u1(t) = 2sin(t)
% u2(t) = 1

dx(1) = x(2);
dx(2) = -x(2)*invRC - x(1)*invLC + 2*cos(t)*invRC + 0 + invLC;

dx = dx';

end