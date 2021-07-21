%% stability eigenvalues for the soap film problem


function EV = catstability(b)
% global b; 
% global L;
% b = 0.8255;
L = 2*sqrt(1-b^2);

if abs(b - 0.866810925772831) < 0.000001
    lambda = -0.9;
else
    lambda = 1;  
end
solinit = bvpinit(linspace(-L/2,L/2,11),@(x)mat4init(x,L),lambda);
options = bvpset('RelTol',1e-6);
sol = bvp4c(@(x,y,lambda) mat4ode(x,y,lambda,b),@mat4bc,solinit,options);

% xvec = sol.x;
% yvec = sol.y(1,:);
% ypvec = sol.y(2,:);

EV = sol.parameters;
% 
% xint = linspace(-L/2,L/2);
% yxint = deval(sol,xint);
% plot(xint,yxint(1,:)); hold on;
% v = linspace(-L/2,L/2);
% plot(v,(1./(v.^2+b^2)-1)*40/46)



end




% ------------------------------------------------------------
function dydx = mat4ode(x,y,lambda,b)
% global b;
dydx = [  y(2)
         -x/(x^2+b^2)*y(2) - 2*b^2/(x^2+b^2)^2*y(1) + lambda*y(1) ];
end
% ------------------------------------------------------------
function res = mat4bc(ya,yb,lambda)
res = [  ya(1) 
         yb(1) 
        ya(2)-1/2 ];
end
% ------------------------------------------------------------
function yinit = mat4init(x,L)
% global L;
yinit = [  -L/pi*cos(x*pi/L)
           sin(x*pi/L) ];
end



