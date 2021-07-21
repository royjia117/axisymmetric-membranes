%% soap film problem with fixed area instead of tension
%% E = mean curvature + area constraint
%% add Lagrange multiplier eta for zs + sin(psi) = 0 constraint
%% combines previous versions into one code

function [] = soapfilm_fixed_area_master()
clc; close all; tic;
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; %surface area
global X; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess; %stores guess from previous iteration
global N; %number of gridpoints
 
% kbarsolutions(); %explores nonzero kappa bar, positive tension (BL)
% kbarsolutions2(); %explores nonzero kappa bar, negative tension (WKB)
% thinsolutions(); %makes figures for A = 2*pi case
% thinsolutions2(); %makes figures for A = 2*pi*1.1 case
% thinsolutions3();  %makes figures for A = 2*pi*1.3 case
% catsolutions2_kbar2()
% thinsolutions2_blah(); %makes figures for A = 2*1.1*pi case
% thinsolutions2_kbar2(); %makes figures for A = 2*1.1*pi, kbar > 0 case
% thinsolutions2_kbar3(); %makes figures for A = 2*1.1*pi, kbar > 0 case
thinsolutions2_kbar4(); %makes figures for A = 2*1.1*pi, kbar > 0 case
% thinsolutions2_kbar5(); %makes figures for A = 2*1.1*pi, kbar > 0 case
% torussolutions(); %parameters that make the Willmore torus
%spheresolutions(); %parameters that make a sphere
% tethertransition(); %for A slightly larger than Ac, no maximal extension
% tethertransition2(); %for A = 10
% tethertransition3(); %for various areas
% tethertransition4_mode1(); %for area A = 2*pi*1.3
% tethertransition4_mode2(); %for area A = 2*pi*1.3
% tethertransition4_mode3(); %for area A = 2*pi*1.3
%  tethertransition4_mode4(); %for area A = 2*pi*1.3
% tethertransition4_mode5(); %for area A = 2*pi*1.3
%variablearea(); %remove area constraint, add stretching modulus
%higherorder(); %explore higher order buckling modes
%otherareas1(); %makes figures for other A = 2*pi*(1/3)
%otherareas2();%makes figures for other A = 2*pi*(2/3)
%otherareas3();%makes figures for other A = 2*pi*1.2
% kbarcomparison();
% kbarcomparison2();
% asymmetricshapes(); %A = 2*1.1*pi
% asymmetricshapes2(); 
% fvhfig();
% scsolutions(); %nonzero spontaneous curvature

% summaryfigure(); % plot F vs z vs A

% evmatrix()

% onecatfig()
% twocatfig()
% nocatfig()

% tocfig()

% close all;
% mergeshapefigs()
% mergeshapefigs2thick() %<----- set colorbar first
%  mergeshapefigs2thin() %<----- set colorbar first
% mergeshapefigs3()
% cbfig();

% coverfig()

toc;

end

function [] = coverfig()

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 8 4];

%% Willmore torus
subplot(2,2,4)
global kappa; kappa = 1;
global kbar; kbar = 0.1;
global A; A = 10.2591; %10.1439; %2*pi*sqrt(2)*2*pi*2-10.1439; %10.1439; %surface area
global X; X = 0.498454; %1/2;
global y_guess;
global N; N = 101; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 2 0.1 0 pi 0 2 0]); %inner solution
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 1.1 0.1 0 pi 0 2*pi*sqrt(2)*.75 0]); %outer solution
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 

%% calculate shape for various widths
Xvec = 0.498454*2;%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
% figure; hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r_w = y(3,:); z_w = y(4,:); psi_s_w = y(2,:); psi_w = y(1,:);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
%     hold on;  grid on;
%     plot(z,r,'LineWidth',2);
%     figure();
%     plot(t,psi_s + sin(psi)./r)
%     snapshotalt(1.5,1.35,psi,psi_s,z,r,X,kappa,N);
%     camzoom(2)
%     wt = centercrop(wt,0.2,0.2); imwrite(wt,'willmore_torus.tif','tif','resolution',600,'compression','none');
    y_guess = y;
end
 
%% tether
%% prepare a solution as initial guess (thin cat to tether)
A = 2*pi*1.2; kbar = 0;
X = 1.5; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi 63.5 2 -71/2/pi]);
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec4 = [1.5 1.35];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r_t = y(3,:); z_t = y(4,:); psi_s_t = y(2,:); psi_t = y(1,:); L = y(8,1);
%     muvec4(j) = y(7,1);
%     Fvec4(j) = -2*pi*y(9,1);
%     Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
%     hold on;  grid on;
    if mod(j,10) ==1
%          plot(z,r,'color',col1,'LineWidth',2);
    end
    if abs(X - 1.3/2) < 0.0001
%         plot(z,r,'color',[col1 6/7],'LineWidth',2);
%        shape1m = snapshot(0.7,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape1m = centercrop(shape1m,0.2,0.2); imwrite(shape1m,'shape1mt.png');
    end

   y_guess = y;
end

%% prepare a solution as initial guess
A = 2*pi*1.2;
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.16 0.1 0 pi 0 2 0]); %m = 5b
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec9 = [1.05:-0.01:0]/2;
muvec9 = zeros(1,length(Xvec9));
Fvec9 = zeros(1,length(Xvec9));
Evec9 = zeros(1,length(Xvec9));
for j = 1:length(Xvec9)
    X = Xvec9(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec9(j) = y(7,1);
    Fvec9(j) = -2*pi*y(9,1);
    Evec9(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col5,'LineWidth',2);
    end
    if 2*X == 0.7
        plot(z,r,'LineWidth',2);
        z_b = z; r_b = r; 
        psi_b = psi; psi_s_b = psi_s; break;
    end

   y_guess = y;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- make figure ----------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% imshow('bg3.jpg','parent',gca); hold on;
% a1=axes('position',[0.3 0.3 .5 .5])


b = 0.54716;
z = linspace(-0.6627,0.6627);
% z = linspace(-0.527,0.527);
% b = 0.8255;
r = b*cosh(z/b);
z = z - z(1);
psi = linspace(0,2*pi);
H_b = psi_s_b + sin(psi_b)./r_b; H_t = psi_s_t + sin(psi_t)./r_t; H_w = psi_s_w + sin(psi_w)./r_w;
z_b = z_b(end:-1:1); z_b = z_b - z_b(1);
z_t = z_t(end:-1:1); z_t = z_t - z_t(1);
z_w = z_w(end:-1:1); z_w = z_w - z_w(1);
eps = 0.025; dd = 3;
h = surf(r_b'*cos(psi),r_b'*sin(psi)+dd,repmat(z_b,100,1)',repmat(H_b,100,1)'); hold on;
lh = surf((1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)+dd,repmat(eps*sin(psi),100,1)'-z_b(1),100*ones(size(psi'*psi)));
rh = surf((1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)+dd,repmat(eps*sin(psi),100,1)'+z_b(end),100*ones(size(psi'*psi)));
% h = surf(repmat(z_b,100,1)',r_b'*cos(psi),r_b'*sin(psi)+dd,repmat(H_b,100,1)'); hold on;
% lh = surf(repmat(eps*sin(psi),100,1)'-z_b(1),(1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)+dd,100*ones(size(psi'*psi)));
% rh = surf(repmat(eps*sin(psi),100,1)'-z_b(end),(1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)+dd,100*ones(size(psi'*psi)));
shading interp
h2 = surf(r'*cos(psi),r'*sin(psi),repmat(z,length(psi),1)',zeros(size(r'*psi))); hold on;
lh2 = surf((1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi),repmat(eps*sin(psi),100,1)'-z(1),100*ones(size(psi'*psi)));
rh2 = surf((1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi),repmat(eps*sin(psi),100,1)'+z(end),100*ones(size(psi'*psi)));
% h2 = surf(repmat(z,length(psi),1)',r'*cos(psi),r'*sin(psi),zeros(size(r'*psi))); hold on;
% lh2 = surf(repmat(eps*sin(psi),100,1)'-z(1),(1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi),100*ones(size(psi'*psi)));
% rh2 = surf(repmat(eps*sin(psi),100,1)'-z(end),(1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi),100*ones(size(psi'*psi)));
shading interp
h3 = surf(r_w'*cos(psi),r_w'*sin(psi)-dd,repmat(z_w,length(psi),1)',repmat(H_w,100,1)'); hold on;
lh3 = surf((1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)-dd,repmat(eps*sin(psi),100,1)'-z_w(1),100*ones(size(psi'*psi)));
rh3 = surf((1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)-dd,repmat(eps*sin(psi),100,1)'+z_w(end),100*ones(size(psi'*psi)));
shading interp
h4 = surf(r_t'*cos(psi),r_t'*sin(psi)-2*dd,repmat(z_t,length(psi),1)',repmat(H_t,100,1)'); hold on;
lh4 = surf((1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)-2*dd,repmat(eps*sin(psi),100,1)'-z_t(1),100*ones(size(psi'*psi)));
rh4 = surf((1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)-2*dd,repmat(eps*sin(psi),100,1)'+z_t(end),100*ones(size(psi'*psi)));
% h3 = surf(repmat(z_t,length(psi),1)',r_t'*cos(psi),r_t'*sin(psi)-dd,repmat(H_t,100,1)'); hold on;
% lh3 = surf(repmat(eps*sin(psi),100,1)'-z_t(1),(1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)-dd,100*ones(size(psi'*psi)));
% rh3 = surf(repmat(eps*sin(psi),100,1)'-z_t(end),(1+eps*cos(psi))'*cos(psi),(1+eps*cos(psi))'*sin(psi)-dd,100*ones(size(psi'*psi)));
shading interp



axis equal;


% h = ezsurf('sin(sqrt(x^2+y^2))/sqrt(x^2+y^2)',[-6*pi,6*pi]);
% view(0,75)
% lightangle(0,-90) %azimuth, elevation
lightangle(-60,20)
h.FaceLighting = 'gouraud'; h2.FaceLighting = 'gouraud'; h3.FaceLighting = 'gouraud'; h4.FaceLighting = 'gouraud';
h.AmbientStrength = 0.3; h2.AmbientStrength = 0.3; h3.AmbientStrength = 0.3; h4.AmbientStrength = 0.3;
h.DiffuseStrength = 0.8; h2.DiffuseStrength = 0.8; h3.DiffuseStrength = 0.8; h4.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9; h2.SpecularStrength = 0.9; h3.SpecularStrength = 0.9; h4.SpecularStrength = 0.9;
h.SpecularExponent = 25; h2.SpecularExponent = 25; h3.SpecularExponent = 25; h4.SpecularExponent = 25;
h.BackFaceLighting = 'unlit'; h2.BackFaceLighting = 'unlit'; h3.BackFaceLighting = 'unlit'; h4.BackFaceLighting = 'unlit';
axis equal;

colormap hot;
caxis([-100 100])
% set(gcf,'color','k'); %,'position',[100 100 2539 1878])
% axis off

fill3([-4 -4 2 2 -4],[5 -8 -8 5 5],-0.04*ones([1,5]),[0.99 0.75 0.0],'linestyle','none')

% for k = 1:50
%     a = 0.25*rand(1)+0.5;
%     x0 = 12*rand(1)-6; y0 = 12*rand(1)-6;
%     c1 = rand(1); c2 = rand(1); c3 = rand(1);
%     fill3(1.5*ones(size(psi)),a*cos(psi)-x0,a*sin(psi)-y0,[0 c2/2 c3],'linestyle','none');
%     
%     
% end

end
 

 
%% y(1) = psi, y(2) = psi_s, y(3) = r
%% y(4) = z, y(5) = gamma, y(6) = area
%% y(7) = mu, y(8) = L, y(9) = eta
function dydx = odesystem(x,y)
global kappa;
 
dydx = [ y(8)*y(2);
    y(8)/(kappa*y(3))*(-kappa*cos(y(1))*y(2) + kappa*sin(y(1))*cos(y(1))/y(3) + y(5)*sin(y(1)) + y(9)*cos(y(1)));
    y(8)*cos(y(1));
    -y(8)*sin(y(1));
    y(8)*(kappa/2*y(2)^2 - kappa/2*sin(y(1))^2/y(3)^2 + y(7));
    y(8)*2*pi*y(3);
    0;
    0;
    0];
end

function dydx = odesystem_sc(x,y)
global kappa;
global c0;
 
dydx = [ y(8)*y(2);
    y(8)/(kappa*y(3))*(-kappa*cos(y(1))*y(2) + kappa*sin(y(1))*cos(y(1))/y(3) + y(5)*sin(y(1)) + y(9)*cos(y(1)));
    y(8)*cos(y(1));
    -y(8)*sin(y(1));
    y(8)*(kappa/2*(y(2)-c0)^2 - kappa/2*sin(y(1))^2/y(3)^2 + y(7));
    y(8)*2*pi*y(3);
    0;
    0;
    0];
end
 
%% a  = right endpoint
%% b = left endpoint
function res = bcs(ya,yb)
global X;
global A;
global kappa;
global kbar;
 
res = [ ya(4) - X; %fixed extension
    yb(4) + X;
    ya(3) - 1; %fixed ring radius
    yb(3) - 1;
    ya(6); %area constraints
    yb(6) - A;
    kappa*(ya(2) + sin(ya(1))/ya(3)) - kbar*sin(ya(1))/ya(3); % no bending moment
    kappa*(yb(2) + sin(yb(1))/yb(3)) - kbar*sin(yb(1))/yb(3);
    ya(3)/2*(ya(2)^2 - sin(ya(1))^2/ya(3)^2) - ya(7)*ya(3) + ya(5)*cos(ya(1)) - ya(9)*sin(ya(1))]; %Hamiltonian = 0
end
 
function res = bcs_sc(ya,yb)
global X;
global A;
global kappa;
global kbar;
global c0;
 
res = [ ya(4) - X; %fixed extension
    yb(4) + X;
    ya(3) - 1; %fixed ring radius
    yb(3) - 1;
    ya(6); %area constraints
    yb(6) - A;
    kappa*(ya(2) + sin(ya(1))/ya(3) - c0) - kbar*sin(ya(1))/ya(3); % no bending moment
    kappa*(yb(2) + sin(yb(1))/yb(3) - c0) - kbar*sin(yb(1))/yb(3);
    ya(3)/2*(ya(2)^2 - (sin(ya(1))/ya(3) - c0)^2) - ya(7)*ya(3) + ya(5)*cos(ya(1)) - ya(9)*sin(ya(1))]; %Hamiltonian = 0
end

%% use previous solution stored in y_guess as initial guess
function v = newguess(q)
global N; 
global y_guess;
 
q = round(q*(N-1));
v = y_guess(:,q+1);
 
end
 
%% ode system with eigenvalue lambda
%% u(1) = u, u(2) = u_s, u(3) = u_ss, u(4) = u_sss, u(5) = p, u(6) = integral(r*2H*u)
function dudt = ams_ode_new(t,u,lambda,ysol,mu,L)

% EE = 1/2;
% FF = -mu/2 + (psi_s(t,ysol)^2 + H(t,ysol)^2 - K(t,ysol));
% G = cos(psi(t,ysol))^2*sin(psi(t,ysol))^2/2/r(t,ysol)^4 + sin(psi(t,ysol))^3*psi_s(t,ysol)/4/r(t,ysol)^3 ...
%     -cos(psi(t,ysol))^2*psi_s(t,ysol)^2/2/r(t,ysol)^2 + sin(psi(t,ysol))^2*psi_s(t,ysol)^2/2/r(t,ysol)^2 ...
%     -3*sin(psi(t,ysol))*psi_s(t,ysol)^3/4/r(t,ysol) -cos(psi(t,ysol))*sin(psi(t,ysol))*psi_ss(t,ysol)/2/r(t,ysol)^2 ...
%     +3*cos(psi(t,ysol))*psi_s(t,ysol)*psi_ss(t,ysol)/2/r(t,ysol) + psi_ss(t,ysol)^2/4 - sin(psi(t,ysol))*psi_sss(t,ysol)/4/r(t,ysol)...
%     +psi_s(t,ysol)*psi_sss(t,ysol)/4;
% GG = mu*K(t,ysol) + 2*(H(t,ysol)^2-K(t,ysol))*(4*H(t,ysol)^2-K(t,ysol)) + G;

AA = 1/2;
BB = -(1/(2*r(t,ysol))*(cos(psi(t,ysol))^2./r(t,ysol) + sin(psi(t,ysol))*psi_s(t,ysol)) + 2*K(t,ysol) + mu/2 + H(t,ysol)*(2*psi_s(t,ysol)-H(t,ysol)));
BBp = -cos(psi(t,ysol))^3/r(t,ysol)^3 + cos(psi(t,ysol))*sin(psi(t,ysol))^2/(2*r(t,ysol)^3) - 2.5*cos(psi(t,ysol))*sin(psi(t,ysol))*psi_s(t,ysol)/(2*r(t,ysol)^2) ...
      +cos(psi(t,ysol))*psi_s(t,ysol)^2/r(t,ysol) + sin(psi(t,ysol))*psi_ss(t,ysol)/r(t,ysol) - 2.5*psi_s(t,ysol)*psi_ss(t,ysol);
BBp = -BBp;
% G = cos(psi(t,ysol))^2*sin(psi(t,ysol))^2/r(t,ysol)^4 - cos(psi(t,ysol))^2*sin(psi(t,ysol))*psi_s(t,ysol)/r(t,ysol)^3 + sin(psi(t,ysol))^3*psi_s(t,ysol)/(2*r(t,ysol)^3) ...
%     -sin(psi(t,ysol))*psi_s(t,ysol)^3/(2*r(t,ysol)) + 2*cos(psi(t,ysol))*psi_s(t,ysol)*psi_ss(t,ysol)/r(t,ysol) + 1.5*psi_ss(t,ysol)^2 - sin(psi(t,ysol))*psi_sss(t,ysol)/(2*r(t,ysol)) + 1.5*psi_s(t,ysol)*psi_sss(t,ysol);
% CC = (mu*K(t,ysol) + 2*(H(t,ysol)^2-K(t,ysol))*(4*H(t,ysol)^2-K(t,ysol)) + G);
CC = cos(psi(t,ysol)).^2.*sin(psi(t,ysol)).^2/2./r(t,ysol).^4 + sin(psi(t,ysol)).^4/2./r(t,ysol).^4 - 3/4*cos(psi(t,ysol)).^2.*sin(psi(t,ysol)).*psi_s(t,ysol)./r(t,ysol).^3 ...
     + mu*sin(psi(t,ysol)).*psi_s(t,ysol)./r(t,ysol) - sin(psi(t,ysol)).^3.*psi_s(t,ysol)/4./r(t,ysol).^3 + cos(psi(t,ysol)).^2.*psi_s(t,ysol).^2/4./r(t,ysol).^2 - sin(psi(t,ysol)).^2.*psi_s(t,ysol).^2/4./r(t,ysol).^2 ...
     - sin(psi(t,ysol)).*psi_s(t,ysol).^3./r(t,ysol) + psi_s(t,ysol).^4/2 + cos(psi(t,ysol)).*sin(psi(t,ysol)).*psi_ss(t,ysol)/4./r(t,ysol).^2 + 2*cos(psi(t,ysol)).*psi_s(t,ysol).*psi_ss(t,ysol)./r(t,ysol) ...
     + 3/2*psi_ss(t,ysol).^2 - sin(psi(t,ysol)).*psi_sss(t,ysol)/2./r(t,ysol) + 3/2*psi_s(t,ysol).*psi_sss(t,ysol);
% BB = 0; BBp = 0;
% CC = 0;

dudt = [  L*u(2);
          L*u(3);
          L*u(4);
          -L/AA*(2*AA*cos(psi(t,ysol))/r(t,ysol)*u(4) + (-AA*sin(psi(t,ysol))*psi_s(t,ysol)/r(t,ysol) + BB)*u(3) + (BBp + cos(psi(t,ysol))/r(t,ysol)*BB)*u(2) + (CC - lambda)*u(1) + 2*H(t,ysol)*u(5));
          0;
          2*pi*L*r(t,ysol)*H(t,ysol)*u(1) ];
end

%% boundary conditions
function res = ams_bc_new(ua,ub,lambda,psi0,psiL)
global kbar;
res = [  ua(1) %variation does not perturb endpoints
         ub(1) 
         ua(3) + (1+kbar)*cos(psi0)*ua(2) %variation does not contribute bending moment
         ub(3) + (1+kbar)*cos(psiL)*ub(2)
         ua(6) %u orthogonal to 2H
         ub(6)
         ua(2)-1]; %set normalization
end

%% initial guess
function uinit = ams_init(t)

n = 2;
% uinit = [  rand
%            rand
%            rand
%            rand
%            rand
%            rand]
uinit = [  cos((t-0.5)*pi*n)
           -n*pi*sin((t-0.5)*pi*n)
           -(pi*n)^2*cos((t-0.5)*pi*n)
           (pi*n)^3*sin((t-0.5)*pi*n)
           -100
           0];
end

%% mean curvature at t
function out = H(t,ysol)
    y = deval(ysol,t);
    out = -(y(2) + sin(y(1))/y(3))/2;
end

%% Gaussian curvature at t
function out = K(t,ysol)
    y = deval(ysol,t);
    out = y(2)*sin(y(1))/y(3);
end

%% radial component
function out = r(t,ysol)
    y = deval(ysol,t);
    out = y(3);
end

%% angle of tangent vector
function out = psi(t,ysol)  
    y = deval(ysol,t);
    out = y(1);
end

%% dpsi/ds
function out = psi_s(t,ysol)
    y = deval(ysol,t);
    out = y(2);
end

%% d^2psi/ds^2
function out = psi_ss(t,ysol)
    y = deval(ysol,t);
    out = (-cos(y(1))*y(2) + sin(y(1))*cos(y(1))/y(3) + y(5)*sin(y(1)) +y(9)*cos(y(1)) )/y(3);
end

%% d^3psi/ds^3
function out = psi_sss(t,ysol)
    y = deval(ysol,t);
    out = -y(9)*cos(y(1))^2/y(3)^2 - 2*cos(y(1))^2*sin(y(1))/y(3)^3 - cos(y(1))*sin(y(1))*y(5)/y(3)^2 ...
        +sin(y(1))*gam_s(t,ysol)/y(3) + 2*cos(y(1))^2*y(2)/y(3)^2 - y(9)*sin(y(1))*y(2)/y(3) ...
        -sin(y(1))^2*y(2)/y(3)^2 + cos(y(1))*y(5)*y(2)/y(3) + sin(y(1))*y(2)^2/y(3) - cos(y(1))*psi_ss(t,ysol)/y(3);
end

%% d(gamma)/ds
function out = gam_s(t,ysol)
    y = deval(ysol,t);
    out = 1/2*y(2)^2 - sin(y(1))^2/2/y(3)^2 + y(7);
end

%% solve for catenoid with given area
%% x = [b,h]
function F = catarea(x)
global A;   

F(1) = A - pi*x(1)^2*(sinh(2*x(2)/x(1)) + 2*x(2)/x(1));
F(2) = 1 - x(1)*cosh(x(2)/x(1));
end

%% initial guess (m = 3, upper branch) %init interval linspace(0,1,9)
function out = guess3b(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -150.2; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(4*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2)+0.2;
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -162; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(4*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2)+0.2;
        -mu*b];
end

%% initial guess (m = 5b)
function out = guess5b(x)
global X;

    b = 0.82; %b = 0.825517;
    mu = -191; %mu = -67.9574;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.005*sin(5*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 3)
function out = guess3(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -120; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2)
function out = guess2(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -49.2; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2), nonzero kbar
function out = guess2_kbar(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -55; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.1*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 3), nonzero kbar
function out = guess3_kbar(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -31; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.1*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 3b), nonzero kbar
function out = guess3b_kbar(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -120; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.02*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2b), nonzero kbar
function out = guess2b_kbar(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -41; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.02*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4), nonzero kbar
function out = guess4_kbar(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = -1200; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.02*sin(4*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% explore higher order buckling modes
function [] = higherorder()

global kappa; kappa = 1;
global kbar; kbar = 0.0;
global A; A = 2*pi; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints

%% prepare a solution as initial guess
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.200001 0.1 0 pi 0 2 0]); %m = 5b
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %m = 5
%solinit = bvpinit(linspace(0,1,5),@guess4); %m = 4
solinit = bvpinit(linspace(0,1,11),@guess3); %m = 3
solinit = bvpinit(linspace(0,1,9),@guess3b); %m = 3b
%solinit = bvpinit(linspace(0,1,11),@guess2); %m = 2
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
%X = 0.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1)
%plot(z,r); hold on;

%% calculate shape for various widths
Xvec = [1.05:-0.05:-0.35]/2; %(2:-0.1:0.1)/2; 
%Xvec = [1.05 1.055]/2
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
%fig = figure('Position', [10 10 1500 750]); %offset from bottom left corner, [width height]
% set(gca,'FontSize',18);
% set(gcf,'color','w');
%hold on; grid on; axis equal; box on;
%v = VideoWriter('fifth_mode.avi');
%v.FrameRate = 10;
%open(v);
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); %r=fliplr(r);
    L = y(8,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    if j == length(Xvec) || j == 4 || j == 12 || j == 18 || j == 22 %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          plot(z,r,'LineWidth',2);    
%         Xvec(j)*2
%         
%         figure();
%         %[T,P] = meshgrid([1,1],linspace(0,2*pi,60));
%         [T,P] = meshgrid(1:N,linspace(0,2*pi,60));
%         C = (psi_s(T) + sin(psi(T))./r(T)).^2*kappa/2;
%         %r = [zeros(60,1) ones(60,1)]; z = zeros(60,2); C = zeros(60,2);
%         %s=surf(z,r.*cos(P),r.*sin(P),log(C)/log(10)); hold on;
%         s=surf(z(T),r(T).*cos(P),r(T).*sin(P),log(C)/log(10)); hold on;
%         cb = colorbar('southoutside'); cb.Label.String = 'log_{10} (\kappa/2)(2H)^2';
%         caxis([-2.1 2.1])
%         s.EdgeColor = 'none'; 
%         %s.FaceColor = 'interp';
%         %X = 0
%          plot3(ones([1,100])*X,cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',4);
%          plot3(-ones([1,100])*X,cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',4);
%         axis equal;
%         axis([-0.6 0.6 -1.35 1.35 -1.35 1.35]); %box on; 
%         axis off;
%         set(gcf,'color','w');set(gca,'FontSize',18);
%         
%         for q = 1:2:60
%             plot3(z,r.*cos(2*pi*q/60),r.*sin(2*pi*q/60),'k');
%             %plot3([0 0],[0,1*cos(2*pi*q/60)],[0,1*sin(2*pi*q/60)],'k')
%         end
%         break
    end
  
%% --------------------------- movie script
%     %% subplot 1: shape
%     subplot(1,3,1);
%     ax = gca;
%     ax.NextPlot = 'replaceChildren';
%     %axis([-0.6 0.6 -1.25 1.25])
%    
%     [T,P] = meshgrid(1:N,linspace(0,2*pi,24*10));
%     C = (psi_s(T) + sin(psi(T))./r(T)).^2*kappa/2;
%     s=surf(z(T),r(T).*cos(P),r(T).*sin(P),C); hold on;
%     %colorbar
%     s.EdgeColor = 'none'; s.FaceColor = 'interp';
%     plot3(ones([1,100])*X,cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k.','LineWidth',4);
%     plot3(-ones([1,100])*X,cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k.','LineWidth',4);
%     axis([-0.6 0.6 -1.25 1.25 -1.25 1.25]); %box on; 
%     axis off;
%     %title('Shape');
%     %% subplot 2: force
%     subplot(1,3,2);
    Fvec(j) = -2*pi*y(9,1);
%     plot(Xvec(1:j),Fvec(1:j),'bo-','LineWidth',2)
%     title('Force');
%     %axis([0 0.6 -420 -150]); %m = 2
%     axis([0 0.6 -3200 -1000]) %m = 5
%     xlabel('$$h/a$$','Interpreter','latex'); ylabel('$$F/\kappa a$$','Interpreter','latex');
%     grid on
%     %b = r((N+1)/2); k1 = psi_s((N+1)/2);
%     %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
%     %% subplot 3: tension
%     subplot(1,3,3);
     muvec(j) = y(7,1);
%     plot(Xvec(1:j),muvec(1:j),'bo-','LineWidth',2);
%     title('Tension');
%     %axis([0 0.6 -30 60]) %m = 2
%     axis([0 0.6 -200 360]); %m = 5
%     xlabel('$$h/a$$','Interpreter','latex'); ylabel('$$\mu a^2/\kappa$$','Interpreter','latex');
%     grid on;
%     Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
%    
%     frame = getframe(fig);
%     writeVideo(v,frame)
 %% --------------------- end movie script
     hold on;  grid on;
    %plot(z,r,'b','LineWidth',2);
    y_guess = y;
   
end
return
X = 1.04/2;
solinit = bvpinit(linspace(0,1,11),@guess3); %m = 3
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %m = 5
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
y = y_guess;
Xvec2 = [1.05:-0.05:0.0]/2; %(2:-0.1:0.1)/2; 
%Xvec = [1.05 1.055]/2
muvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    L = y(8,1);
    Fvec2(j) = -2*pi*y(9,1);
    muvec2(j) = y(7,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
     hold on;  grid on;
    %plot(z,r,'r','LineWidth',2);
    y_guess = y;
   
end



% close(v);
% title('Force vs. extension')
% xlabel('z'); ylabel('r');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% bcat = 0.825517; v = linspace(-0.527697,0.527697);
% plot(v,bcat*cosh(v/bcat),'k--','LineWidth',2)

figure();
subplot(1,3,1)
plot([ Xvec],[Fvec],'bo-','LineWidth',2); hold on;
plot([ Xvec2],[Fvec2],'ro-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18);
 
subplot(1,3,2)
plot([Xvec],[muvec],'bo-','LineWidth',2); hold on;
plot([Xvec2],[muvec2],'ro-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
%plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
%plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
%plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
%plot([0.527697 0.527697],[-50 200],'k--');
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

end

%% remove area constraint, add stretching modulus
function [] = variablearea()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A0; A0 = 2*pi; %preferred area, does not have to match
global A; A = 2*pi*2; %initial area, needs to be big enough so that surface exists given X
global M; M = 10000; %stretching modulus
global X; X = 0.55; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints

%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 2 0]);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec = 0.45:0.005:0.6;
Evec = zeros(size(Xvec));
for jj = 1:length(Xvec)
%bcat = 0.825517; Xcat = 0.527697;
X = Xvec(jj);

b = fzero(@(b) b*cosh(X/b) - 1,1);
A_lb = pi*b^2*(sinh(2*X/b) + 2*X/b);
[Amin,Emin] = fmincon(@(A) minFixedA(A) + 0.5*M*(A-A0)^2/A0^2,A0*2,[],[],[],[],A_lb,100);
Evec(jj) = Emin;
end

plot(Xvec,Evec); grid on; hold on;
title(['M = ' num2str(M)]);
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18);
end

%% auxiliary solver
function bendingE = minFixedA(A_in)
global A; A = A_in;
global kappa;
global kbar;
global N;
global M;
global X;
global y_guess;
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    bendingE = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r);
    
end

%% creates a figure showing the thin solutions (continuation of
%% TRP's solutions for smaller separations)
function [] = thinsolutions()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
 
%% prepare a solution as initial guess
%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %wacky
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi 0 2 0]); %also wacky, A = 1.11*2*pi
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]);

% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.200001 0.1 0 pi 0 2 0]); %m = 5b
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %m = 5
% solinit = bvpinit(linspace(0,1,5),@guess4); %m = 4
% solinit = bvpinit(linspace(0,1,11),@guess3); %m = 3
% solinit = bvpinit(linspace(0,1,9),@guess3b); %m = 3b
% solinit = bvpinit(linspace(0,1,11),@guess2); %m = 2
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
X = 0.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
% X = 0.25; solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 0.08 0.1 0 pi 0 1.2 -80]); 
% X = 0.25; solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 0.09 0.1 0 pi 0 1.2 -80]); 
% X = 0.25; solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 0.30002 0.1 0 pi 5 1.2 15]); 
% X = 0.25; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.13 0.1 0 pi 50 1 15]); 
% X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 -20 0.11 0.1 0 pi 70 1.2 70]); 
% X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 -20 0.21 0.1 0 pi 70 1.2 70]); 
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
gam = y(5,:);
mu = y(7,1)
L = y(8,1)
eta = y(9,1)
% plot(z,r); return;
% plot(z,psi_s + sin(psi)./r); return;

% psi_ss = 1./r.*(-cos(psi).*psi_s + sin(psi).*cos(psi)./r + gam.*sin(psi) + eta*cos(psi));
% C = psi_ss.^2 + psi_s.^4/4 - mu*psi_s.^2;
% plot(t,C)
% %plot(z,r);
% return;



%% calculate shape for various widths
% Xvec = [1.04:-0.05:0.2]/2; %(2:-0.1:0.1)/2; 
Xvec = [ 1.05:-0.05:0.5]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
evvec = zeros(1,length(Xvec));
hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1); %b = r((N+1)/2);
    muvec(j) = mu;
     %k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,2) == 1
%     plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'b-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'b-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-.','linewidth',2); hold on;
        end
    end
end
 

%% prepare a solution as initial guess
X = 0.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
% X = 0.52; solinit = bvpinit(linspace(0,1,11),@guess2); %m = 2
% X = 0.52; solinit = bvpinit(linspace(0,1,9),@guess3b); %m = 3b
% solinit = bvpinit(linspace(0,1,11),@guess3); %m = 3
% solinit = bvpinit(linspace(0,1,5),@guess4); %m = 4
% X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.200001 0.1 0 pi 0 2 0]); %m = 5b
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %m = 5
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec2 =  [ 1.05:-0.05:0.0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
evvec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
%     if mod(j,2) ==1
%     plot(z,r)
%     plot(t,psi_s,t(1:end-1),diff(psi)./diff(t)/L,'o');
%       plot(t,cos(psi),t(1:end-1),diff(r)./diff(t)/L,'o');
%     end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec2(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end
subplot(1,2,1);
set(gca,'fontsize',18);
set(gcf,'color','w');
title('Shapes');
grid on;
xlabel('z/a'); ylabel('r/a');
axis equal;

subplot(1,2,2);
set(gca,'fontsize',18);
set(gcf,'color','w');
title('Minimizing eigenfunctions');
xlabel('t'); ylabel('u/a');

figure();
subplot(1,3,1);
plot(Xvec,Fvec,'bo-','linewidth',2); grid on; hold on;
plot(Xvec2,Fvec2,'ro-','linewidth',2);
set(gca,'fontsize',18);
set(gcf,'color','w');
ylabel('Force F')
subplot(1,3,2);
plot(Xvec,muvec,'bo-','linewidth',2); grid on; hold on;
plot(Xvec2,muvec2,'ro-','linewidth',2);
set(gca,'fontsize',18)
ylabel('Tension \mu')
subplot(1,3,3);
plot(Xvec,evvec,'bo-','linewidth',2); grid on; hold on;
plot(Xvec2,evvec2,'ro-','linewidth',2);
set(gca,'fontsize',18)
ylabel('\lambda_m_i_n');

return

% muvec'
% muvec2'
% Fvec'
% Fvec2'
% 
% %% prepare a solution as initial guess
% solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
% options = bvpset('RelTol',1e-5);
% sol = bvp4c(@odesystem,@bcs,solinit,options);
% t = linspace(0,1,N);
% y_guess = deval(sol,t);
%  
% Xvec3 = [1.05 1.0543 1.0539 1.054 1.0543 1.05433 1.05434 1.05435 1.05436]/2;
% muvec3 = zeros(1,length(Xvec3));
% Fvec3 = zeros(1,length(Xvec3));
% for j = 1:length(Xvec3)
%     X = Xvec3(j);
%     solinit = bvpinit(linspace(0,1,11),@newguess);
%     sol = bvp4c(@odesystem,@bcs,solinit);
%     t = linspace(0,1,N);
%     y = deval(sol,t);
%     %r = y(3,:); z = y(4,:); psi_s = y(2,:);
%     muvec3(j) = y(7,1);
%     %b = r((N+1)/2); k1 = psi_s((N+1)/2);
%     %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
%     Fvec3(j) = -2*pi*y(9,1);
%     hold on;  grid on;
%    % plot(z,r,'LineWidth',2);
%    % y_guess = y;
% end

set(gca,'FontSize',18);
set(gcf,'color','w');
 
%% catenoid solution
bcat = 0.825517; v = linspace(-0.527697,0.527697);
plot(v,bcat*cosh(v/bcat),'k','LineWidth',3)
% plot(v./(v.*v+bcat^2*cosh(v/bcat).^2)*(v(1)^2+bcat^2*cosh(v(1)/bcat)^2),...
%     bcat*cosh(v/bcat)./(v.*v+bcat^2*cosh(v/bcat).^2)*(v(1)^2+bcat^2*cosh(v(1)/bcat)^2),'k--');
xlabel('z'); ylabel('r');
title('Shapes (n = 4)');

figure();
subplot(1,2,1)
plot(Xvec,Fvec,'bo-','linewidth',2); hold on;
% plot(Xvec2,Fvec2,'ro-','linewidth',2); hold on;
% plot([Xvec3(end) Xvec],[Fvec3(end) Fvec],'bo-','LineWidth',2); hold on;
% plot([Xvec3(end) Xvec2],[Fvec3(end) Fvec2],'ro-','LineWidth',2); hold on;
% plot(Xvec3(end),Fvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18);
 
subplot(1,2,2)
plot(Xvec,muvec,'bo-','linewidth',2); hold on;
% plot(Xvec2,muvec2,'ro-','linewidth',2); hold on;
% plot([Xvec3(end) Xvec],[muvec3(end) muvec],'bo-','LineWidth',2); hold on;
% plot([Xvec3(end) Xvec2],[muvec3(end) muvec2],'ro-','LineWidth',2); hold on;
% plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

% subplot(1,3,3)
% plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
% plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3(end),0,'kp','MarkerSize',18,'MarkerFaceColor','k'); hold on;
%plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
%plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
%plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
% %plot([0.527697 0.527697],[-50 200],'k--');
% title('Energy vs. extension')
% xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
end

%% shapes for area = 2*pi*1.1
function [] = thinsolutions2()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
 
%% prepare a solution as initial guess
% X = 0.2/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether9);
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]);
% X = 1;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.22 0.1 0 pi 600
% 3.7 -40]); %A = 1.05*2*pi
% X = 1;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.32 0.1 0 pi 100
% 3.7 -40]); %tether
X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess4_tether9); %n =4
% N = 1001; X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
gam = y(5,:);
mu = y(7,1)
L = y(8,1)
eta = y(9,1)

plot(z,r)

% return

%% some tether branches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xvec = [1:0.01:0.57*2 0.585*2 0.59*2 0.5901*2]/2;
% Xvec = [0.45:0.005:0.59];
% Xvec = [1:-0.05:0.55];
Xvec = [0.59:-0.005:0];
% Xvec = [0.59:-0.005:0];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
evvec = zeros(1,length(Xvec));
hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1); %b = r((N+1)/2);
    muvec(j) = mu;
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    plot(z,r)
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec(j) = eval;
    if mod(j,20) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare a solution as initial guess
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether9); %n = 2, thick
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.215 0.1 0 pi -10 2 0]);
N = 1001; X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess4_tether10); %n =4
% N = 1001; X = 0.466; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths
% Xvec2 =  [0.45:-0.005:0];
% Xvec2 = [0.52:-0.005:0.42];
Xvec2 = [0.52:-0.005:0.295];
% Xvec2 = [0.52:-0.005:0.0];
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
evvec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;

    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec2(j) = eval;
    if mod(j,20) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'b-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'b-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-.','linewidth',2); hold on;
        end
    end
end

%% prepare a solution as initial guess
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]);
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether12); %n = 2, thin
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]);
N = 1001; X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.352 0.1 0 pi -10 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
% Xvec3 = [0.5901 0.59:-0.005:0.525 0.522];
% Xvec3 = [0.45:-0.005:0.33];
Xvec3 = [0.59:-0.005:0];
% Xvec3 = [0.52:-0.005:0.73/2];
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
evvec3 = zeros(1,length(Xvec3));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;

    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec3(j) = eval;
    if mod(j,20) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'m-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'m-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'m-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'m-.','linewidth',2); hold on;
        end
    end
end

%% prepare a solution as initial guess (thin cat to tether)
% X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi 100 2 0]);
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether12); %n = 2, thin
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.113 0.1 0 pi -10 2 0]);
N = 1001; X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.142 0.1 0 pi -150 2 70]); 
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
% Xvec4 = [0.523 0.525:0.005:1];
Xvec4 = [0.45:0.005:0.52];
% Xvec4 = [0.52:-0.005:0];
% Xvec4 = [0.59:-0.005:0];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
evvec4 = zeros(1,length(Xvec4));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec4(j) = eval;
    if mod(j,20) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'g-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'g-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'g-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'g-.','linewidth',2); hold on;
        end
    end
end



subplot(1,2,1);
set(gca,'fontsize',18);
set(gcf,'color','w');
title('Shapes');
grid on;
xlabel('z/a'); ylabel('r/a');
axis equal;

subplot(1,2,2);
set(gca,'fontsize',18);
set(gcf,'color','w');
title('Minimizing eigenfunctions');
xlabel('t'); ylabel('u/a');

figure();
subplot(1,3,1);
plot(Xvec,Fvec,'ro-','linewidth',2); grid on; hold on;
plot(Xvec2,Fvec2,'bo-','linewidth',2);
plot(Xvec3,Fvec3,'mo-','linewidth',2);
plot(Xvec4,Fvec4,'go-','linewidth',2);
set(gca,'fontsize',18);
set(gcf,'color','w');
ylabel('Force F')
subplot(1,3,2);
plot(Xvec,muvec,'ro-','linewidth',2); grid on; hold on;
plot(Xvec2,muvec2,'bo-','linewidth',2);
plot(Xvec3,muvec3,'mo-','linewidth',2);
plot(Xvec4,muvec4,'go-','linewidth',2);
set(gca,'fontsize',18)
ylabel('Tension \mu')
subplot(1,3,3);
plot(Xvec,evvec,'ro-','linewidth',2); grid on; hold on;
plot(Xvec2,evvec2,'bo-','linewidth',2);
plot(Xvec3,evvec3,'mo-','linewidth',2);
plot(Xvec4,evvec4,'go-','linewidth',2);
set(gca,'fontsize',18)
ylabel('\lambda_m_i_n');
end

%% shapes for area = 2*pi*1.3
function [] = thinsolutions3()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
 
%% prepare a solution as initial guess
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 0 1.7 0]);
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6);
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
% X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether6); %n =4
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); 
    y = deval(sol,t);
%     r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
%     H = (psi_s + sin(psi)./r)/2;
%     muvec = y(7,1);
%     L = y(8,1);
%     Fvec = -2*pi*y(9,1);
%     indices = find([0 diff(sign(H))]~=0);
%     plot(z,r,'b','linewidth',2);
%     Xvec = 0.59;
%     Evec = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
%         %% stability analysis
%     [eval,efun]= evmatrix(sol);
%     evvec = eval;
%         if eval > 0
%         subplot(1,2,1)
%         plot(z,r,'r-','linewidth',2); hold on;
%         subplot(1,2,2)
%         plot(t,efun,'r-','linewidth',2); hold on;
%         else
%         subplot(1,2,1)
%         plot(z,r,'r-.','linewidth',2); hold on;
%         subplot(1,2,2)
%         plot(t,efun,'r-.','linewidth',2); hold on;
%         end
% options = bvpset('RelTol',1e-3);
% sol = bvp4c(@odesystem,@bcs,solinit,options);
% %Xvec = [0.58 0.56 0.54];
% Xvec = [Xvec 0.61 0.65 0.68 0.7 0.705 0.706 0.707 0.708 0.709];

%% calculate shape for various widths
% Xvec = [0:0.01:1.5];
% Xvec = [0.6:0.01:0.7 0.701:0.001:0.712];
Xvec = [1.04:0.01:1.4 1.401:0.001:1.431 1.4312]/2;
% Xvec = [1.05:0.01:1.40 1.401:0.001:1.415 1.4161]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
evvec = zeros(1,length(Xvec));
% muvec = [muvec zeros(1,length(Xvec)-1)];
% Fvec = [Fvec zeros(1,length(Xvec)-1)];
% Evec = [Evec zeros(1,length(Xvec)-1)];
% evvec = [evvec zeros(1,length(Xvec)-1)];
hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1); %b = r((N+1)/2);
    muvec(j) = mu;
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;

    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec(j) = eval;
    if mod(j,10) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end
 

%% prepare a solution as initial guess
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6);
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
% X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether6); %n =4
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths
% Xvec2 = [0.6:-0.01:0.47];
Xvec2 = [1.04:-0.01:0]/2;
% Xvec2 = [0.58 0.57 0.56];
% Xvec2 = [1.05:-0.01:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
evvec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;

    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec2(j) = eval;
    if mod(j,10) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end

%% prepare a solution as initial guess
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5);
X = 1.24/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.1 0.11 0 pi -3 1.9 1.4]); %n=3T
% X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether5); %n=4
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
% Xvec3 = [0.6:0.01:0.7 0.701:0.001:0.712];
Xvec3 = [1.04:0.01:1.41 1.411 1.412 1.413 1.414 1.4155]/2;
% Xvec3 = [0.59:0.01:0.7 0.705 0.706 0.707 0.708 0.7085 0.7088];
% Xvec3 = [1.04:0.01:1.4 1.401:0.001:1.417 1.4172 1.4174 1.4176 1.4178]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
evvec3 = zeros(1,length(Xvec3));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;

    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec3(j) = eval;
    if mod(j,10) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'g-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'g-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'g-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'g-.','linewidth',2); hold on;
        end
    end
end

%% prepare a solution as initial guess (thin cat to tether)
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5);
X = 1.24/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.1 0.11 0 pi -3 1.9 1.4]); %n=3T
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
% Xvec4 = [0.6:-0.01:0];
Xvec4 = [1.04:-0.01:0]/2;
% Xvec4 = [0.59:-0.01:0];
% Xvec4 = [1.04:-0.01:0]/2;
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
evvec4 = zeros(1,length(Xvec4));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec4(j) = eval;
    if mod(j,10) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'g-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'g-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'g-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'g-.','linewidth',2); hold on;
        end
    end
end

%% prepare a solution as initial guess (thin cat to tether)
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5);
X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.14 0.1 0 pi 0 1.7 0]); %n=3L
% X = 1.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.13 0.1 0 pi -70 1.3 90]); %n=5S
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
% Xvec5 = [0.6:-0.01:0];
Xvec5 = [1.04:0.01:1.4 1.401:0.001:1.415 1.4155]/2;
% Xvec5 = [1.4178 1.4174 1.417 1.415:-0.001:1.401 1.4:-0.01:1.3 1.295:-0.005:1.225 1.2]/2;
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
Evec5 = zeros(1,length(Xvec5));
evvec5 = zeros(1,length(Xvec5));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec5)
    X = Xvec5(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec5(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec5(j) = -2*pi*y(9,1);
    Evec5(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec5(j) = eval;
    if mod(j,10) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'b-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'b-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-.','linewidth',2); hold on;
        end
    end
end

%% prepare a solution as initial guess (thin cat to tether)
% % X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5);
X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.14 0.1 0 pi 0 1.7 0]); %n=3L
% X = 1.4/2; solinit = bvpinit(linspace(0,1,9),[pi/2 2 0.14 0.1 0 pi -70 1.3 90]); %n=5T
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
% Xvec6 = [0.6:-0.01:0];
% Xvec6 = [1.04:-0.01:0]/2;
Xvec6 = [1.415:-0.001:1.4 1.395:-0.01:1.315 1.3149 1.3148 ]/2;
muvec6 = zeros(1,length(Xvec6));
Fvec6 = zeros(1,length(Xvec6));
Evec6 = zeros(1,length(Xvec6));
evvec6 = zeros(1,length(Xvec6));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec6)
    X = Xvec6(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec6(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec6(j) = -2*pi*y(9,1);
    Evec6(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec6(j) = eval;
    if mod(j,10) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'m-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'m-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'m-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'m-.','linewidth',2); hold on;
        end
    end
end
 
%% prepare a solution as initial guess (thin cat to tether)
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5);
X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -1 1.7 1]); %n=3S
% X = 1.4/2; solinit = bvpinit(linspace(0,1,9),[pi/2 2 0.14 0.1 0 pi -70 1.3 90]); %n=5T
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
% Xvec4 = [0.6:-0.01:0];
Xvec7 = [1.4311 1.43:-0.001:1.401 1.4:-0.01:1.28]/2;
% Xvec7 = [1.415 1.2989]/2;
muvec7 = zeros(1,length(Xvec7));
Fvec7 = zeros(1,length(Xvec7));
Evec7 = zeros(1,length(Xvec7));
evvec7 = zeros(1,length(Xvec7));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec7)
    X = Xvec7(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec7(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec7(j) = -2*pi*y(9,1);
    Evec7(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec7(j) = eval;
    if mod(j,10) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'m-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'m-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'m-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'m-.','linewidth',2); hold on;
        end
    end
end


subplot(1,2,1);
set(gca,'fontsize',18);
set(gcf,'color','w');
title('Shapes');
grid on;
xlabel('z/a'); ylabel('r/a');
axis equal;

subplot(1,2,2);
set(gca,'fontsize',18);
set(gcf,'color','w');
title('Minimizing eigenfunctions');
xlabel('t'); ylabel('u/a');

figure();
subplot(1,3,1);
plot(Xvec,Fvec,'ro-','linewidth',2); grid on; hold on;
plot(Xvec2,Fvec2,'ro-','linewidth',2);
plot(Xvec3,Fvec3,'go-','linewidth',2);
plot(Xvec4,Fvec4,'go-','linewidth',2);
plot(Xvec5,Fvec5,'bo-','linewidth',2);
plot(Xvec6,Fvec6,'mo-','linewidth',2);
plot(Xvec7,Fvec7,'mo-','linewidth',2);
set(gca,'fontsize',18);
set(gcf,'color','w');
ylabel('Force F')
subplot(1,3,2);
plot(Xvec,muvec,'ro-','linewidth',2); grid on; hold on;
plot(Xvec2,muvec2,'ro-','linewidth',2);
plot(Xvec3,muvec3,'go-','linewidth',2);
plot(Xvec4,muvec4,'go-','linewidth',2);
plot(Xvec5,muvec5,'bo-','linewidth',2);
plot(Xvec6,muvec6,'mo-','linewidth',2);
plot(Xvec7,muvec7,'mo-','linewidth',2);
set(gca,'fontsize',18)
ylabel('Tension \mu')
subplot(1,3,3);
plot(Xvec,evvec,'ro-','linewidth',2); grid on; hold on;
plot(Xvec2,evvec2,'ro-','linewidth',2);
plot(Xvec3,evvec3,'go-','linewidth',2);
plot(Xvec4,evvec4,'go-','linewidth',2);
plot(Xvec5,evvec5,'bo-','linewidth',2);
plot(Xvec6,evvec6,'mo-','linewidth',2);
plot(Xvec7,evvec7,'mo-','linewidth',2);
set(gca,'fontsize',18)
ylabel('\lambda_m_i_n');
end

%% creates a figure showing the thin solutions (area = 2*pi*1.1)
function [] = thinsolutions2_blah()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
 
%% prepare a solution as initial guess (1)
%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1, upper
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %n =7, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi 0 2 0]); %also wacky, A = 1.11*2*pi
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %<--- n=1, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.37 0.1 0 pi -10 2 0]); %n = 7, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 10 2 0]); %tether
%X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 10 2 0]); %tether
%X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi 100 2 0]); %tether
% X = 0.55; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 100 2 0]); %3rd, big
%X = 0.2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 100 2 0]); %negative extension?
% X = 0.525; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.39 0.1 0 pi -18.5 1.93 15.4]);
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -18.5 1.93 15.4]); %double kink
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.51 0.1 0 pi -18.5 2 15.4]); %double kink
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether9); %n = 2,thick
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether10); %upper kink, increase mu by 1 for lower
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether11);
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether12); %n = 2, thin
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess4_tether9); %n =4
% N = 1001; X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess4_tether10); %n =4

X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n= 5, lower
%  X = 0.466; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.514; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.352 0.1 0 pi -10 2 0]); %n = 5, upper
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.142 0.1 0 pi -150 2 70]); %n = 5 lower
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -150 2 70]); %n = 6
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.117 0.1 0 pi -150 2 70]); %n = 8
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.124 0.1 0 pi -150 2 70]); %n = 12

 options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
y(7,1)
y(8,1)
y(9,1)
H = (psi_s + sin(psi)./r)/2;
indices = find([0 diff(sign(H))]~=0);

plot(z,r); hold on; plot(z(indices),r(indices),'o'); return;

%% calculate shape for various widths (catenoid lower)
disp('branch 1')
% Xvec = [0.45:0.005:0.524];
Xvec = [0.466:0.005:0.59 0.5905];
%Xvec =  [0.4:0.005:0.52]; 
%Xvec = [0.5:0.0005:0.5905]; %(2:-0.1:0.1)/2; 
%Xvec = [0.56*2:-0.001:0.52*2]/2;
%Xvec = [1.0:-0.01:0.58];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(); hold on; axis equal; box on;
v1 = linspace(-0.590668192994804,0.590668192994804);
v2 = linspace(-0.521378393180893,0.521378393180893);
b1 = 0.757966837819559;
b2 = 0.255417422349254;
plot(v1,b1*cosh(v1/b1),'g--','linewidth',3); hold on;
plot(v2,b2*cosh(v2/b2),'g--','linewidth',3);
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
        plot(z,r,'c','LineWidth',2);
    end
    y_guess = y;
end
%return

%% prepare a solution as initial guess (2)
X = 0.2/2;
solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths (catenoid upper)
disp('branch 2')
Xvec2 = [1:0.01:0.57*2 0.585*2 0.59*2 0.5901*2]/2; %[ 1.15:-0.05:0]/2;
%Xvec2 = [1.049:-0.05:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,2) == 0
        plot(z,r,'r','LineWidth',2);
    end
%     z2= z./(r.*r+z.*z)*(1+z(end)^2);
%     r2 = r./(r.*r+z.*z)*(1+z(end)^2);
%     plot(z2,r2,'r--');
    %Psi = acos(((r.*r+z.*z).*cos(psi) - 2*r.*(r.*cos(psi) - z.*sin(psi)))./(r.*r+z.*z).^2);
    %plot(t,Psi);
   % y_guess = y;
end

%% prepare a solution as initial guess (3)
disp('branch 3')
X = 1.04/2;
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec3 = [1.0:0.001:0.59*2]/2; %(2:-0.1:0.1)/2; ; %[1.0:-0.05:0.5];
%Xvec3 = [1.0:-0.01:0]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,40) == 0
    plot(z,r,'k','LineWidth',2);
    end
    y_guess = y;
end

%% prepare a solution as initial guess (4)
disp('branch 4')
X = 1.04/2;
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec4 = [1.0:-0.01:0]/2; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 0
    plot(z,r,'m','LineWidth',2);
    end
    y_guess = y;
end

%% prepare a solution as initial guess (5)
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.35 0.1 0 pi -10 2 0]);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
disp('branch 5')
Xvec5 = [0.54:-0.005:0]; %[0.5:-0.01:0]; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
%Xvec5 = [0.55:0.005:0.59]
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
Evec5 = zeros(1,length(Xvec5));
hold on; axis equal; box on;
for j = 1:length(Xvec5)
    X = Xvec5(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec5(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec5(j) = -2*pi*y(9,1);
    Evec5(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 0
    plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
end

%% prepare a solution as initial guess (6)
%X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.27 0.1 0 pi -19 2 2]);
 X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether11);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
disp('branch 6')
Xvec6 = [0.52:0.005:0.59 0.5905]; %[0.5:-0.01:0]; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
muvec6 = zeros(1,length(Xvec6));
Fvec6 = zeros(1,length(Xvec6));
Evec6 = zeros(1,length(Xvec6));
hold on; axis equal; box on;
for j = 1:length(Xvec6)
    X = Xvec6(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec6(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec6(j) = -2*pi*y(9,1);
    Evec6(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 0
        plot(z,r,'y','LineWidth',2);
    end
    y_guess = y;
end

%% prepare a solution as initial guess (7)
X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi 100 2 0]); %tether
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
disp('branch 7')
Xvec7 = [0.4:0.005:0.52]; %[0.5:-0.01:0]; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
muvec7 = zeros(1,length(Xvec7));
Fvec7 = zeros(1,length(Xvec7));
Evec7 = zeros(1,length(Xvec7));
hold on; axis equal; box on;
for j = 1:length(Xvec7)
    X = Xvec7(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec7(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec7(j) = -2*pi*y(9,1);
    Evec7(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 0
    plot(z,r,'LineWidth',2,'color',[0.4660, 0.6740, 0.1880]);
    end
    y_guess = y;
end

%% prepare a solution as initial guess (7)
X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi 100 2 0]); %tether
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
disp('branch 8')
Xvec8 = [0.525:0.005:0.8]; %[0.5:-0.01:0]; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
muvec8 = zeros(1,length(Xvec8));
Fvec8 = zeros(1,length(Xvec8));
Evec8 = zeros(1,length(Xvec8));
hold on; axis equal; box on;
for j = 1:length(Xvec8)
    X = Xvec8(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec8(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec8(j) = -2*pi*y(9,1);
    Evec8(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 0
    plot(z,r,'LineWidth',2,'color',[0.9290, 0.6940, 0.1250]);
    end
    y_guess = y;
end

set(gca,'FontSize',18);
set(gcf,'color','w');

figure();
subplot(1,3,1)
plot([Xvec],[ Fvec],'co-','LineWidth',2); hold on;
plot([Xvec2],[Fvec2],'ro-','LineWidth',2); hold on;
plot(Xvec3,Fvec3,'ko-','linewidth',2)
plot(Xvec4,Fvec4,'mo-','linewidth',2)
plot(Xvec5,Fvec5,'bo-','linewidth',2)
plot(Xvec6,Fvec6,'yo-','linewidth',2)
plot(Xvec7,Fvec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
plot(Xvec8,Fvec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18);
 
subplot(1,3,2)
plot([ Xvec],[muvec],'co-','LineWidth',2); hold on;
plot([ Xvec2],[muvec2],'ro-','LineWidth',2); hold on;
plot(Xvec3,muvec3,'ko-','linewidth',2)
plot(Xvec4,muvec4,'mo-','linewidth',2)
plot(Xvec5,muvec5,'bo-','linewidth',2)
plot(Xvec6,muvec6,'yo-','linewidth',2)
plot(Xvec7,muvec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
plot(Xvec8,muvec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

subplot(1,3,3)
plot(Xvec,Evec,'co-','LineWidth',2); hold on;
plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
plot(Xvec3,Evec3,'ko-','linewidth',2); hold on;
plot(Xvec4,Evec4,'mo-','linewidth',2); hold on;
plot(Xvec5,Evec5,'bo-','linewidth',2); hold on;
plot(Xvec6,Evec6,'yo-','linewidth',2); hold on;
plot(Xvec7,Evec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
plot(Xvec8,Evec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
%plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
%plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
%plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
%plot([0.527697 0.527697],[-50 200],'k--');
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');
end

function [] = thinsolutions2_kbar2()
 
global kappa; kappa = 1;
global kbar; kbar = 0.1;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
 
%% prepare a solution as initial guess 
% X = 0.5; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.113 0.1 0 pi -10 2 0]); %n = 3
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.14 0.1 0 pi -16 2.5 3]); %n = 3
% X = 0.55; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 100 2 0]); %n = 5
% X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.12 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.12 0.1 0 pi -16 2.5 3]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.11 0.1 0 pi -16 2.5 3]); %n = 7
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.08 0.1 0 pi -16 2.5 3]); %n = 11

% X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess4_tether12);
% X = 0.55; solinit = bvpinit(linspace(0,1,6),@guess4_tether11);
% X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess2_tether13);
X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess2_tether14);

 options = bvpset('RelTol',1e-8);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
y(7,1)
y(8,1)
y(9,1)
H = (psi_s + sin(psi)./r)/2;
indices = find([0 diff(sign(H))]~=0);

% plot(z,r); hold on; plot(z(indices),r(indices),'o'); return;

%% calculate shape for various widths (catenoid lower)
disp('branch 1')
Xvec = [0.5 0.52:0.001:0.5905 0.5909];
% Xvec = [0.5:0.005:0.52 0.521 0.5213 0.52131];
% Xvec = [0.46:0.005:0.52 0.521 0.5211 0.5212 0.5213 0.52134 ]; 
%Xvec =  [0.4:0.005:0.52]; 
%Xvec = [0.5:0.0005:0.5905]; %(2:-0.1:0.1)/2; 
%Xvec = [0.56*2:-0.001:0.52*2]/2;
%Xvec = [1.0:-0.01:0.58];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(); hold on; axis equal; box on;
v1 = linspace(-0.590668192994804,0.590668192994804);
v2 = linspace(-0.521378393180893,0.521378393180893);
b1 = 0.757966837819559;
b2 = 0.255417422349254;
plot(v1,b1*cosh(v1/b1),'g--','linewidth',3); hold on;
plot(v2,b2*cosh(v2/b2),'g--','linewidth',3);
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    H = (psi_s + sin(psi)./r)/2;
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
        plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
end

figure
plot(Xvec,Fvec,'o-'); grid on;
return

%% prepare a solution as initial guess (2)
X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess2_tether14);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths (catenoid upper)
disp('branch 2')
Xvec2 = [1:0.01:0.58*2 0.585*2 0.588*2 0.59*2  0.5902*2:0.000001:0.59061*2]/2; %[ 1.15:-0.05:0]/2;
% Xvec2 = [0.5:-0.005:0.03 0.0029];
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    H = (psi_s + sin(psi)./r)/2;
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 1
        plot(z,r,'r','LineWidth',2);
    end
%     z2= z./(r.*r+z.*z)*(1+z(end)^2);
%     r2 = r./(r.*r+z.*z)*(1+z(end)^2);
%     plot(z2,r2,'r--');
    %Psi = acos(((r.*r+z.*z).*cos(psi) - 2*r.*(r.*cos(psi) - z.*sin(psi)))./(r.*r+z.*z).^2);
    %plot(t,Psi);
   % y_guess = y;
end

figure
plot(Xvec,muvec,'o-',Xvec2,muvec2,'o-'); grid on;
return

%% prepare a solution as initial guess (3)
disp('branch 3')
X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.12 0.1 0 pi -16 2.5 3]); %n = 5
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec3 = [0.46:0.005:0.52 0.521 0.5211 0.5212 0.5213 0.52137 ]; %(2:-0.1:0.1)/2; ; %[1.0:-0.05:0.5];
%Xvec3 = [1.0:-0.01:0]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
    plot(z,r,'k','LineWidth',2);
    end
    y_guess = y;
end

%% prepare a solution as initial guess (4)
disp('branch 4')
 X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.12 0.1 0 pi -10 2 0]); %n = 5
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec4 = [0.46:0.005:0.52 0.521 0.5211 0.5212 0.5213 0.52136 ]; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
    plot(z,r,'m','LineWidth',2);
    end
    y_guess = y;
end

figure()
plot(Xvec,muvec,'o-','linewidth',2); hold on; grid on; box on;
plot(Xvec2,muvec2,'o-','linewidth',2)
plot(Xvec3,muvec3,'o-','linewidth',2)
plot(Xvec4,muvec4,'o-','linewidth',2)
plot([0.521378393180893 0.521378393180893],[-40 -160],'k--')
plot([0.590668192994804 0.590668192994804],[-40 -160],'k--')
set(gca,'fontsize',18)
set(gcf,'color','w')
title(['Tension vs. Extension (5th mode, \kappabar/\kappa = ' num2str(kbar/kappa) ')']);
xlabel('Half-width h/2')
ylabel('Tension \mu')


% figure();
% subplot(1,3,1)
% plot([Xvec],[ Fvec],'co-','LineWidth',2); hold on;
% plot([Xvec2],[Fvec2],'ro-','LineWidth',2); hold on;
% plot(Xvec3,Fvec3,'ko-','linewidth',2)
% plot(Xvec4,Fvec4,'mo-','linewidth',2)
% plot(Xvec5,Fvec5,'bo-','linewidth',2)
% plot(Xvec6,Fvec6,'yo-','linewidth',2)
% plot(Xvec7,Fvec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,Fvec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% title('Force vs. extension')
% xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
% set(gca,'FontSize',18);
%  
% subplot(1,3,2)
% plot([ Xvec],[muvec],'co-','LineWidth',2); hold on;
% plot([ Xvec2],[muvec2],'ro-','LineWidth',2); hold on;
% plot(Xvec3,muvec3,'ko-','linewidth',2)
% plot(Xvec4,muvec4,'mo-','linewidth',2)
% plot(Xvec5,muvec5,'bo-','linewidth',2)
% plot(Xvec6,muvec6,'yo-','linewidth',2)
% plot(Xvec7,muvec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,muvec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% title('Tension vs. extension')
% xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% 
% subplot(1,3,3)
% plot(Xvec,Evec,'co-','LineWidth',2); hold on;
% plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Evec3,'ko-','linewidth',2); hold on;
% plot(Xvec4,Evec4,'mo-','linewidth',2); hold on;
% plot(Xvec5,Evec5,'bo-','linewidth',2); hold on;
% plot(Xvec6,Evec6,'yo-','linewidth',2); hold on;
% plot(Xvec7,Evec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,Evec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% %plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
% %plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
% %plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
% %plot([0.527697 0.527697],[-50 200],'k--');
% title('Energy vs. extension')
% xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
end

%% creates a figure showing the thin solutions
function [] = thinsolutions2_kbar3()
 
global kappa; kappa = 1;
global kbar; kbar = 0.1;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
 
%% prepare a solution as initial guess 
% X = 0.5; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.215 0.1 0 pi -10 2 0]); %n = 3, upper
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.113 0.1 0 pi -10 2 0]); %n = 3
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.14 0.1 0 pi -16 2.5 3]); %n = 3
% X = 0.55; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 100 2 0]); %n = 5
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.12 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.12 0.1 0 pi -16 2.5 3]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.11 0.1 0 pi -16 2.5 3]); %n = 7
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.08 0.1 0 pi -16 2.5 3]); %n = 11



%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1, upper
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %n =7, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi 0 2 0]); %also wacky, A = 1.11*2*pi
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %<--- n=1, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.37 0.1 0 pi -10 2 0]); %n = 7, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 10 2 0]); %tether
%X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 10 2 0]); %tether
%X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi 100 2 0]); %tether
% X = 0.55; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 100 2 0]); %3rd, big
%X = 0.2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 100 2 0]); %negative extension?
% X = 0.525; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.39 0.1 0 pi -18.5 1.93 15.4]);
%X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -18.5 1.93 15.4]); %double kink
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.51 0.1 0 pi -18.5 2 15.4]); %double kink
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether9); %n = 2
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether10); %upper kink, increase mu by 1 for lower
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether11);
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether12);
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess4_tether9); %n =4
% N = 1001; X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess4_tether10); %n =4
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n= 5, lower
%  X = 0.466; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.514; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.352 0.1 0 pi -10 2 0]); %n = 5, upper
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.142 0.1 0 pi -150 2 70]); %n = 5 lower
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -150 2 70]); %n = 6
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.117 0.1 0 pi -150 2 70]); %n = 8
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.124 0.1 0 pi -150 2 70]); %n = 12

 options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
y(7,1)
y(8,1)
y(9,1)
H = (psi_s + sin(psi)./r)/2;
indices = find([0 diff(sign(H))]~=0);

% plot(z,r); hold on; plot(z(indices),r(indices),'o'); return;

%% calculate shape for various widths (catenoid lower)
disp('branch 1')
Xvec = [0.45:0.005:0.524];
% Xvec = [0.466:0.005:0.59 0.5901 0.5905 0.5906];
%Xvec =  [0.4:0.005:0.52]; 
%Xvec = [0.5:0.0005:0.5905]; %(2:-0.1:0.1)/2; 
%Xvec = [0.56*2:-0.001:0.52*2]/2;
%Xvec = [1.0:-0.01:0.58];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
evvec = zeros(1,length(Xvec));
figure(); hold on; axis equal; box on;
v1 = linspace(-0.590668192994804,0.590668192994804);
v2 = linspace(-0.521378393180893,0.521378393180893);
b1 = 0.757966837819559;
b2 = 0.255417422349254;
plot(v1,b1*cosh(v1/b1),'g--','linewidth',3); hold on;
plot(v2,b2*cosh(v2/b2),'g--','linewidth',3);
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 1
        plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
    
     %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec(j) = eval;
    if mod(j,2) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'b-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'b-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-.','linewidth',2); hold on;
        end
    end
end


%% prepare a solution as initial guess (2)
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths (catenoid upper)
disp('branch 2')
Xvec2 = [1:0.01:0.58*2 0.585*2 0.588*2 0.59*2  0.5902*2 0.5905*2 0.5906*2]/2; %[ 1.15:-0.05:0]/2;
%Xvec2 = [1.049:-0.05:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
evvec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
        plot(z,r,'r','LineWidth',2);
    end
%     z2= z./(r.*r+z.*z)*(1+z(end)^2);
%     r2 = r./(r.*r+z.*z)*(1+z(end)^2);
%     plot(z2,r2,'r--');
    %Psi = acos(((r.*r+z.*z).*cos(psi) - 2*r.*(r.*cos(psi) - z.*sin(psi)))./(r.*r+z.*z).^2);
    %plot(t,Psi);
    y_guess = y;
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec2(j) = eval;
    if mod(j,2) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end


%% prepare a solution as initial guess (3)
disp('branch 3')
%X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.14 0.1 0 pi -16 2.5 3]); %n = 3
X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec3 = [0.52:0.005:0.59]; %(2:-0.1:0.1)/2; ; %[1.0:-0.05:0.5];
%Xvec3 = [1.0:-0.01:0]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
evvec3 = zeros(1,length(Xvec3));
hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
    plot(z,r,'k','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec3(j) = eval;
    if mod(j,2) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'g-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'g-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'g-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'g-.','linewidth',2); hold on;
        end
    end
end

%% prepare a solution as initial guess (4)
disp('branch 4')
X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.113 0.1 0 pi -10 2 0]); %n = 3
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec4 = [0.45:0.005:0.52 0.521 0.5211 0.5212 0.5213 0.52136 ]; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
evvec4 = zeros(1,length(Xvec3));
hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
    plot(z,r,'m','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec4(j) = eval;
    if mod(j,2) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'m-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'m-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'m-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'m-.','linewidth',2); hold on;
        end
    end
end

figure()
plot(Xvec,muvec,'o-','linewidth',2); hold on; grid on; box on;
plot(Xvec2,muvec2,'o-','linewidth',2)
plot(Xvec3,muvec3,'o-','linewidth',2)
plot(Xvec4,muvec4,'o-','linewidth',2)
plot([0.521378393180893 0.521378393180893],[0 -70],'k--')
plot([0.590668192994804 0.590668192994804],[0 -70],'k--')
set(gca,'fontsize',18)
set(gcf,'color','w')
title(['Tension vs. Extension (3rd mode, \kappabar/\kappa = ' num2str(kbar/kappa) ')']);
xlabel('Half-width h/2')
ylabel('Tension \mu')

figure();
plot(Xvec,evvec,'bo-','linewidth',2); hold on; grid on;
plot(Xvec2,evvec2,'ro-','linewidth',2)
plot(Xvec3,evvec3,'go-','linewidth',2)
plot(Xvec4,evvec4,'mo-','linewidth',2)
set(gca,'fontsize',18)
set(gcf,'color','w')

% figure();
% subplot(1,3,1)
% plot([Xvec],[ Fvec],'co-','LineWidth',2); hold on;
% plot([Xvec2],[Fvec2],'ro-','LineWidth',2); hold on;
% plot(Xvec3,Fvec3,'ko-','linewidth',2)
% plot(Xvec4,Fvec4,'mo-','linewidth',2)
% plot(Xvec5,Fvec5,'bo-','linewidth',2)
% plot(Xvec6,Fvec6,'yo-','linewidth',2)
% plot(Xvec7,Fvec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,Fvec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% title('Force vs. extension')
% xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
% set(gca,'FontSize',18);
%  
% subplot(1,3,2)
% plot([ Xvec],[muvec],'co-','LineWidth',2); hold on;
% plot([ Xvec2],[muvec2],'ro-','LineWidth',2); hold on;
% plot(Xvec3,muvec3,'ko-','linewidth',2)
% plot(Xvec4,muvec4,'mo-','linewidth',2)
% plot(Xvec5,muvec5,'bo-','linewidth',2)
% plot(Xvec6,muvec6,'yo-','linewidth',2)
% plot(Xvec7,muvec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,muvec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% title('Tension vs. extension')
% xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% 
% subplot(1,3,3)
% plot(Xvec,Evec,'co-','LineWidth',2); hold on;
% plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Evec3,'ko-','linewidth',2); hold on;
% plot(Xvec4,Evec4,'mo-','linewidth',2); hold on;
% plot(Xvec5,Evec5,'bo-','linewidth',2); hold on;
% plot(Xvec6,Evec6,'yo-','linewidth',2); hold on;
% plot(Xvec7,Evec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,Evec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% %plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
% %plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
% %plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
% %plot([0.527697 0.527697],[-50 200],'k--');
% title('Energy vs. extension')
% xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
end

%% creates a figure showing the thin solutions
function [] = thinsolutions2_kbar4()
 
global kappa; kappa = 1;
global kbar; kbar = 0.1;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
 
%% prepare a solution as initial guess 
% X = 0.5; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.113 0.1 0 pi -10 2 0]); %n = 3
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.14 0.1 0 pi -16 2.5 3]); %n = 3
% X = 0.55; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 100 2 0]); %n = 5
% X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.12 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.12 0.1 0 pi -16 2.5 3]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.11 0.1 0 pi -16 2.5 3]); %n = 7
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.08 0.1 0 pi -16 2.5 3]); %n = 11

X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess4_tether12);
% X = 0.55; solinit = bvpinit(linspace(0,1,6),@guess4_tether13);

options = bvpset('RelTol',1e-8);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
y(7,1)
y(8,1)
y(9,1)
H = (psi_s + sin(psi)./r)/2;
indices = find([0 diff(sign(H))]~=0);

% plot(z,r); hold on; plot(z(indices),r(indices),'o'); return;

%% calculate shape for various widths (catenoid lower)
disp('branch 1')
% Xvec = [0.45:0.005:0.524];
Xvec = [0.5:0.005:0.52 0.521 0.5213];
%Xvec =  [0.4:0.005:0.52]; 
%Xvec = [0.5:0.0005:0.5905]; %(2:-0.1:0.1)/2; 
%Xvec = [0.56*2:-0.001:0.52*2]/2;
%Xvec = [1.0:-0.01:0.58];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(); hold on; axis equal; box on;
v1 = linspace(-0.590668192994804,0.590668192994804);
v2 = linspace(-0.521378393180893,0.521378393180893);
b1 = 0.757966837819559;
b2 = 0.255417422349254;
plot(v1,b1*cosh(v1/b1),'g--','linewidth',3); hold on;
plot(v2,b2*cosh(v2/b2),'g--','linewidth',3);
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 1
        plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
end


%% prepare a solution as initial guess (2)
%X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5 
X = 0.55; solinit = bvpinit(linspace(0,1,6),@guess4_tether13);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths (catenoid upper)
disp('branch 2')
Xvec2 = [0.55 0.56 0.57 0.58 0.59 0.5901:0.0001:0.5906 0.590601:0.000001:0.590654 0.5906545];
% Xvec2 = [1:0.01:0.58*2 0.585*2 0.588*2 0.59*2  0.5902*2 0.5905*2 0.5906*2]/2; %[ 1.15:-0.05:0]/2;
%Xvec2 = [1.049:-0.05:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
        plot(z,r,'r','LineWidth',2);
    end
%     z2= z./(r.*r+z.*z)*(1+z(end)^2);
%     r2 = r./(r.*r+z.*z)*(1+z(end)^2);
%     plot(z2,r2,'r--');
    %Psi = acos(((r.*r+z.*z).*cos(psi) - 2*r.*(r.*cos(psi) - z.*sin(psi)))./(r.*r+z.*z).^2);
    %plot(t,Psi);
   % y_guess = y;
end


%% prepare a solution as initial guess (3)
disp('branch 3')
X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.12 0.1 0 pi -16 2.5 3]); %n = 5
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec3 = [0.46:0.005:0.52 0.521 0.5211 0.5212 0.5213 0.52137 ]; %(2:-0.1:0.1)/2; ; %[1.0:-0.05:0.5];
%Xvec3 = [1.0:-0.01:0]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
    plot(z,r,'k','LineWidth',2);
    end
    y_guess = y;
end

%% prepare a solution as initial guess (4)
disp('branch 4')
 X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.12 0.1 0 pi -10 2 0]); %n = 5
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec4 = [0.46:0.005:0.52 0.521 0.5211 0.5212 0.5213 0.52136 ]; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
    plot(z,r,'m','LineWidth',2);
    end
    y_guess = y;
end

figure()
plot(Xvec,muvec,'o-','linewidth',2); hold on; grid on; box on;
plot(Xvec2,muvec2,'o-','linewidth',2)
plot(Xvec3,muvec3,'o-','linewidth',2)
plot(Xvec4,muvec4,'o-','linewidth',2)
plot([0.521378393180893 0.521378393180893],[-40 -160],'k--')
plot([0.590668192994804 0.590668192994804],[-40 -160],'k--')
set(gca,'fontsize',18)
set(gcf,'color','w')
title(['Tension vs. Extension (5th mode, \kappabar/\kappa = ' num2str(kbar/kappa) ')']);
xlabel('Half-width h/2')
ylabel('Tension \mu')


% figure();
% subplot(1,3,1)
% plot([Xvec],[ Fvec],'co-','LineWidth',2); hold on;
% plot([Xvec2],[Fvec2],'ro-','LineWidth',2); hold on;
% plot(Xvec3,Fvec3,'ko-','linewidth',2)
% plot(Xvec4,Fvec4,'mo-','linewidth',2)
% plot(Xvec5,Fvec5,'bo-','linewidth',2)
% plot(Xvec6,Fvec6,'yo-','linewidth',2)
% plot(Xvec7,Fvec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,Fvec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% title('Force vs. extension')
% xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
% set(gca,'FontSize',18);
%  
% subplot(1,3,2)
% plot([ Xvec],[muvec],'co-','LineWidth',2); hold on;
% plot([ Xvec2],[muvec2],'ro-','LineWidth',2); hold on;
% plot(Xvec3,muvec3,'ko-','linewidth',2)
% plot(Xvec4,muvec4,'mo-','linewidth',2)
% plot(Xvec5,muvec5,'bo-','linewidth',2)
% plot(Xvec6,muvec6,'yo-','linewidth',2)
% plot(Xvec7,muvec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,muvec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% title('Tension vs. extension')
% xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% 
% subplot(1,3,3)
% plot(Xvec,Evec,'co-','LineWidth',2); hold on;
% plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Evec3,'ko-','linewidth',2); hold on;
% plot(Xvec4,Evec4,'mo-','linewidth',2); hold on;
% plot(Xvec5,Evec5,'bo-','linewidth',2); hold on;
% plot(Xvec6,Evec6,'yo-','linewidth',2); hold on;
% plot(Xvec7,Evec7,'o-','linewidth',2,'color',[0.4660, 0.6740, 0.1880])
% plot(Xvec8,Evec8,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
% %plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
% %plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
% %plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
% %plot([0.527697 0.527697],[-50 200],'k--');
% title('Energy vs. extension')
% xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
end

%% creates a figure showing the thin solutions
function [] = thinsolutions2_kbar5()
 
global kappa; kappa = 1;
global kbar; kbar = 0.1;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
 
%% prepare a solution as initial guess 
% X = 0.5; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.113 0.1 0 pi -10 2 0]); %n = 3
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.14 0.1 0 pi -16 2.5 3]); %n = 3
X = 0.55; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 100 2 0]); %n = 5
% X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.12 0.1 0 pi -10 2 0]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.12 0.1 0 pi -16 2.5 3]); %n = 5
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.11 0.1 0 pi -16 2.5 3]); %n = 7
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.08 0.1 0 pi -16 2.5 3]); %n = 11


%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1, upper
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %n =7, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi 0 2 0]); %also wacky, A = 1.11*2*pi
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %<--- n=1, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.37 0.1 0 pi -10 2 0]); %n = 7, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 10 2 0]); %tether
%X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 10 2 0]); %tether
%X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi 100 2 0]); %tether
% X = 0.55; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 100 2 0]); %3rd, big
%X = 0.2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 100 2 0]); %negative extension?
% X = 0.525; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.39 0.1 0 pi -18.5 1.93 15.4]);
%X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -18.5 1.93 15.4]); %double kink
% X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.51 0.1 0 pi -18.5 2 15.4]); %double kink
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether9); %n = 2
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether10); %upper kink, increase mu by 1 for lower
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether11);
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether12);
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess4_tether9); %n =4
% N = 1001; X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess4_tether10); %n =4
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n= 5, lower
%  X = 0.466; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.514; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.352 0.1 0 pi -10 2 0]); %n = 5, upper
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.142 0.1 0 pi -150 2 70]); %n = 5 lower
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -150 2 70]); %n = 6
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.117 0.1 0 pi -150 2 70]); %n = 8
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.124 0.1 0 pi -150 2 70]); %n = 12

 options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
y(7,1)
y(8,1)
y(9,1)
H = (psi_s + sin(psi)./r)/2;
indices = find([0 diff(sign(H))]~=0);

% plot(z,r); hold on; plot(z(indices),r(indices),'o'); return;

%% calculate shape for various widths (catenoid lower)
disp('branch 1')
% Xvec = [0.45:0.005:0.524];
Xvec = [0.466:0.005:0.59 0.5901 0.5905 0.5906];
%Xvec =  [0.4:0.005:0.52]; 
%Xvec = [0.5:0.0005:0.5905]; %(2:-0.1:0.1)/2; 
%Xvec = [0.56*2:-0.001:0.52*2]/2;
%Xvec = [1.0:-0.01:0.58];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(); hold on; axis equal; box on;
v1 = linspace(-0.590668192994804,0.590668192994804);
v2 = linspace(-0.521378393180893,0.521378393180893);
b1 = 0.757966837819559;
b2 = 0.255417422349254;
plot(v1,b1*cosh(v1/b1),'g--','linewidth',3); hold on;
plot(v2,b2*cosh(v2/b2),'g--','linewidth',3);
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 1
        plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
end


%% prepare a solution as initial guess (2)
X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths (catenoid upper)
disp('branch 2')
Xvec2 = [1:0.01:0.58*2 0.585*2 0.588*2 0.59*2  0.5902*2 0.5905*2 0.5906*2]/2; %[ 1.15:-0.05:0]/2;
%Xvec2 = [1.049:-0.05:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
        plot(z,r,'r','LineWidth',2);
    end
%     z2= z./(r.*r+z.*z)*(1+z(end)^2);
%     r2 = r./(r.*r+z.*z)*(1+z(end)^2);
%     plot(z2,r2,'r--');
    %Psi = acos(((r.*r+z.*z).*cos(psi) - 2*r.*(r.*cos(psi) - z.*sin(psi)))./(r.*r+z.*z).^2);
    %plot(t,Psi);
   % y_guess = y;
end


%% prepare a solution as initial guess (3)
disp('branch 3')
X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.12 0.1 0 pi -16 2.5 3]); %n = 5
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec3 = [0.46:0.005:0.52 0.521 0.5211 0.5212 0.5213 0.52137 ]; %(2:-0.1:0.1)/2; ; %[1.0:-0.05:0.5];
%Xvec3 = [1.0:-0.01:0]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
    plot(z,r,'k','LineWidth',2);
    end
    y_guess = y;
end

%% prepare a solution as initial guess (4)
disp('branch 4')
 X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.12 0.1 0 pi -10 2 0]); %n = 5
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shapes for various widths (tethers)
Xvec4 = [0.46:0.005:0.52 0.521 0.5211 0.5212 0.5213 0.52136 ]; %(2:-0.1:0.1)/2;  %[1.0:-0.05:0.5];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 0
    plot(z,r,'m','LineWidth',2);
    end
    y_guess = y;
end

figure()
plot(Xvec,muvec,'o-','linewidth',2); hold on; grid on; box on;
plot(Xvec2,muvec2,'o-','linewidth',2)
plot(Xvec3,muvec3,'o-','linewidth',2)
plot(Xvec4,muvec4,'o-','linewidth',2)
plot([0.521378393180893 0.521378393180893],[-40 -160],'k--')
plot([0.590668192994804 0.590668192994804],[-40 -160],'k--')
set(gca,'fontsize',18)
set(gcf,'color','w')
title(['Tension vs. Extension (5th mode, \kappabar/\kappa = ' num2str(kbar/kappa) ')']);
xlabel('Half-width h/2')
ylabel('Tension \mu')



end

%% creates a figure showing the thin solutions
function [] = catsolutions2_kbar2()
 
global kappa; kappa = 1;
global kbar; kbar = 0.1;
global A; A = 2*pi; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
 
%% prepare a solution as initial guess 
X = 0.5; solinit = bvpinit(linspace(0,1,11),@guess2_tether15); 


 options = bvpset('RelTol',1e-10);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
y(7,1)
y(8,1)
y(9,1)
H = (psi_s + sin(psi)./r)/2;
indices = find([0 diff(sign(H))]~=0);

% plot(z,r); hold on; plot(z(indices),r(indices),'o'); return;

%% calculate shape for various widths (catenoid lower)
disp('branch 1')
% Xvec = [0.45:0.005:0.524];
Xvec = [0.466:0.005:0.52 0.521 0.522 0.523 0.524 0.525 0.526 0.527 0.5271 0.5272 0.5273 0.5274 0.5275 0.5276 0.52765];
%Xvec =  [0.4:0.005:0.52]; 
%Xvec = [0.5:0.0005:0.5905]; %(2:-0.1:0.1)/2; 
%Xvec = [0.56*2:-0.001:0.52*2]/2;
%Xvec = [1.0:-0.01:0.58];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(); hold on; axis equal; box on;
v1 = linspace(-0.590668192994804,0.590668192994804);
v2 = linspace(-0.521378393180893,0.521378393180893);
b1 = 0.757966837819559;
b2 = 0.255417422349254;
plot(v1,b1*cosh(v1/b1),'g--','linewidth',3); hold on;
plot(v2,b2*cosh(v2/b2),'g--','linewidth',3);
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,4) == 1
        plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
end


%% prepare a solution as initial guess (2)
X = 0.52; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi -10 2 0]); %n = 5
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
figure()
plot(Xvec,muvec,'o-','linewidth',2); hold on; grid on; box on;
set(gca,'fontsize',18)
set(gcf,'color','w')
title(['Tension vs. Extension (2nd mode, \kappabar/\kappa = ' num2str(kbar/kappa) ')']);
xlabel('Half-width h/2')
ylabel('Tension \mu')





end

%% finds the Willmore torus (known genus 1 minimizer)
%% force on outer = - force on inner, tension equals zero
function [] = torussolutions()
 
global kappa; kappa = 1;
global kbar; kbar = 0.1;
global A; A = 10.2591; %10.1439; %2*pi*sqrt(2)*2*pi*2-10.1439; %10.1439; %surface area
global X; X = 0.498454; %1/2;
global y_guess;
global N; N = 101; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 2 0.1 0 pi 0 2 0]); %inner solution
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 1.1 0.1 0 pi 0 2*pi*sqrt(2)*.75 0]); %outer solution
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 

%% calculate shape for various widths
Xvec = 0.498454*2;%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
% figure; hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
%     hold on;  grid on;
%     plot(z,r,'LineWidth',2);
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    wt = snapshot(1.5,1.35,psi,psi_s,z,r,X,kappa,N);
    wt = centercrop(wt,0.2,0.2); imwrite(wt,'willmore_torus.tif','tif','resolution',600,'compression','none');
    y_guess = y;
end
 


% plot3(ones([1,100])*X,cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',4);
% plot3(-ones([1,100])*X,cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',4);


% %% torus solution
% hold on;
% rw = 1.49278; %sqrt(2);
% v = rw*linspace(-1,1,1001);
% tp = rw*sqrt(2)+sqrt(rw^2-v.*v);
% tn = rw*sqrt(2)-sqrt(rw^2-v.*v);
% plot(v,tp,'k--',v,tn,'k--');
% set(gcf,'Color','w');
% set(gca,'FontSize',18);
% axis([-1.5 1.5 0 3.5]); axis equal; 
% wt = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N); hold on;
% wt = centercrop(wt,0.2,0.2);  imwrite(wt,'willmore_torus.png');
% 
return
figure();
subplot(1,2,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
%plot(Xvec2,Fvec2,'bo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
 
subplot(1,2,2)
plot(Xvec,muvec,'ro-','LineWidth',2); hold on;
%plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
%plot(Xvec2,Fvec2)
 
 
end

%% finds the sphere (kappa = kbar/2)
function [] = spheresolutions()
 
global kappa; kappa = 1;
global kbar; kbar = 2;
global X; X = 0.68; %1/2;
rbar = sqrt(1 + X^2)
psibar = asin(1/rbar)
global A; A = 4*pi*rbar^2*sin((pi/2-asin(1/rbar))) %10.1439; %2*pi*sqrt(2)*2*pi*2-10.1439; %10.1439; %surface area
global y_guess;
global N; N = 101; %number of gridpoints
global c0; c0 = 0;


 
%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 0.2 0.0 0 pi/2 0 2*rbar*abs(pi/2-psibar) -10]); %inner solution
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 1.1 0.1 0 pi 0 2*pi*sqrt(2)*.75 0]); %outer solution
%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %wacky
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 2 0]);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
y = y_guess; r = y(3,:); z = y(4,:);
mu = y(7,1)
F = -2*pi*y(9,1)
plot(z,r); hold on;
v = linspace(-X,X,101);
plot(v,sqrt(rbar^2-v.^2),'o'); axis equal
return;

%% calculate shape for various widths
Xvec = 0.498454*2;%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
figure; hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    hold on;  grid on;
    plot(z,r,'LineWidth',2);
    psi_s
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end
 
%% torus solution
hold on;
rw = 1.49278; %sqrt(2);
v = rw*linspace(-1,1,1001);
tp = rw*sqrt(2)+sqrt(rw^2-v.*v);
tn = rw*sqrt(2)-sqrt(rw^2-v.*v);
plot(v,tp,'k--',v,tn,'k--');
set(gcf,'Color','w');
set(gca,'FontSize',18);
axis([-1.5 1.5 0 3.5]); axis equal; 

figure();
subplot(1,2,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
%plot(Xvec2,Fvec2,'bo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
 
subplot(1,2,2)
plot(Xvec,muvec,'ro-','LineWidth',2); hold on;
%plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
%plot(Xvec2,Fvec2)
 
 
end

%% summary figure 3D surface
function [] = summaryfigure()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
 
%% prepare a solution as initial guess
options = bvpset('RelTol',1e-3);
t = linspace(0,1,N);

figure();
xvals = []; xvals2 = []; xvals3 = []; xvals4 = []; xvals5 = []; xvals6 = [];
Avals = []; Avals2 = []; Avals3 = []; Avals4 = []; Avals5 = []; Avals6 = [];
Fvals = []; Fvals2 = []; Fvals3 = []; Fvals4 = []; Fvals5 = []; Fvals6 = [];

% A = 2*pi*1.01;
% Abar = A/2/pi
% bhans =  fsolve(@catarea,[0.1,1]);
% bcat = bhans(1)
% xcat = bhans(2)
% if Abar > 1 && Abar < 1.2 %in this regime, there are two catenoids
%     bhans2 = fsolve(@catarea,[0.1,0.1]);
%     bcat2 = bhans2(1)
%     xcat2 = bhans2(2)
% end
% X = 0.25;  solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 0.16 0.1 0 pi 200 3.73 -100]); %A = 1.05*2*pi
% sol = bvp4c(@odesystem,@bcs,solinit,options);
% y_guess = deval(sol,t); y = y_guess;
% plot(y(4,:),y(3,:))
% mu = y(7,1)
% L = y(8,1)
% eta = y(9,1)
% % return
% 
% %% calculate shape for various widths (lower branch)
% Xvec3 = [0.225 0.22 0.215 0.21 0.205];
% Fvec3 = zeros(1,length(Xvec3));
% for j = 1:length(Xvec3)
%     X = Xvec3(j)
%     solinit = bvpinit(linspace(0,1,11),@newguess);
%     sol = bvp4c(@odesystem,@bcs,solinit);
%     t = linspace(0,1,N);
%     y = deval(sol,t);
%     r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
%     mu = y(7,1) %b = r((N+1)/2);
%     Fvec3(j) = -2*pi*y(9,1)
%     plot(z,r); hold on;
%     y_guess = y;
% end
% 
% return

%% case 1: Abar < 1 (upper and lower branches)
Abarvec = 0.90:0.01:1.0;% 0.1:0.1:1.3;
for Ael = Abarvec*2*pi
    A = Ael;
    Abar = Ael/2/pi
    if Ael < 1.2*2*pi
    bhans =  fsolve(@catarea,[0.1,1]);
    bcat = bhans(1)
    xcat = bhans(2)
    else
        bcat = 0.9;
        xcat = 0.6;
    end
    X = xcat - 0.01;
    solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 bcat 0.1 0 pi -10 2 0]); %m = 1
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    y_guess = deval(sol,t); y = y_guess;

    %% calculate shape for various widths (lower branch)
    Xvec = [ 2*xcat - 0.01:-0.025:0.0 0]/2;
    Fvec = zeros(1,length(Xvec));
    for j = 1:length(Xvec)
        X = Xvec(j);
        solinit = bvpinit(linspace(0,1,11),@newguess);
        sol = bvp4c(@odesystem,@bcs,solinit);
        t = linspace(0,1,N);
        y = deval(sol,t);
        r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
        mu = y(7,1); %b = r((N+1)/2);
        Fvec(j) = -2*pi*y(9,1);
    
        y_guess = y;  
    end
    if A < 1.2*2*pi %add catenoid data
        Fcat = 2*pi*catstability(bcat)*bcat; hold on; grid on; box on;
    end
    xvals = [xvals Xvec xcat]; Avals = [Avals A*ones(size(Xvec))/2/pi A/2/pi]; Fvals = [Fvals Fvec Fcat]; 
    


    X = xcat - 0.01;
    solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 bcat/2 0.1 0 pi -10 2 0]); %m = 1
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    y_guess = deval(sol,t); y = y_guess;

    %% calculate shape for various widths (upper branch)
    if A == 2*pi*1
        Xvec2 = [2*xcat - 0.01:-0.025:0.3]/2;
    else
        Xvec2 = [2*xcat - 0.01:-0.025:0.0 0]/2;
    end
    Fvec2 = zeros(1,length(Xvec2));
    for j = 1:length(Xvec2)
        X = Xvec2(j);
        solinit = bvpinit(linspace(0,1,11),@newguess);
        sol = bvp4c(@odesystem,@bcs,solinit);
        t = linspace(0,1,N);
        y = deval(sol,t);
        r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
        mu = y(7,1); %b = r((N+1)/2);
        Fvec2(j) = -2*pi*y(9,1);

        y_guess = y;
    end
    if A == 2*pi*1 %manually fill in points when Abar = 1 (critical case)
        xvals2 = [xvals2 0 0.05 0.1]; Fvals2 = [Fvals2 0 -0.0001 -0.001]; Avals2 = [Avals2 1 1 1];
    end  
    xvals2 = [xvals2 Xvec2 xcat]; Avals2 = [Avals2 A*ones(size(Xvec2))/2/pi A/2/pi]; Fvals2 = [Fvals2 Fvec2 Fcat]; 

end

%% case 2: Abar > 1.2 (only one branch)
Abarvec = 1.2:0.01:1.4;
for Ael = Abarvec*2*pi
    A = Ael;
    Abar = Ael/2/pi
    if Ael < 1.2*2*pi
    bhans =  fsolve(@catarea,[0.1,1]);
    bcat = bhans(1)
    xcat = bhans(2)
    else
        bcat = 0.9;
        xcat = 0.6;
    end
    X = xcat - 0.01;
    solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 bcat 0.1 0 pi -10 2 0]); %m = 1
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    y_guess = deval(sol,t); y = y_guess;

    %% calculate shape for various widths (lower branch)
    Xvec = [ 2*xcat - 0.01:-0.025:0.0 0]/2;
    muvec = zeros(1,length(Xvec));
    Fvec = zeros(1,length(Xvec));
    %Evec = zeros(1,length(Xvec));
    for j = 1:length(Xvec)
        X = Xvec(j);
        solinit = bvpinit(linspace(0,1,11),@newguess);
        sol = bvp4c(@odesystem,@bcs,solinit);
        t = linspace(0,1,N);
        y = deval(sol,t);
        r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
        mu = y(7,1); %b = r((N+1)/2);
        muvec(j) = mu;
        Fvec(j) = -2*pi*y(9,1);
        %Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    
        y_guess = y;  
    end
    xvals3 = [xvals3 Xvec]; Avals3 = [Avals3 A*ones(size(Xvec))/2/pi]; Fvals3 = [Fvals3 Fvec]; 
    
    %% calculate shape for various widths (lower branch, pt 2)
    X = xcat - 0.01;
    solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 bcat 0.1 0 pi -10 2 0]); %m = 1
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    y_guess = deval(sol,t); y = y_guess;
    if A >= 1.2*2*pi    
    Xvec2 = [ 2*xcat:0.005:2]/2;
    muvec2 = zeros(1,length(Xvec2));
    Fvec2 = zeros(1,length(Xvec2));
    %Evec = zeros(1,length(Xvec));
    for j = 1:length(Xvec2)
        X = Xvec2(j);
        solinit = bvpinit(linspace(0,1,11),@newguess);
        sol = bvp4c(@odesystem,@bcs,solinit);
        t = linspace(0,1,N);
        y = deval(sol,t);
        r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
        mu = y(7,1); %b = r((N+1)/2);
        muvec(j) = mu;
        Fvec2(j) = -2*pi*y(9,1);
        %Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));

        y_guess = y;  
    end
    xvals3 = [xvals3 Xvec2]; Avals3 = [Avals3 A*ones(size(Xvec2))/2/pi]; Fvals3 = [Fvals3 Fvec2]; 
    end
    
    
end

%% case 3: 1 < Abar < 1.2 (three branches)
Abarvec = 1.00:0.01:1.2;% 0.1:0.1:1.3;
for Ael = Abarvec*2*pi
    A = Ael;
    Abar = Ael/2/pi
    
    bhans =  fsolve(@catarea,[0.1,1]);
    bcat = bhans(1)
    xcat = bhans(2)
    if Abar > 1 && Abar < 1.2 %in this regime, there are two catenoids
    bhans2 = fsolve(@catarea,[0.1,0.1]);
    bcat2 = bhans2(1)
    xcat2 = bhans2(2)
    end

    X = xcat - 0.01;
    solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 bcat 0.1 0 pi -10 2 0]); %m = 1
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    y_guess = deval(sol,t); y = y_guess;

    %% calculate shape for various widths (lower branch)
    Xvec = [ 2*xcat - 0.01:-0.025:0.0 0]/2;
    Fvec = zeros(1,length(Xvec));
    for j = 1:length(Xvec)
        X = Xvec(j);
        solinit = bvpinit(linspace(0,1,11),@newguess);
        sol = bvp4c(@odesystem,@bcs,solinit);
        t = linspace(0,1,N);
        y = deval(sol,t);
        r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
        mu = y(7,1); %b = r((N+1)/2);
        Fvec(j) = -2*pi*y(9,1);
    
        y_guess = y;  
    end
    if A < 1.2*2*pi %add catenoid data
        Fcat = 2*pi*catstability(bcat)*bcat;
        if A > 2*pi
            Fcat2 = 2*pi*catstability(bcat2)*bcat2;
        end %overwrite small b cases
        if abs(Abar - 1.01) < 0.001
             Fcat2 = 2*pi*166.712*bcat2;
        elseif abs(Abar - 1.02) < 0.001
             Fcat2 = 2*pi*73.6746*bcat2;
        elseif abs(Abar - 1.03) < 0.001
             Fcat2 = 2*pi*44.5071*bcat2;
        elseif abs(Abar - 1.04) < 0.001
             Fcat2 = 2*pi*30.7885*bcat2;
        elseif abs(Abar - 1.05) < 0.001
             Fcat2 = 2*pi*22.9503*bcat2;
        elseif abs(Abar - 1.06) < 0.001
             Fcat2 = 2*pi*17.9306*bcat2;
        elseif abs(Abar - 1.07) < 0.001
             Fcat2 = 2*pi*14.4596*bcat2;
        elseif abs(Abar - 1.08) < 0.001
             Fcat2 = 2*pi*11.921*bcat2;
        elseif abs(Abar - 1.09) < 0.001
             Fcat2 = 2*pi*9.98208*bcat2;
        elseif abs(Abar - 1.10) < 0.001
             Fcat2 = 2*pi*8.44836*bcat2;
        elseif abs(Abar - 1.11) < 0.001
             Fcat2 = 2*pi*7.19878*bcat2;
        elseif abs(Abar - 1.12) < 0.001
             Fcat2 = 2*pi*6.15405*bcat2;
        elseif abs(Abar - 1.13) < 0.001
             Fcat2 = 2*pi*5.25976*bcat2;
        end   
    elseif abs(A - 1.2*2*pi) < 0.001
        Fcat = 7.247;
    end
    xvals4 = [xvals4 Xvec xcat]; Avals4 = [Avals4 A*ones(size(Xvec))/2/pi A/2/pi]; Fvals4 = [Fvals4 Fvec Fcat]; 
    


    %% calculate shape for various widths (upper branch between catenoids)
    Xvec2 = []; Fvec2 = [];
    if Abar <= 1.13
        X = xcat - 0.01;
        solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 bcat/2 0.1 0 pi -10 2 0]); %m = 1
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;

        if A == 2*pi*1
            Xvec2 = [2*xcat - 0.01:-0.025:0.3]/2;
        else
            Xvec2 = [2*xcat - 0.01:-0.025:2*xcat2 + 0.01]/2;
        end
        Fvec2 = zeros(1,length(Xvec2));
        for j = 1:length(Xvec2)
            X = Xvec2(j);
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec2(j) = -2*pi*y(9,1);

            y_guess = y;
        end
    end
    if A == 2*pi*1 %manually fill in points when Abar = 1 (critical case)
        xvals5 = [xvals5 0 0.05 0.1 xcat Xvec2]; Fvals5 = [Fvals5 0 -0.0001 -0.001 Fcat Fvec2]; Avals5 = [Avals5 1 1 1 1 A*ones(size(Xvec2))/2/pi];
    else
%         xvals5 = [xvals5 Xvec2 xcat]; Avals5 = [Avals5 A*ones(size(Xvec2))/2/pi A/2/pi]; Fvals5 = [Fvals5 Fvec2 Fcat]; 
        xvals5 = [xvals5 Xvec2 xcat xcat2]; Avals5 = [Avals5 A*ones(size(Xvec2))/2/pi A/2/pi A/2/pi]; Fvals5 = [Fvals5 Fvec2 Fcat Fcat2]; 
    end
    
    %% tether branch
    Xvec3 = []; Fvec3 = [];
    if Abar > 1.05 %tether branch
        X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 20 2 -40]);
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        
        %% calculate shape for various widths (lower branch)
        Xvec3 = [2:-0.05:2*xcat2+0.01]/2;
        Fvec3 = zeros(1,length(Xvec3));
        for j = 1:length(Xvec3)
            X = Xvec3(j);
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3(j) = -2*pi*y(9,1);
            y_guess = y;  
        end
    elseif Abar == 1.05 
        X = 1;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.23 0.1 0 pi 600 3.7 -40]); %A = 1.05*2*pi
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        
        %% calculate shape for various widths (lower branch)
        Xvec3 = [0.9:-0.025:0.4];%Abar = 1.05
        Fvec3 = zeros(1,length(Xvec3));
        for j = 1:length(Xvec3)
            X = Xvec3(j);
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3(j) = -2*pi*y(9,1)
            y_guess = y;
        end
        X = 1;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.23 0.1 0 pi 600 3.7 -40]); %A = 1.05*2*pi
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        Xvec3b = [0.95:0.01:1];%Abar = 1.05
        Fvec3b = zeros(1,length(Xvec3b));
        for j = 1:length(Xvec3b)
            X = Xvec3b(j);
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3b(j) = -2*pi*y(9,1);
            y_guess = y;
        end
        Xvec3 = [Xvec3 Xvec3b]; Fvec3 = [Fvec3 Fvec3b];
    elseif Abar == 1.04
        X = 1;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.23 0.1 0 pi 600 3.7 -40]); %A = 1.05*2*pi
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        
        %% calculate shape for various widths (lower branch)
        Xvec3 = [0.8:-0.025:0.4];%Abar = 1.04
        Fvec3 = zeros(1,length(Xvec3));
        for j = 1:length(Xvec3)
            X = Xvec3(j);
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3(j) = -2*pi*y(9,1);
            y_guess = y;
        end
        X = 1;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.23 0.1 0 pi 600 3.7 -40]); %A = 1.05*2*pi
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        
        %% calculate shape for various widths (lower branch)
        Xvec3b = [0.85 0.9];%Abar = 1.04
        Fvec3b = zeros(1,length(Xvec3b));
        for j = 1:length(Xvec3b)
            X = Xvec3b(j);
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3b(j) = -2*pi*y(9,1);
            y_guess = y;
        end
        Xvec3 = [Xvec3 Xvec3b 1 0.96 0.93]; Fvec3 = [Fvec3 Fvec3b 2*pi*48.976 2*pi*46.932 2*pi*45.397];
    elseif Abar == 1.03
        X = 0.7;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.195 0.1 0 pi 300 2.0 -0]); %A = 1.05*2*pi
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        
        %% calculate shape for various widths (lower branch)
        Xvec3 = [0.7:-0.025:0.325 ];
        Fvec3 = zeros(1,length(Xvec3));
        for j = 1:length(Xvec3)
            X = Xvec3(j);
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3(j) = -2*pi*y(9,1);
            y_guess = y;
        end
        X = 0.7;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.195 0.1 0 pi 300 2.0 -0]); %A = 1.05*2*pi
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        
        %% calculate shape for various widths (lower branch)
        Xvec3b = [0.7 0.725 0.75];
        Fvec3b = zeros(1,length(Xvec3b));
        for j = 1:length(Xvec3b)
            X = Xvec3b(j)
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3b(j) = -2*pi*y(9,1);
            y_guess = y;
        end
        Xvec3 = [Xvec3 Xvec3b 1 0.91 0.85 0.8]; Fvec3 = [Fvec3 Fvec3b 2*pi*65.65 2*pi*59.547 2*pi*55.461 2*pi*52.051];
    elseif Abar == 1.02
        X = 0.6;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.195 0.1 0 pi 300 2.0 -0]); %A = 1.05*2*pi
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        
        %% calculate shape for various widths (lower branch)
        Xvec3 = [0.6:-0.025:0.275]
        Fvec3 = zeros(1,length(Xvec3));
        for j = 1:length(Xvec3)
            X = Xvec3(j)
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3(j) = -2*pi*y(9,1);
            y_guess = y;
        end
        X = 0.6;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.195 0.1 0 pi 300 2.0 -0]); %A = 1.05*2*pi
        sol = bvp4c(@odesystem,@bcs,solinit,options);
        y_guess = deval(sol,t); y = y_guess;
        
        %% calculate shape for various widths (lower branch)
        Xvec3b = [0.6:0.025:0.7];
        Fvec3b = zeros(1,length(Xvec3b));
        for j = 1:length(Xvec3b)
            X = Xvec3b(j)
            solinit = bvpinit(linspace(0,1,11),@newguess);
            sol = bvp4c(@odesystem,@bcs,solinit);
            t = linspace(0,1,N);
            y = deval(sol,t);
            r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
            mu = y(7,1); %b = r((N+1)/2);
            Fvec3b(j) = -2*pi*y(9,1);
            y_guess = y;
        end
        Xvec3 = [Xvec3 Xvec3b 1 0.91 0.8 0.75]; Fvec3 = [Fvec3 Fvec3b 2*pi*98.986 2*pi*89.887 2*pi*78.724 2*pi*73.644];
    elseif Abar == 1.01
       Xvec3 = [Xvec3  0.205 0.25 0.275 0.3 0.325 0.35 0.45 0.49]; 
       Fvec3 = [Fvec3 98.55 287 320.9 354.7 388 421 551.1 602.6];
       %1e2*[2.459586499368660   2.244610335905828   1.938658977640152   1.529812827481815 0.985492465677431]
    end
    if A > 2*pi
        xvals6 = [xvals6 Xvec3 xcat2]; Avals6 = [Avals6 A*ones(size(Xvec3))/2/pi A/2/pi]; Fvals6 = [Fvals6 Fvec3 Fcat2]; 
%         xvals6 = [xvals6 Xvec3 ]; Avals6 = [Avals6 A*ones(size(Xvec3))/2/pi ]; Fvals6 = [Fvals6 Fvec3 ]; 
    end
end

Abarsamp1 = 0.90:0.01:1.0; % case 1 (small Abar)
xcatsamp1 = zeros(size(Abarsamp1));
Fcatsamp1 = zeros(size(Abarsamp1));
xsamplen1 = 40;
sampmat1 = zeros(length(Abarsamp1),xsamplen1);
for n =1:length(Abarsamp1)
    if Abarsamp1(n) < 1.2
        A = Abarsamp1(n)*2*pi;
        bhans =  fsolve(@catarea,[0.1,1]); %xcatsamp(n) = bhans(2)
        sampmat1(n,:) = linspace(0,bhans(2),xsamplen1);
    else
        sampmat1(n,:) = linspace(0,1,xsamplen1);
    end
end
yq1 = sampmat1; xq1 = repmat(Abarsamp1,[xsamplen1,1])';

Abarsamp2 = 1.0:0.01:1.2; %case 3 (medium Abar)
xcatsamp2 = zeros(size(Abarsamp2));
Fcatsamp2 = zeros(size(Abarsamp2));
xsamplen2 = 40;
sampmat2 = zeros(length(Abarsamp2),xsamplen2);
for n =1:length(Abarsamp2)
    A = Abarsamp2(n)*2*pi;
    bhans =  fsolve(@catarea,[0.1,1]); %xcatsamp(n) = bhans(2)
    bhans2 = fsolve(@catarea,[0.1,0.1]);
    if A == 2*pi
        sampmat2(n,:) = linspace(0,bhans(2),xsamplen2);
    elseif abs(A - 1.2*2*pi) < 0.001
        sampmat2(n,:) = linspace(0,0.6627,xsamplen2);
    else
        sampmat2(n,:) = linspace(bhans2(2),bhans(2),xsamplen2);
    end
end
yq2 = sampmat2; xq2 = repmat(Abarsamp2,[xsamplen2,1])';

Abarsamp3 = 1.2:0.01:1.4; %case 2 (large Abar)
xcatsamp3 = zeros(size(Abarsamp3));
Fcatsamp3 = zeros(size(Abarsamp3));
xsamplen3 = 80;
sampmat3 = zeros(length(Abarsamp3),xsamplen3);
for n =1:length(Abarsamp3)
    if Abarsamp3(n) < 1.2
        A = Abarsamp3(n)*2*pi;
        bhans =  fsolve(@catarea,[0.1,1]); %xcatsamp(n) = bhans(2)
        sampmat3(n,:) = linspace(0,bhans(2),xsamplen3);
    else
        sampmat3(n,:) = linspace(0,1,xsamplen3);
    end
end
yq3 = sampmat3; xq3 = repmat(Abarsamp3,[xsamplen3,1])';

Abarsamp4 = Abarsamp2; %case 3 (medium Abar)
xcatsamp4 = zeros(size(Abarsamp4));
Fcatsamp4 = zeros(size(Abarsamp4));
xsamplen4 = 50;
sampmat4 = zeros(length(Abarsamp4),xsamplen4);
for n =1:length(Abarsamp4)
    A = Abarsamp4(n)*2*pi;
    bhans =  fsolve(@catarea,[0.1,0.1]); %xcatsamp4(n) = bhans(2)
    sampmat4(n,:) = linspace(bhans(2),1,xsamplen4);
    if abs(Abarsamp4(n) - 1.2) < 0.0001
        sampmat4(n,:) = linspace(0.6627,1,xsamplen4);
    end
end
yq4 = sampmat4; xq4 = repmat(Abarsamp4,[xsamplen4,1])';

Abarsamp5 = Abarsamp2; %case 3 (medium Abar)
xcatsamp5 = zeros(size(Abarsamp5));
Fcatsamp5 = zeros(size(Abarsamp5));
xsamplen5 = 40;
sampmat5 = zeros(length(Abarsamp5),xsamplen5);
for n =1:length(Abarsamp5)
    if Abarsamp5(n) < 1.1997
        A = Abarsamp5(n)*2*pi;
        bhans =  fsolve(@catarea,[0.1,1]); %xcatsamp(n) = bhans(2)
        sampmat5(n,:) = linspace(0,bhans(2),xsamplen5);
    elseif abs(Abarsamp5(n) - 1.2) < 0.0001
        sampmat5(n,:) = linspace(0,0.6627,xsamplen5);
    end
end
yq5 = sampmat5; xq5 = repmat(Abarsamp5,[xsamplen5,1])';

[xq,yq,vq] = griddata(Avals3',xvals3',Fvals3',xq3,yq3);
k = 0;
for j =1:length(xq(:,1))
    xvec = xq(j,:); yvec = yq(j,:); zvec = vq(j,:);
    if j== 11 || j == 1%Abar = 1.3
        plot3(xvec,yvec,zvec,'k-','linewidth',3);
    elseif mod(k,2) == 0
        plot3(xvec,yvec,zvec,'k-','linewidth',0.5);
    end
    k = k + 1;
end

[xq,yq,wq] = griddata(Avals',xvals',Fvals',xq1,yq1);
k = 0;
for j =1:length(xq(:,1))
    xvec = xq(j,:); yvec = yq(j,:); zvec = wq(j,:);
    if j== length(xq(:,1)) %Abar = 1
        plot3(xvec,yvec,zvec,'k-','linewidth',3);
    elseif mod(k,2) == 0
        plot3(xvec,yvec,zvec,'k-','linewidth',0.5);
    end
    k = k+1;
end

[xq,yq,uq] = griddata(Avals2',xvals2',Fvals2',xq1,yq1);
k = 0;
for j =1:length(xq(:,1))
    xvec = xq(j,:); yvec = yq(j,:); zvec = uq(j,:);
    if j== length(xq(:,1)) %Abar = 1.2
        plot3(xvec,yvec,zvec,'k-','linewidth',3);
    elseif mod(k,2) == 0
        plot3(xvec,yvec,zvec,'k-','linewidth',0.5);
    end
    k = k+1;
end

[xq,yq,tq] = griddata(Avals4',xvals4',Fvals4',xq5,yq5);
k = 0;
for j =1:length(xq(:,1))
    xvec = xq(j,:); yvec = yq(j,:); zvec = tq(j,:);
    if j== 11 || j == 1
        plot3(xvec,yvec,zvec,'k-','linewidth',3);
    elseif mod(k,2) == 0
        plot3(xvec,yvec,zvec,'k-','linewidth',0.5);
    end
    k = k+1;
end

[xq,yq,sq] = griddata(Avals5',xvals5',Fvals5',xq2,yq2);
k = 0;
for j =1:length(xq(:,1))
    xvec = xq(j,:); yvec = yq(j,:); zvec = sq(j,:);
    if j== 11
        plot3(xvec,yvec,zvec,'k-','linewidth',3);
    elseif mod(k,2) == 0
        plot3(xvec,yvec,zvec,'k-','linewidth',0.5);
    end
    k = k+1;
end

[xq,yq,rq] = griddata(Avals6',xvals6',Fvals6',xq4,yq4);
k = 0;
for j =1:length(xq(:,1))
    xvec = xq(j,:); yvec = yq(j,:); zvec = rq(j,:);
    if j== 11|| j == 1
        plot3(xvec,yvec,zvec,'k-','linewidth',3);
    elseif mod(k,2) == 0
        plot3(xvec,yvec,zvec,'k-','linewidth',0.5);
    end
    k = k+1;
end


s1 = surf(xq3,yq3,vq,'facealpha',0.8); hold on; grid on; box on;
s2 = surf(xq1,yq1,wq,'facealpha',0.8);
s3 = surf(xq1,yq1,uq,'facealpha',0.8);
s4 = surf(xq5,yq5,tq,'facealpha',0.8);
s5 = surf(xq2,yq2,sq,'facealpha',0.8);
s6 = surf(xq4,yq4,rq,'facealpha',0.8);
s7 = surf([1 1.01; 1 1.01],[0 0.2054; 0 0.2216],[0 60.78; 250 143],'facealpha',0.8); s7.FaceColor = 'interp';
s8 = surf([1 1.01; 1 1.01],[0 0.2216; 0 0.2378],[250 143; 400 225.1],'facealpha',0.8); s8.FaceColor = 'interp';
s9 = surf([1.19 1.19; 1.19 1.2],[0.6547 0.6547; 0.6531 0.6627],[-4.1 -4.1; 3.687 6.762],'facealpha',0.8); s9.FaceColor = 'interp';
s1.EdgeColor = 'none'; s2.EdgeColor = 'none'; 
s3.EdgeColor = 'none'; s4.EdgeColor = 'none';
s5.EdgeColor = 'none'; s6.EdgeColor = 'none';
s7.EdgeColor = 'none'; s8.EdgeColor = 'none';
s9.EdgeColor = 'none';
s1.FaceColor = 'interp'; s2.FaceColor = 'interp';
s3.FaceColor = 'interp'; s4.FaceColor = 'interp';
s5.FaceColor = 'interp'; s6.FaceColor = 'interp';
s7.FaceColor = 'interp'; s8.FaceColor = 'interp';
s9.FaceColor = 'interp';
xlabel('$\bar{A}$','interpreter','latex'); zlabel('$Fa/\kappa$','interpreter','latex'); ylabel('$h/a$','interpreter','latex');
set(gca, 'YDir', 'reverse');
set(gcf,'color','w');
set(gca,'fontsize',18);
zlim([-230 130]); caxis([-230 180]); colormap(flipud(cool));



end

%% nonzero spontaneous curvature
function [] = scsolutions()
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
global c0; c0 = 1; %spontaneous curvature
 
%% prepare a solution as initial guess
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.200001 0.1 0 pi 0 2 0]); %m = 5b
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %m = 5
% solinit = bvpinit(linspace(0,1,5),@guess4); %m = 4
% % solinit = bvpinit(linspace(0,1,11),@guess3); %m = 3
% solinit = bvpinit(linspace(0,1,9),@guess3b); %m = 3b
% solinit = bvpinit(linspace(0,1,11),@guess2); %m = 2
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
X = 0.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
gam = y(5,:);
mu = y(7,1)
L = y(8,1)
eta = y(9,1)

% plot(z,r); hold on;

sol = bvp4c(@odesystem_sc,@bcs_sc,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r2 = y(3,:); z2 = y(4,:);
% plot(z2,r2);

Xvec = [ 1.0 1.05 1.055 0.5276*2 0.52765*2]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
evvec = zeros(1,length(Xvec));
hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem_sc,@bcs_sc,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1); %b = r((N+1)/2);
    muvec(j) = mu;
     %k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r - c0).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
%     if mod(j,2) == 1
         plot(z,r,'b','LineWidth',2);
%     end
    y_guess = y;
    
%     %% stability analysis
%     [eval,efun]= evmatrix(sol);
%     evvec(j) = eval;
%     if mod(j,4) == 1
%         if eval > 0
%         subplot(1,2,1)
%         plot(z,r,'b-','linewidth',2); hold on;
%         subplot(1,2,2)
%         plot(t,efun,'b-','linewidth',2); hold on;
%         else
%         subplot(1,2,1)
%         plot(z,r,'b-.','linewidth',2); hold on;
%         subplot(1,2,2)
%         plot(t,efun,'b-.','linewidth',2); hold on;
%         end
%     end
end
Fvec

end

%% shows no maximal extension for A > Ac
function [] = tethertransition()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A =  7.59; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess (A = 7.59)
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3, upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3, upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi 0 2 0]); %<--- n=1 upper branch
% solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1 upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]);%n=3 lower branch
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]); %n=1 lower branch (A=7.59)
% X = 1.34/2; 
%solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.3 0.1 0 pi -13 2.22 2.84]); %n=1 lower branch (A=10)
% solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.162 0.1 0 pi 0 1.2 2.84]); 
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.6 0]); %5th order, upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.11 0.1 0 pi -200 1.5 200]); %7th order, lower branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.16 0.1 0 pi -200 1.5 200]); %11th order, lower branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.16 0.1 0 pi -300 1.5 300]); %9th order, lower branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.17 0.1 0 pi -300 1.5 300]); %7th order, upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.10 0.1 0 pi -120 1.5 300]); %9th order, upper branch
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.17 0.1 0 pi -120 1.5 200]); %11th order, upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.165 0.1 0 pi -120 1.5 300]); %5th order, special
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.2 0.18001 0.1 0 pi -80 1.2 100]); %7th order, special
% solinit = bvpinit(linspace(0,1,11),[0.1 0.2 0.35 0.1 0 pi 40 1.2 -50]); %???
% solinit = bvpinit(linspace(0,1,9),@guess2_tether2); 
% solinit = bvpinit(linspace(0,1,6),@guess2_tether6); 
% solinit = bvpinit(linspace(0,1,5),@guess4_tether3); 
% solinit = bvpinit(linspace(0,1,5),@guess2_tether3); %2nd order <----
% solinit = bvpinit(linspace(0,1,5),@guess4_tether2); 
% solinit = bvpinit(linspace(0,1,6),@guess4_tether2); 
% X = 0.65; solinit = bvpinit(linspace(0,1,6),@guess4_tether4); 
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.15 0.1 0 pi 100 2 -100]) ;

options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); H = (psi_s + sin(psi)./r)/2;
    y(7,1)
    y(8,1)
    y(9,1)
    indices = find([0 diff(sign(H))]~=0);
    
% plot(t,H); return;
plot(t,r); hold on; plot(z(indices),r(indices),'o'); return

Xvec0 = [1:0.01:1.33 1.333 1.3333]/2;
%Xvec0 = [1.33:-0.01:1.05]/2;
%Xvec0 = [1.0:-0.01:0.6]/2;%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
muvec0 = zeros(1,length(Xvec0));
Fvec0 = zeros(1,length(Xvec0));
figure; hold on; axis equal; box on;
for j = 1:length(Xvec0)
    X = Xvec0(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec0(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec0(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2; K = psi_s.*sin(psi)./r;
    if mod(j,5) == 1
    %plot(z,r,'k','linewidth',2);
    end
   
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end

disp('branch 1')
X = 1.04/2;
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3, upper branch
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec = [1.0:0.01:1.33 1.333]/2;
%Xvec = [1.0:-0.01:0]/2;%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
%Xvec = [1.4:0.1:3]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    %plot(z,r,'b','LineWidth',2);
    end
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end

disp('branch 2')
X = 1.04/2;
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]);%n=3 lower branch
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec2 = [1.0:0.01:1.33 1.333]/2;
%Xvec2 = [1.0:-0.01:0]/2;%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    %plot(z,r,'r','linewidth',2);
    end
   
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end

disp('branch 3')
X = 1.04/2; 
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]); %n=1 lower branch (A=7.59)
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%Xvec3 = [1.0:-0.01:0]/2;
Xvec3 = [1.0:0.01:1.33 1.333]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec3(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    %plot(z,r,'g','linewidth',2);
    end
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end

disp('branch 4')
X = 1.04/2; 
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.6 0]); %5th order, upper branch
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%Xvec4 = [1.0:-0.01:0]/2;
Xvec4 = [1.0:0.01:1.33 1.332]/2;
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec4(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    %plot(z,r,'r','linewidth',2);
    end

    y_guess = y;
end

disp('branch 5')
X = 1.04/2; 
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%Xvec5 = [1.0]/2;
Xvec5 = [1.0:0.01:1.33 1.332]/2;
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec5)
    X = Xvec5(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec5(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec5(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    %plot(z,r,'r','linewidth',2);
    end

    y_guess = y;
end

disp('branch 6')
X = 1.04/2; 
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.165 0.1 0 pi -120 1.5 300]); %5th order, special
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec6 = [1.0:0.01:1.33 1.332]/2;
%Xvec6 = [1.0:-0.01:0]/2;
muvec6 = zeros(1,length(Xvec6));
Fvec6 = zeros(1,length(Xvec6));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec6)
    X = Xvec6(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec6(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec6(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    %plot(z,r,'r','linewidth',2);
    end
 
    y_guess = y;
end

disp('branch 7')
X = 1.04/2; 
solinit = bvpinit(linspace(0,1,9),@guess2_tether2); 
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec7 = [1.0:0.01:1.33 1.332 1.3333]/2;
%Xvec7 = [1.0:-0.01:0]/2;
muvec7 = zeros(1,length(Xvec7));
Fvec7 = zeros(1,length(Xvec7));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec7)
    X = Xvec7(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec7(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec7(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2; K = psi_s.*sin(psi)./r;
    if mod(j,5) == 1
%     plot(fliplr(z),r,'m','linewidth',2);
    end

    y_guess = y;
end

disp('branch 8')
X = 0.65; solinit = bvpinit(linspace(0,1,6),@guess4_tether4); 
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec8 = [1.333 1.33:-0.01:1.05]/2;
%Xvec8 = [1.0:-0.01:0]/2;
muvec8 = zeros(1,length(Xvec8));
Fvec8 = zeros(1,length(Xvec8));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec8)
    X = Xvec8(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec8(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec8(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    %plot(z,r,'m','linewidth',2);
    end

    y_guess = y;
end

disp('branch 9')
X = 0.5; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi 0 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec9 = [1:0.01:1.4]/2;
%Xvec9 = [1.0:-0.01:0]/2;
muvec9 = zeros(1,length(Xvec9));
Fvec9 = zeros(1,length(Xvec9));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec9)
    X = Xvec9(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec9(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec9(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    %plot(z,r,'m','linewidth',2);
    end

    y_guess = y;
end

disp('branch 10')
X = 0.52; solinit = bvpinit(linspace(0,1,6),@guess4_tether2);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec10 = [1:0.01:1.33 1.333]/2;
%Xvec10 = [1.0:-0.01:0]/2;
muvec10 = zeros(1,length(Xvec10));
Fvec10 = zeros(1,length(Xvec10));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec10)
    X = Xvec10(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec10(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec10(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
    plot(z,r,'m','linewidth',2);
    end

    y_guess = y;
end


%% catenoid solution
% bcat = 0.825517; v = linspace(-0.527697,0.527697);
% plot(v,bcat*cosh(v/bcat),'k')
% v = linspace(-0.670909,0.670909,1001);
% bcat = 0.559241;
% plot(v,bcat*cosh(v/bcat),'k--','linewidth',2)

xlabel('z'); ylabel('r');
title(['Tethers: A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure();
subplot(1,2,1)
plot(Xvec0,Fvec0,'o-','LineWidth',2,'color',[0.8500, 0.3250, 0.0980]); hold on;%%
plot(Xvec,Fvec,'o-','LineWidth',2,'color',[0.929 0.694 0.125]); hold on;%
plot(Xvec2,Fvec2,'o-','LineWidth',2,'color',[0.929 0.694 0.125]); hold on;%
plot(Xvec3,Fvec3,'o-','LineWidth',2,'color',[0.9290, 0.6940, 0.1250]); hold on;%
plot(Xvec4,Fvec4,'o-','LineWidth',2,'color',[0.4660, 0.6740, 0.1880]); hold on;
plot(Xvec5,Fvec5,'o-','LineWidth',2,'color',[0.4660, 0.6740, 0.1880]); hold on;
plot(Xvec6,Fvec6,'o-','LineWidth',2,'color',[0.4660, 0.6740, 0.1880]); hold on;
plot(Xvec7,Fvec7,'o-','LineWidth',2,'color',[0.850 0.325 0.098]); hold on;%%
plot(Xvec8,Fvec8,'o-','LineWidth',2,'color',[0.9290, 0.6940, 0.1250]); hold on;%
plot(Xvec9,Fvec9,'o-','LineWidth',2,'color',[0, 0.4470, 0.7410]); hold on;
plot(Xvec10,Fvec10,'o-','LineWidth',2,'color',[0.4940, 0.1840, 0.5560]); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18); 

subplot(1,2,2)
plot(Xvec0,muvec0,'o-','LineWidth',2,'color',[0.8500, 0.3250, 0.0980]); hold on;
plot(Xvec,muvec,'o-','LineWidth',2,'color',[0.929 0.694 0.125]); hold on;
plot(Xvec2,muvec2,'o-','LineWidth',2,'color',[0.929 0.694 0.125]); hold on;
plot(Xvec3,muvec3,'o-','LineWidth',2,'color',[0.9290, 0.6940, 0.1250]); hold on;
plot(Xvec4,muvec4,'o-','LineWidth',2,'color',[0.4660, 0.6740, 0.1880]); hold on;
plot(Xvec5,muvec5,'o-','LineWidth',2,'color',[0.4660, 0.6740, 0.1880]); hold on;
plot(Xvec6,muvec6,'o-','LineWidth',2,'color',[0.4660, 0.6740, 0.1880]); hold on;
plot(Xvec7,muvec7,'o-','LineWidth',2,'color',[0.850 0.325 0.098]); hold on;
plot(Xvec8,muvec8,'o-','LineWidth',2,'color',[0.9290, 0.6940, 0.1250]); hold on;
plot(Xvec9,muvec9,'o-','LineWidth',2,'color',[0, 0.4470, 0.7410]); hold on;
plot(Xvec10,muvec10,'o-','LineWidth',2,'color',[0.4940, 0.1840, 0.5560]); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');
%plot(Xvec2,Fvec2)

end

function [] = tethertransition2()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A =  10; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess (A = 7.59)
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3, upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi 0 2 0]); %<--- n=1 upper branch
% solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1 upper branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.11 0.1 0 pi 0 2 0]);%n=7 lower branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]);%n=5 lower branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 2 0]);%n=7 upper branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 2 0]);%n=7 upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]); %n=1 lower branch (A=7.59)
%   X = 1.4/2;
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.3 0.1 0 pi -13 2.22 2.84]); %n=3
X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.09 0.1 0 pi -13 2.22 2.84]); %n=5L
X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.12 0.1 0 pi -13 2.22 2.84]); %n=5U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.14 0.1 0 pi -13 2.22 2.84]); %n=3U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.22 0.1 0 pi -13 2.22 2.84]); %n=3L
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.1 0 pi -13 2.22 2.84]); %n=7L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.12 0.1 0 pi -13 2.22 2.84]); %n=7U
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.13 0.11 0 pi -13 2.22 2.84]); %n=3T
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 4 0.08 0.11 0 pi -13 2.22 2.84]); %n=9L
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.15 0.1 0 pi -13 2.22 2.84]); %n=3U
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.22 0.11 0 pi -13 2.22 2.84]); %n=3S
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.15 0.11 0 pi -3 3.52 1.84]); %n=3S
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.21 0.11 0 pi -3 3.52 1.84]); %n=5T
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.28 0.11 0 pi -7.4 3.54 2.51]); %n=??
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.27 0.11 0 pi -7.4 3.54 2.51]); %n=??
X = 1.66/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.21 0.11 0 pi -3 3.52 1.84]); %n=5S
%X = 1.7/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.11 0 pi -13 2.22 2.84]); %n=9L
X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -13 2.22 2.84]); %n=5T
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.25 0.11 0 pi -13 2.22 2.84]); %n=3L
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.23 0.11 0 pi -13 2.22 2.84]); %n=3S

%% prepare a solution as initial guess (A = 2*pi*2)
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3, upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1 upper branch

options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    K = psi_s.*sin(psi)./r;
    y(7,1)
    y(8,1)
    y(9,1)
    indices = find([0 diff(sign(H))]~=0);
    
% plot(z,r); hold on; grid on; plot(z(indices),r(indices),'o'); return
% plot(t,H);  hold on; plot(t,-H);  grid on; return;
  
%% calculate shape for various widths
Xvec = [1.8:0.01:1.88]/2;
%Xvec = [1.5:0.01:1.68 1.6893]/2;
%Xvec = [1.64:0.01:1.68 1.6893 1.7:0.01:1.85 1.8525 1.855]/2;%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
%Xvec = [1.4:0.1:3]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
figure; hold on;  box on; grid on; %axis equal;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,2) == 1
        plot(z,r,'b','LineWidth',2);
    end
%     if X == Xvec(end)
%         m = sqrt(abs(mu));
%         s =L*(t-0.5);
%         d = sqrt(2/abs(mu));
%         b  = d;
%         %plot(t,r,'LineWidth',2); hold on;
%         vmax = fzero(@(z) 1 - b*cosh(z/b),1);
%         v= linspace(-vmax,vmax);
%         %plot(v-z(end)-v(end),b*cosh(v/b))        
%     end
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end


%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.1 0 pi -13 2.22 2.84]); %n=7L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.09 0.1 0 pi -13 2.22 2.84]); %n=5L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.22 0.1 0 pi -13 2.22 2.84]); %n=3L
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.25 0.11 0 pi -13 2.22 2.84]); %n=3S
X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.14 0.1 0 pi -13 2.22 2.84]); %n=3U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.12 0.1 0 pi -13 2.22 2.84]); %n=5U
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.1 0 pi -13 2.22 2.84]); %n=7L
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec2 = [1.5:0.01:1.68 1.6893]/2;
%Xvec2 = [1.8:0.01:1.85 1.855]/2;
%Xvec2 = [1.5:0.01:1.85 1.851 1.852 1.853 1.8533 1.85335]/2;%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,2) == 1
      plot(z,r,'r','linewidth',2);
    end
    if X == Xvec(end)
        mu
        plot(z,r,'k','LineWidth',2) ;
    end
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end
 
%% catenoid solution
% bcat = 0.825517; v = linspace(-0.527697,0.527697);
% plot(v,bcat*cosh(v/bcat),'k')
% v = linspace(-0.670909,0.670909,1001);
% bcat = 0.559241;
% plot(v,bcat*cosh(v/bcat),'k--','linewidth',2)

xlabel('z'); ylabel('r');
title(['Tethers: A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure();
subplot(1,2,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18); 

subplot(1,2,2)
plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');
%plot(Xvec2,Fvec2)

end

function [] = tethertransition3()
 
global kappa; kappa = 1;
global kbar; kbar = 2;
global A; A =  2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
global c0; c0 = 0;


%% prepare a solution as initial guess (A = 2*pi*4)
% X = 1.5; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.12 0.1 0 pi 0 2 0]);%n=5U
% X = 2; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.12 0.1 0 pi 0 2 0]);%n=5U
% X = 2.5; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.15 0.1 0 pi 0 5.3 0]);%n=11U
% X = 3; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.1 0.1 0 pi 0 4 0]);%n=11T
% X = 3; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.16 0.1 0 pi -100 4 20]);%n=13T
% X = 3; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.17 0.1 0 pi -10 5 20]);%n=17T
% X = 1.5; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.16 0.1 0 pi -10 5 20]);%n=13U
% X = 1.5; solinit = bvpinit(linspace(0,1,11),[pi/2 1 0.12 0.1 0 pi -10 4 20]);%n=13L
% X = 1.5; solinit = bvpinit(linspace(0,1,11),[pi/2 1 0.15 0.1 0 pi -10 4 20]);%n=7L
% X = 1.2; solinit = bvpinit(linspace(0,1,11),[pi/2 1 0.18 0.1 0 pi -10 4 20]);%n=5L
% X = 1.2; solinit = bvpinit(linspace(0,1,11),[pi/2 1 0.17 0.1 0 pi -1 4 2]);%n=7U
% X = 1.2; solinit = bvpinit(linspace(0,1,11),[pi/2 1 0.19 0.1 0 pi -1 4 2]);%n=5U
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 1 0.11 0.1 0 pi -1 4 2]);%n=9U
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 1 0.19 0.1 0 pi -1 4 2]);%n=9U
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.1 0.2 0.1 0 pi -1 5 -1]);%n=3U
% X = 0.7; solinit = bvpinit(linspace(0,1,5),@guess2_tether5);

%% prepare a solution as initial guess (A = 2*pi*3)
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 -2 0.29 0.1 0 pi 0 2.6 0]);%n=3T
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.17 0.1 0 pi 0 2.6 0]);%n=5U
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi 0 2.6 0]);%n=3U
% % X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 3 0.11 0.1 0 pi -1 2.8 0]);%n=9U
% X = 1.2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -1 2.8 0]);%n=3L


%% prepare a solution as initial guess (A = 2*pi*2)
%  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.18 0.1 0 pi -3 2.7 5]); %n=3L
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.26 0.1 0 pi -3 2.7 5]); %n=3T
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.34 0.1 0 pi -3 2.7 5]); %n=2
%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1U
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]);%n=5L
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.11 0.1 0 pi 0 2 0]);%n=3L
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.09 0.1 0 pi 0 2 0]);%n=5T
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.16 0.1 0 pi 0 2 0]);%n=5S
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 -2 0.29 0.1 0 pi 0 2.6 0]);%n=3T

%% prepare a solution as initial guess (A = 2*pi*1.4)
% X = 0.7; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.18 0.1 0 pi 0 1.7 0]); %n=3U
% X = 0.79; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.09 0.1 0 pi -5 2.5 -10]); %n=9L
% X = 0.79; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.11 0.1 0 pi -5 2.5 -10]); %n=7L
% X = 0.79; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.2 0.1 0 pi -3 2.5 -5]); %n=3S

%% prepare a solution as initial guess (A = 2*pi*1.5)
% X = 0.8; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.18 0.1 0 pi 0 1.7 0]); %n=3U
% X = 0.8; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.2 0.1 0 pi -3 2.5 -5]); %n=3S

%% prepare a solution as initial guess (A = 2*pi*1.6)
% X = 0.8; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.18 0.1 0 pi 0 1.7 0]); %n=3U
% X = 0.81; solinit = bvpinit(linspace(0,1,11),[pi/2 2 Xg 0.1 0 pi -3 2.5 -5]); %n=3U
% X = 0.79; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.09 0.1 0 pi -5 2.5 -10]); %n=9L
% X = 0.79; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.11 0.1 0 pi -5 2.5 -10]); %n=7L

%% prepare a solution as initial guess (A = 2*pi*1.3)
X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess4_tether_kbar);
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether_kbar2);
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5);
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6;
%  X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether5);
%  X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether6);


options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); H = (psi_s + sin(psi)./r)/2;
    y(7,1)
    y(8,1)
    y(9,1)
    indices = find([0 diff(sign(H))]~=0);
    
%  plot(z,r); hold on; plot(z(indices),r(indices),'o'); grid on; return
% plot(t,H); return

%% calculate shape for various widths
options = bvpset('RelTol',1e-8);
Xvec = [0.6:0.01:0.69 0.695 0.6955 0.6958];
% Xvec = [2:0.01:2.3];%Xvec = [1.4:0.1:3]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure; hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,2) == 1
        plot(t,H,'b','LineWidth',2);
    end
    if X == Xvec(end)
        plot(t,H,'k','LineWidth',2);
    end
%     if X == Xvec(end)
%         m = sqrt(abs(mu));
%         s =L*(t-0.5);
%         d = sqrt(2/abs(mu));
%         b  = d;
%         plot(z,r,'LineWidth',2); hold on;
%         vmax = fzero(@(z) 1 - b*cosh(z/b),1);
%         v= linspace(-vmax,vmax);
%         %plot(v-z(end)-v(end),b*cosh(v/b))        
%     end
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end


% X = 1.2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -1 2.8 0]);%n=3L
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.11 0.1 0 pi 0 2 0]);%n=3L
% X = 0.79; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.2 0.1 0 pi -3 2.5 -5]); %n=3S
% X = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]);%n=5L
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6);
X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess4_tether_kbar2);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec2 = [0.6:0.01:0.7 0.701:0.0001:0.7027 0.70275 0.7028];
%Xvec2 = [1.2:0.01:1.5 1.49:-0.01:1.34];%[4.23:-0.05:1.10 1.05:-0.05:0.10]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
for j = 1:length(Xvec2)
    X = Xvec2(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    plot(z,r,'r','linewidth',2);
%     if X == Xvec(end)
%         mu
%         s =L*(t-0.5);
%         d = sqrt(mu);
%         sqrt(mu/2)
%         %plot(t-0.5,-1/2/b*exp(-s/d).*(-exp(L/d) + exp(2*s/d))); hold on;
%         plot(t-0.5,H,'LineWidth',2) ;
%     end
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end
 
%% catenoid solution
% bcat = 0.825517; v = linspace(-0.527697,0.527697);
% plot(v,bcat*cosh(v/bcat),'k')
% v = linspace(-0.670909,0.670909,1001);
% bcat = 0.559241;
% plot(v,bcat*cosh(v/bcat),'k--','linewidth',2)

xlabel('z'); ylabel('r');
title(['Tethers: A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure();
subplot(1,3,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18); 

subplot(1,3,2)
plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');


end

function [] = tethertransition4_mode1()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A =  2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 0 1.7 0]);
X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether7); %kink
% X = 0.46; solinit = bvpinit(linspace(0,1,6),@guess4_tether8); %double kink
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    indices = find([0 diff(sign(H))]~=0);
    y(7,1)
    y(8,1)
    y(9,1)
    
% plot(z,r); hold on; grid on; plot(z(indices),r(indices),'o'); return
% plot(t,H);  hold on;  grid on; return;
  
%% ========== colors ============
%[0, 0.4470, 0.7410]
%[0.8500, 0.3250, 0.0980]
%[0.9290, 0.6940, 0.1250] 	
%[0.4940, 0.1840, 0.5560]
%[0.4660, 0.6740, 0.1880]

%% ------------------ 1st branch --------------------
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
%Xvec = [1.04:0.01:3]/2;
 Xvec = [ 0.58:-0.01:0.12 0.115 0.11];
%  Xvec = [0.45 0.445 0.435 0.40:-0.01:0.29 0.271:-0.01:-0.11 -0.15 -0.25 -0.3 -0.35 -0.4 -0.45 -0.475...
%      -0.45 -0.4 -0.35 -0.25 -0.21 -0.15 -0.1 -0.05 0 0.1];
%  Xvec = [0.45 0.445 0.435 0.40:-0.01:0.29 0.271:-0.01:-0.11 0];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
        plot(z,r,'LineWidth',2);
    end
    if X == Xvec(end)
        plot(z,r,'k','LineWidth',2);
    end

    y_guess = y;
end

%% ------------------ 2nd branch --------------------
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 0 1.7 0]);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
%Xvec = [1.04:0.01:3]/2;
 Xvec2 = [ 0.11:0.02:1];
%  Xvec = [0.45 0.445 0.435 0.40:-0.01:0.29 0.271:-0.01:-0.11 -0.15 -0.25 -0.3 -0.35 -0.4 -0.45 -0.475...
%      -0.45 -0.4 -0.35 -0.25 -0.21 -0.15 -0.1 -0.05 0 0.1];
%  Xvec = [0.45 0.445 0.435 0.40:-0.01:0.29 0.271:-0.01:-0.11 0];
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
        %plot(z,r,'r','LineWidth',2);
    end


    y_guess = y;
end



figure(1)
xlabel('z'); ylabel('r'); axis equal; box on;
title(['1st mode, A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure();
subplot(1,3,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Fvec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,Fvec4,'mo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18); 

subplot(1,3,2)
plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,muvec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,muvec4,'mo-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Evec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,Evec4,'mo-','LineWidth',2); hold on;
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');


end

function [] = tethertransition4_mode2()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A =  2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess (A = 7.59)
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.7 0]); %n=1U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
%  X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.07 0.1 0 pi 0 1.7 0]); %n=3S
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.05 0.1 0 pi 0 1.7 0]); %n=3L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.95 0.1 0 pi 0 1.7 0]); %n=5L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.081 0.11 0 pi -3 1.9 1.4]); %n=7L
% X = 1.24/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.1 0.11 0 pi -3 1.9 1.4]); %n=3T
% X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.7 0]); %n=5S
% X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.14 0.1 0 pi 0 1.7 0]); %n=3L
% X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -1 1.7 1]); %n=3S

 %X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.25 0.11 0 pi -1 1.7 1]);

%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi 0 2 0]); %<--- n=1 upper branch
% solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1 upper branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.11 0.1 0 pi 0 2 0]);%n=7 lower branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]);%n=5 lower branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 2 0]);%n=7 upper branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 2 0]);%n=7 upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]); %n=1 lower branch (A=7.59)
%   X = 1.4/2;
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.3 0.1 0 pi -13 2.22 2.84]); %n=3
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.09 0.1 0 pi -13 2.22 2.84]); %n=5L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.12 0.1 0 pi -13 2.22 2.84]); %n=5U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.14 0.1 0 pi -13 2.22 2.84]); %n=3U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.22 0.1 0 pi -13 2.22 2.84]); %n=3L
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.1 0 pi -13 2.22 2.84]); %n=7L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.12 0.1 0 pi -13 2.22 2.84]); %n=7U
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.13 0.11 0 pi -13 2.22 2.84]); %n=3T
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 4 0.08 0.11 0 pi -13 2.22 2.84]); %n=9L
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.15 0.1 0 pi -13 2.22 2.84]); %n=3U
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.22 0.11 0 pi -13 2.22 2.84]); %n=3S
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.15 0.11 0 pi -3 3.52 1.84]); %n=3S
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.21 0.11 0 pi -3 3.52 1.84]); %n=5T
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.28 0.11 0 pi -7.4 3.54 2.51]); %n=??
% X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.27 0.11 0 pi -7.4 3.54 2.51]); %n=??
%X = 1.66/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.21 0.11 0 pi -3 3.52 1.84]); %n=5S
%X = 1.7/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.11 0 pi -13 2.22 2.84]); %n=9L
%X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -13 2.22 2.84]); %n=5T
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.25 0.11 0 pi -13 2.22 2.84]); %n=3L
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.23 0.11 0 pi -13 2.22 2.84]); %n=3S
%  X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether5); %n=4
% X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether6); %n =4
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5); %n=2
X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6); %n=2

options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    y(7,1)
    y(8,1)
    y(9,1)
    indices = find([0 diff(sign(H))]~=0);
    
% plot(z,r); hold on; grid on; plot(-z,r); plot(z(indices),r(indices),'o'); return
% plot(t,H);  hold on; grid on; return;
  
%% ========== colors ============
%[0, 0.4470, 0.7410]
%[0.8500, 0.3250, 0.0980]
%[0.9290, 0.6940, 0.1250] 	
%[0.4940, 0.1840, 0.5560]
%[0.4660, 0.6740, 0.1880]

%% ------------------ 1st branch --------------------
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
Xvec = [0.6:0.01:0.7 0.71 0.712];
% Xvec = [1.04:-0.01:0]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
        plot(-z,r,'b','LineWidth',2);
    end
    if X == Xvec(end)
        plot(-z,r,'k','LineWidth',2);
    end
%     if X == Xvec(end)
%         m = sqrt(abs(mu));
%         s =L*(t-0.5);
%         d = sqrt(2/abs(mu));
%         b  = d;
%         %plot(t,r,'LineWidth',2); hold on;
%         vmax = fzero(@(z) 1 - b*cosh(z/b),1);
%         v= linspace(-vmax,vmax);
%         %plot(v-z(end)-v(end),b*cosh(v/b))        
%     end
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end

%% ------------------ 2nd branch --------------------
X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5); %n=2
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
figure(1)
Xvec2 = [0.6:0.01:0.71 0.712];
% Xvec2 = [1.04:-0.01:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
      plot((z),r,'r','linewidth',2);
    end
    if X == Xvec2(end)
       plot((z),r,'k','linewidth',2); 
    end
    y_guess = y;
end
 

figure(1)
xlabel('z'); ylabel('r'); axis equal; box on;
title(['2nd mode, A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure();
subplot(1,3,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Fvec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,Fvec4,'mo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18); 

subplot(1,3,2)
plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,muvec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,muvec4,'mo-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Evec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,Evec4,'mo-','LineWidth',2); hold on;
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');


end

function [] = tethertransition4_mode3()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A =  2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess (A = 7.59)
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.7 0]); %n=1U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
%  X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.07 0.1 0 pi 0 1.7 0]); %n=3S
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.05 0.1 0 pi 0 1.7 0]); %n=3L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.95 0.1 0 pi 0 1.7 0]); %n=5L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.081 0.11 0 pi -3 1.9 1.4]); %n=7L
% X = 1.24/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.1 0.11 0 pi -3 1.9 1.4]); %n=3T
% X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.7 0]); %n=5S
% X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.14 0.1 0 pi 0 1.7 0]); %n=3L
% X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -1 1.7 1]); %n=3S

 %X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.25 0.11 0 pi -1 1.7 1]);

%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi 0 2 0]); %<--- n=1 upper branch
% solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1 upper branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.11 0.1 0 pi 0 2 0]);%n=7 lower branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]);%n=5 lower branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 2 0]);%n=7 upper branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 2 0]);%n=7 upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]); %n=1 lower branch (A=7.59)
%   X = 1.4/2;
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.3 0.1 0 pi -13 2.22 2.84]); %n=3
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.09 0.1 0 pi -13 2.22 2.84]); %n=5L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.12 0.1 0 pi -13 2.22 2.84]); %n=5U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.14 0.1 0 pi -13 2.22 2.84]); %n=3U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.22 0.1 0 pi -13 2.22 2.84]); %n=3L
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.1 0 pi -13 2.22 2.84]); %n=7L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.12 0.1 0 pi -13 2.22 2.84]); %n=7U
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.13 0.11 0 pi -13 2.22 2.84]); %n=3T
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 4 0.08 0.11 0 pi -13 2.22 2.84]); %n=9L
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.15 0.1 0 pi -13 2.22 2.84]); %n=3U
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.22 0.11 0 pi -13 2.22 2.84]); %n=3S
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.15 0.11 0 pi -3 3.52 1.84]); %n=3S
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.21 0.11 0 pi -3 3.52 1.84]); %n=5T
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.28 0.11 0 pi -7.4 3.54 2.51]); %n=??
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.27 0.11 0 pi -7.4 3.54 2.51]); %n=??
%X = 1.66/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.21 0.11 0 pi -3 3.52 1.84]); %n=5S
%X = 1.7/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.11 0 pi -13 2.22 2.84]); %n=9L
%X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -13 2.22 2.84]); %n=5T
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.25 0.11 0 pi -13 2.22 2.84]); %n=3L
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.23 0.11 0 pi -13 2.22 2.84]); %n=3S
%  X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether5); %n=4
% X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether6); %n =4
%X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5); %n=2
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6; %n=2

X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    indices = find([0 diff(sign(H))]~=0);
    
% plot(z,r); hold on; grid on; plot(z(indices),r(indices),'o'); return
% plot(t,H);  hold on; plot(t,-H);  grid on; return;
  
%% ========== colors ============
%[0, 0.4470, 0.7410]
%[0.8500, 0.3250, 0.0980]
%[0.9290, 0.6940, 0.1250] 	
%[0.4940, 0.1840, 0.5560]
%[0.4660, 0.6740, 0.1880]

%% ------------------ 1st branch --------------------
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
Xvec = [1.04:0.01:1.43 1.4312]/2;
% Xvec = [1.04:-0.01:0]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
        plot(z,r,'b','LineWidth',2);
    end
    if X == Xvec(end)
        plot(z,r,'k','LineWidth',2);
    end
%     if X == Xvec(end)
%         m = sqrt(abs(mu));
%         s =L*(t-0.5);
%         d = sqrt(2/abs(mu));
%         b  = d;
%         %plot(t,r,'LineWidth',2); hold on;
%         vmax = fzero(@(z) 1 - b*cosh(z/b),1);
%         v= linspace(-vmax,vmax);
%         %plot(v-z(end)-v(end),b*cosh(v/b))        
%     end
    
%     figure();
%     plot(t,psi_s + sin(psi)./r)
    y_guess = y;
end

%% ------------------ 2nd branch --------------------
X = 1.24/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.1 0.11 0 pi -3 1.9 1.4]); %n=3T
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

figure(2)
Xvec2 = [1.04:0.01:1.41 1.4155]/2;
% Xvec2 = [1.04:-0.01:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
      plot(z,r,'r','linewidth',2);
    end
    if abs(X - Xvec2(end)) < 0.0001
       plot(z,r,'k','linewidth',2); 
    end
    y_guess = y;
end
 
%% ------------------ 3rd branch --------------------
X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.14 0.1 0 pi 0 1.7 0]); %n=3L
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

figure(2);
Xvec3 = [1.04:0.01:1.41 1.415 1.4155]/2;
% Xvec3 = [1.04:-0.01:0]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec3(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
      plot(z,r,'g','linewidth',2);
    end
    y_guess = y;
end

%% ------------------ 4th branch --------------------
X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -1 1.7 1]); %n=3S
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

figure(1);
Xvec4 = [1.4311 1.43:-0.01:1.28]/2;
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec4(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
      plot(z,r,'m','linewidth',2);
    end
    y_guess = y;
end

figure(1)
xlabel('z'); ylabel('r'); axis equal; box on;
title(['3rd mode, A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure(2)
xlabel('z'); ylabel('r'); axis equal; box on;
title(['3rd mode, A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure();
subplot(1,3,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
plot(Xvec3,Fvec3,'go-','LineWidth',2); hold on;
plot(Xvec4,Fvec4,'mo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18); 

subplot(1,3,2)
plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
plot(Xvec3,muvec3,'go-','LineWidth',2); hold on;
plot(Xvec4,muvec4,'mo-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
plot(Xvec3,Evec3,'go-','LineWidth',2); hold on;
plot(Xvec4,Evec4,'mo-','LineWidth',2); hold on;
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');


end

function [] = tethertransition4_mode4()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A =  2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess (A = 7.59)
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.7 0]); %n=1U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
%  X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.07 0.1 0 pi 0 1.7 0]); %n=3S
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.05 0.1 0 pi 0 1.7 0]); %n=3L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.95 0.1 0 pi 0 1.7 0]); %n=5L
% X = 1.04/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.081 0.11 0 pi -3 1.9 1.4]); %n=7L
% X = 1.24/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.1 0.11 0 pi -3 1.9 1.4]); %n=3T
% X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.7 0]); %n=5S
% X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.14 0.1 0 pi 0 1.7 0]); %n=3L
% X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -1 1.7 1]); %n=3S

 %X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.25 0.11 0 pi -1 1.7 1]);

%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi 0 2 0]); %<--- n=1 upper branch
% solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %n=1 upper branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.11 0.1 0 pi 0 2 0]);%n=7 lower branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]);%n=5 lower branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 2 0]);%n=7 upper branch
%X = 1.54/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 2 0]);%n=7 upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 2 0]); %n=1 lower branch (A=7.59)
%   X = 1.4/2;
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.3 0.1 0 pi -13 2.22 2.84]); %n=3
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.09 0.1 0 pi -13 2.22 2.84]); %n=5L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.12 0.1 0 pi -13 2.22 2.84]); %n=5U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.14 0.1 0 pi -13 2.22 2.84]); %n=3U
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -0.01 0.22 0.1 0 pi -13 2.22 2.84]); %n=3L
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.1 0 pi -13 2.22 2.84]); %n=7L
%X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.12 0.1 0 pi -13 2.22 2.84]); %n=7U
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.13 0.11 0 pi -13 2.22 2.84]); %n=3T
%  X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 4 0.08 0.11 0 pi -13 2.22 2.84]); %n=9L
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.15 0.1 0 pi -13 2.22 2.84]); %n=3U
% X = 1.54/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.22 0.11 0 pi -13 2.22 2.84]); %n=3S
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.15 0.11 0 pi -3 3.52 1.84]); %n=3S
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.21 0.11 0 pi -3 3.52 1.84]); %n=5T
%X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.28 0.11 0 pi -7.4 3.54 2.51]); %n=??
% X = 1.64/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.27 0.11 0 pi -7.4 3.54 2.51]); %n=??
%X = 1.66/2; solinit = bvpinit(linspace(0,1,10),[pi/2 -2 0.21 0.11 0 pi -3 3.52 1.84]); %n=5S
%X = 1.7/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.08 0.11 0 pi -13 2.22 2.84]); %n=9L
%X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -13 2.22 2.84]); %n=5T
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.25 0.11 0 pi -13 2.22 2.84]); %n=3L
% X = 1.8/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.23 0.11 0 pi -13 2.22 2.84]); %n=3S
%  X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether5); %n=4
X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether6); %n =4
%X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5); %n=2
% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6; %n=2

options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    muvec = y(7,1)
    L = y(8,1)
    Fvec = -2*pi*y(9,1)
    indices = find([0 diff(sign(H))]~=0);
    plot(z,r,'b','linewidth',2);
    Xvec = 0.59;
    Evec = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    
    
%  plot(z,r); hold on; grid on; plot(z(indices),r(indices),'o'); return
% plot(t,H);  hold on; grid on; return;
  
%% ========== colors ============
%[0, 0.4470, 0.7410]
%[0.8500, 0.3250, 0.0980]
%[0.9290, 0.6940, 0.1250] 	
%[0.4940, 0.1840, 0.5560]
%[0.4660, 0.6740, 0.1880]

%% ------------------ 1st branch --------------------
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
%Xvec = [0.58 0.56 0.54];
Xvec = [Xvec 0.61 0.65 0.68 0.705 0.708 0.709];
% Xvec = [0.59:0.01:0.7 0.705];
% Xvec = [1.04:-0.01:0]/2;
muvec = [muvec zeros(1,length(Xvec)-1)];
Fvec = [Fvec zeros(1,length(Xvec)-1)];
Evec = [Evec zeros(1,length(Xvec)-1)];
figure(1); hold on;  box on; grid on; axis equal;
for j = 2:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,1) == 0
        plot(z,r,'b','LineWidth',2);
    end
    if X == Xvec(end)
        plot(z,r,'k','LineWidth',2);
    end
 
    y_guess = y;
end


%% ------------------ 2nd branch --------------------
X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether5); %n=4
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
figure(1)
Xvec2 = [0.59:0.01:0.7 0.705 0.7093];
%  Xvec2 = [0.59:-0.01:0];
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
      plot(z,r,'r','linewidth',2);
    end
    if abs(X - Xvec2(end)) < 0.0001
       plot(z,r,'k','linewidth',2); 
    end
    
    y_guess = y;
end
 

figure(1)
xlabel('z'); ylabel('r'); axis equal; box on;
title(['4th mode, A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure();
subplot(1,3,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Fvec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,Fvec4,'mo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18); 

subplot(1,3,2)
plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,muvec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,muvec4,'mo-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
% plot(Xvec3,Evec3,'go-','LineWidth',2); hold on;
% plot(Xvec4,Evec4,'mo-','LineWidth',2); hold on;
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');


end

function [] = tethertransition4_mode5()
 
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A =  2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
global c0; c0 = 0;
 
%% prepare a solution as initial guess
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.7 0]); %n=1U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
% X = 1.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.13 0.1 0 pi -70 1.3 90]); %n=5S
% X = 1.4/2; solinit = bvpinit(linspace(0,1,9),[pi/2 2 0.14 0.1 0 pi -70 1.3 90]); %n=5T
% 
 X = 0.3; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric16); %n=5S
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    y(7,1)
    y(8,1)
    y(9,1)
    indices = find([0 diff(sign(H))]~=0);
    
plot(z,r); hold on; grid on; plot(z(indices),r(indices),'o'); return
% % % plot(t,H);  hold on; plot(t,-H);  grid on; return;
  
%% ========== colors ============
%[0, 0.4470, 0.7410]
%[0.8500, 0.3250, 0.0980]
%[0.9290, 0.6940, 0.1250] 	
%[0.4940, 0.1840, 0.5560]
%[0.4660, 0.6740, 0.1880]

%% ------------------ 1st branch --------------------
disp('branch 1')
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
Xvec = [1.05:0.01:1.40 1.415 1.4161]/2;
%  Xvec = [0.1 0.2 0.3 0.32 0.34]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,1) == 0
        plot(z,r,'b','LineWidth',2);
    end
    if X == Xvec(end)
        plot(z,r,'k','LineWidth',2);
    end

    y_guess = y;
end



%% ------------------ 2nd branch --------------------
disp('branch 2')
 X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

figure(2)
Xvec2 = [1.04:0.01:1.41 1.415 1.4172]/2;
% Xvec2 = [1.04:-0.01:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
      plot(z,r,'r','linewidth',2);
    end
    if abs(X - Xvec2(end)) < 0.0001
       plot(z,r,'k','linewidth',2); 
    end
    y_guess = y;
end
 
%% ------------------ 3rd branch --------------------
disp('branch 3')
X = 1.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.13 0.1 0 pi -70 1.3 90]); %n=5S
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

figure(2);
Xvec3 = [1.4174 1.415 1.4:-0.01:1.3 1.295:-0.005:1.225 1.2]/2;
% Xvec3 = [1.04:-0.01:0]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec3(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,3) == 1
      plot(z,r,'g','linewidth',2);
    end
    y_guess = y;
    
end


%% ------------------ 4th branch --------------------
disp('branch 4')
X = 1.4/2; solinit = bvpinit(linspace(0,1,9),[pi/2 2 0.14 0.1 0 pi -70 1.3 90]); %n=5T
options = bvpset('RelTol',1e-9);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

figure(1);
Xvec4 = [1.415 1.2989]/2;
% Xvec4 = [1.415 1.395:-0.01:1.315 1.3149 1.3148]/2;
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j)
    solinit = bvpinit(linspace(0,1,9),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec4(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,1) == 0 || j == length(Xvec4)
       plot(z,r,'m','linewidth',2);
    end
    y_guess = y;
    
    
end

figure(1)
xlabel('z'); ylabel('r'); axis equal; box on;
title(['5th mode, A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure(2)
xlabel('z'); ylabel('r'); axis equal; box on;
title(['5th mode, A = 2\pi*' num2str(A/2/pi)]);
set(gcf,'Color','w');
set(gca,'FontSize',18);

figure();
subplot(1,3,1)
plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
plot(Xvec3,Fvec3,'go-','LineWidth',2); hold on;
plot(Xvec4,Fvec4,'mo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18); 

subplot(1,3,2)
plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
plot(Xvec3,muvec3,'go-','LineWidth',2); hold on;
plot(Xvec4,muvec4,'mo-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
plot(Xvec3,Evec3,'go-','LineWidth',2); hold on;
plot(Xvec4,Evec4,'mo-','LineWidth',2); hold on;
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18); 
set(gcf,'Color','w');


end

%% explores nonzero kbar, mu > 0 (with BL)
function [] = kbarsolutions()
 
global kappa; kappa = 1;
global kbar; kbar = 0.5;
global A; A = 2*pi*1; %surface area
global X; X = 0.5; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
 
%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %wacky
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 2 0]); %lower branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1 0]); %lower branch
options = bvpset('RelTol',1e-8);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 y  = y_guess; r = y(3,:); z = y(4,:); %plot(z,r)
 y(7,1)
%return;
 

%% calculate shape for various widths
bcat = 0.8255174536525; Lcat = 2*sqrt(1-bcat^2); Xcat = 0.527697396963; %v = linspace(-0.527697,0.527697);
%lower branch
%Xvec = [0.45:0.01:0.52 0.521:0.001:0.523 0.5231:0.0001:0.5238 0.52381:0.00001:0.52386 0.5239 0.524 0.526 0.527 0.5265 0.5264 0.5263 0.526];
%upper branch
Xvec = [0.45:0.01:0.52 0.521:0.001:0.523 0.5231:0.0001:0.5238 0.52381:0.00001:0.52386 0.5239 0.524 0.526 0.527 0.5273 0.5274 0.5275 0.5276 0.52765 0.52769 0.527693 0.527697 0.52769739];%Xvec = [1.053:-0.01:1 0.95:-0.05:0.2]/2; %(2:-0.1:0.1)/2; 
%Xvec = [1.053/2 0.527 0.5271:0.0001:0.5276 0.52765 0.52766 0.52767 0.52768 0.52769 0.527695 0.527696 0.5276965 0.5276969 0.5276971]
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
Avec = zeros(1,length(Xvec));
Lvec = zeros(1,length(Xvec));
Avec2 = zeros(1,length(Xvec));
qvec = zeros(1,length(Xvec));
figure; hold on; %axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N); s = Lcat*(t-0.5);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); gam = y(5,:); ar = y(6,:);
    mu = y(7,1);  eta = y(9,1); g1 = gam(1); psi1 = psi(1); L = y(8,1)
    muvec(j) = y(7,1); L = y(8,1);
    H = (psi_s+sin(psi)./r)/2;
    K = psi_s.*sin(psi)./r;
    Ham = r/2.*(psi_s.^2-(sin(psi)./r).^2) - mu*r + gam.*cos(psi) - eta*sin(psi);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Lvec(j) = L;
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) =  2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r)+2*pi*kbar*(0-2*cos(psi(1)));
    hold on;  grid on;

    eps = abs(L-Lcat)/Lcat
    eps2 = kbar*bcat^2/mu
    dz = abs(Xcat-X)/Xcat
    
    rcat = sqrt((Lcat*(t-0.5)).^2 + bcat^2);
    zcat = -bcat*atanh(s./sqrt(s.^2+bcat^2));

    Avec(j) = trapz(t,r-rcat)*L;
    
    C = (kbar*bcat^2/mu + eps*Lcat^2/4)/eps; %matching constant
    L0=2*sqrt(1-bcat*bcat);
    router = -4*(1-bcat^2)*(t-1/2).^2./rcat + C./rcat;
    rinner = (-kbar*bcat^2/mu*exp(-sqrt(mu)*(Lcat/2-s))+kbar*bcat^2/mu./rcat)/eps;
    rinner = [fliplr(rinner(502:end)) rinner(501:end)]; %use even symmetry
    rcomp = router + rinner - kbar*bcat^2/mu/eps;
    Avec2(j) = trapz(t,r - rcat)*Lcat;

    qvec(j) = kbar*bcat/mu^(3/2)*2;
    
   if X == Xvec(end) 
       plot(t-0.5,H,'b','linewidth',2); hold on;
       plot(t-0.5,kbar*bcat/2*cosh(s*sqrt(mu))./cosh(Lcat/2*sqrt(mu)),'r','linewidth',2); hold on; grid on;
       set(gcf,'color','w');
       set(gca,'fontsize',18);
       box on; grid on; xlabel('$t$','interpreter','latex');ylabel('$Ha$','interpreter','latex');
       legend('numerical sol.','cosh approx.')
       %        plot(t-0.5,kbar*bcat/2*exp((s-L/2)*sqrt(mu)),'o')

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% positive tension case (upper branch)
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure();
    mu
%     plot(t-0.5,-(r-rcat)/eps,'b','linewidth',2); hold on; grid on;%perturbation to catenoid
    C = (kbar*bcat^2/mu + eps*Lcat^2/4)/eps; %matching constant
    
    delta = 1/sqrt(abs(mu));
    router = -4*(1-bcat^2)*(t-1/2).^2./rcat + C./rcat;
    %      router = (-eps*s.^2 + C)./rcat;
    %     rinner = (-bcat^2*kbar*exp(exp(-(Lcat/2-s)*sqrt(mu)) ))/eps/mu + kbar*bcat^2/mu/eps;
    rinner = (-kbar*bcat^2/mu*exp(-sqrt(mu)*(Lcat/2-s))+kbar*bcat^2/mu)/eps;
    rinner = [fliplr(rinner(502:end)) rinner(501:end)]; %use even symmetry
%     rcomp = router + rinner - kbar*bcat^2/mu/eps;
    rcomp = (eps*(s.^2-Lcat^2/4)./rcat + kbar*bcat^2/mu*(exp(-(s+Lcat/2)/delta)+exp((s-Lcat/2)/delta) - 1./rcat))/eps;
    s2 = Lcat*linspace(0,1,N);
    %   rinner_trp = -kbar*bcat^2/mu*(exp(-s2*sqrt(mu)) + exp(-(Lcat-s2)*sqrt(mu)) - 1) ;
    % 	router_trp = -eps*(s2-Lcat/2).^2./sqrt((s2-Lcat/2).^2 + bcat^2) + eps*C./sqrt((s2-Lcat/2).^2 + bcat^2);
    %   router_trp = -eps.*s2.*(s2-Lcat)./sqrt((s2-Lcat/2).^2 + bcat^2) + kbar*bcat^2/mu./sqrt(bcat^2 + (s2-Lcat/2).^2);
    %     plot(t-0.5,(kbar*bcat^2/mu*(exp(-s2*sqrt(mu)) + exp(-(Lcat-s2)*sqrt(mu)) - 1./sqrt((s2-Lcat/2).^2 + bcat^2)) + eps*s2.*(s2-Lcat)./sqrt((s2-Lcat/2).^2 + bcat^2))/eps,'p')
%     plot(t-0.5,-router,'r--','linewidth',2);
%     plot(t-0.5,-rinner,'g--','linewidth',2);
%     plot(t-0.5,rcomp/eps,'ko','linewidth',2);

%     plot(t(1:end-1)-0.5,diff(z-zcat)./diff(t)/L); hold on;
%     plot(t(1:end-1)-0.5,Lcat/bcat*(t(1:end-1)-0.5).*diff(r-rcat)./diff(t)/L+eps*rcat(1:end-1)/bcat)

    plot(t-0.5,(z - zcat)/eps,'linewidth',2); hold on;
    plot(t-0.5,s/bcat.*(r-rcat)/eps - 1/(1+eps)*L*cumtrapz(t-0.5,r-rcat)/bcat/eps + (Xcat - X) + L*cumtrapz(t-0.5,rcat)/bcat,'linewidth',2);
%     plot(t-0.5,1/4/bcat*(-((2*bcat^2+Lcat^2)*eps*atanh(s./rcat)) + (s.*(4*bcat^2*kbar+(2*bcat^2+Lcat^2-2*s.^2)*eps*mu)./rcat - 4*bcat^2*kbar*(log(s + rcat)+2*exp(-Lcat/2/delta).*(s.*cosh(s/delta) - delta*sinh(s/delta))))/mu),'linewidth',2)
    %      plot(t-0.5,(kbar*bcat^2/mu*(exp(-s2*sqrt(mu)) + exp(-(Lcat-s2)*sqrt(mu)) - 1./sqrt((s2-Lcat/2).^2 + bcat^2)) + eps*s2.*(s2-Lcat)./sqrt((s2-Lcat/2).^2 + bcat^2))/eps,'mp')
    %      plot(t-0.5,rinner_trp/eps,'linewidth',2);
    %      plot(t-0.5,router_trp/eps,'linewidth',2);
    box on; grid on;
    legend('num. solution','outer approx.','inner approx.','TRP solution')
    xlabel('$t$','interpreter','latex'); ylabel('$z_1/a$','interpreter','latex');
    set(gca,'FontSize',18);
    set(gcf,'color','w');
    
   end
   

    y_guess = y;
end


%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec2 = [];
%Xvec2 = [ 1.054:-0.01:1 0.95:-0.05:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    muvec2(j) = y(7,1); L = y(8,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    H = (psi_s+sin(psi)./r)/2;
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r)+2*pi*kbar*(0-2*cos(psi(1)));
    %check = 2*pi*L*trapz(t,psi_s.*sin(psi)) - 2*pi*(0+2*cos(psi(1)))
    hold on;  grid on;
    plot(z,H,'r','LineWidth',2);
%     z2= z./(r.*r+z.*z)*(1+z(end)^2);
%     r2 = r./(r.*r+z.*z)*(1+z(end)^2);
%     plot(z2,r2,'r--');
    %Psi = acos(((r.*r+z.*z).*cos(psi) - 2*r.*(r.*cos(psi) - z.*sin(psi)))./(r.*r+z.*z).^2);
    %plot(t,Psi);
   % y_guess = y;
end

% %% prepare a solution as initial guess
% solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
% sol = bvp4c(@odesystem,@bcs,solinit);
% t = linspace(0,1,N);
% y_guess = deval(sol,t);
%  
% Xvec3 = [1.05 0.527*2]/2;
% muvec3 = zeros(1,length(Xvec3));
% Fvec3 = zeros(1,length(Xvec3));
% for j = 1:length(Xvec3)
%     X = Xvec3(j);
%     solinit = bvpinit(linspace(0,1,11),@newguess);
%     sol = bvp4c(@odesystem,@bcs,solinit);
%     t = linspace(0,1,N);
%     y = deval(sol,t);
%     %r = y(3,:); z = y(4,:); psi_s = y(2,:);
%     muvec3(j) = y(7,1);
%     %b = r((N+1)/2); k1 = psi_s((N+1)/2);
%     %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
%     Fvec3(j) = -2*pi*y(9,1);
%     hold on;  grid on;
%    % plot(z,r,'LineWidth',2);
%    % y_guess = y;
% end


 
%% catenoid solution
% bcat = 0.825517; v = linspace(-0.527697,0.527697);
%plot(v,bcat*cosh(v/bcat),'k','LineWidth',3)
%% plot(v./(v.*v+bcat^2*cosh(v/bcat).^2)*(v(1)^2+bcat^2*cosh(v(1)/bcat)^2),...
%%     bcat*cosh(v/bcat)./(v.*v+bcat^2*cosh(v/bcat).^2)*(v(1)^2+bcat^2*cosh(v(1)/bcat)^2),'k--');
% xlabel('$t$','interpreter','latex'); ylabel('$r_1/a$','interpreter','latex');
% title(['Matched asymptotic expansion of r_1 (kbar/\kappa = ' num2str(kbar/kappa) ')']);


figure();
% subplot(1,3,1)
% plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
% plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
% %plot(Xvec3(end),Fvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
% plot([0.527697 0.527697],[-250 250],'k--');
% title('Force vs. extension')
% xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
% set(gca,'FontSize',18);
%  
% subplot(1,3,2)
% plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
% plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
% %plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
% plot([0.527697 0.527697],[-50 200],'k--');
% title('Tension vs. extension')
% xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% 
% subplot(1,3,3)
%plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
%plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
%plot(0.527697,-4*pi*kbar*cos(atan(sinh(0.527697/0.825517))+pi/2),'kp','MarkerSize',18,'MarkerFaceColor','k')
plot(log(abs(Lcat-Lvec)/Lcat)/log(10),log(2*abs(Xcat-Xvec))/log(10),'bo-','linewidth',2); hold on;
plot(log(abs(Lcat-Lvec)/Lcat)/log(10),log(abs(qvec))/log(10),'ro-','linewidth',2); hold on;
% plot([-4 -8],[2 (2+2*4/3)],'linewidth',2)
% title('Tension vs. |h-h_m_a_x| (loglog)');
xlabel('$\log(|L-L_0|/L_0)$','interpreter','latex'); ylabel('$\log(\mu a^2/\kappa)$','interpreter','latex'); grid on; box on;
plot([-5 -2],[4 1],'--','linewidth',2)
%%plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2)% factor of 2
%%plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2) %from extension/2
%%plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
%%plot([0.527697 0.527697],[-50 200],'k--');
%%title('Energy vs. extension')
%%xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

end

%% explores nonzero kbar, mu < 0 (with WKB)
function [] = kbarsolutions2()
 
global kappa; kappa = 1;
global kbar; kbar = 0.5;
global A; A = 2*pi*1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints
 
%% prepare a solution as initial guess
%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %upper branch
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %wacky
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 2 0]); %lower %branch
% solinit = bvpinit(linspace(0,1,5),@guess4_kbar); %m = 4
% X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
X = 1.04/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.09 0.11 0 pi -1300 2.22 2909]); %n=9L
X = 1.04/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.085 0.11 0 pi -1300 2.22 2909]); %n=9L
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
y = y_guess;
plot(y(4,:),y(3,:));
mu = y(7,1)
delta = 1/sqrt(abs(mu))
% return


%% calculate shape for various widths
bcat = 0.8255174536525; Lcat = 2*sqrt(1-bcat^2); Xcat = 0.527697396963; %v = linspace(-0.527697,0.527697);
%lower branch
Xvec = [0.45:0.01:0.52 0.521:0.001:0.523 0.5231:0.0001:0.5238 0.52381:0.00001:0.52386 0.5239 0.524 0.526 0.5265 0.5267 0.5268 ];%0.5265 0.5264 0.5263 0.526 0.52 0.5 0.46];
%3rd mode
%Xvec = [0.51 0.52 0.526 0.527 0.5273 0.5274 0.52742 0.52745 0.52746:0.000001:0.527463];
%use the vector below with guess4_kbar
Xvec = [0.523 0.526 0.5267 0.527 0.5272 0.5274 0.527425 0.52744 0.52745 0.5275 0.5276 0.52763 0.527635 0.527638 0.5276385 0.5276386 0.5276387];
Xvec = [0.523 0.526 0.5267 0.527 0.5272 0.5274]% 0.527425 0.52744 0.52745 0.5275 0.5276 0.52764  0.52765 0.52767 ];
% Xvec = [0.51 0.52 0.522 0.5274]; %Xcat = 0.527697396963

muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
Avec = zeros(1,length(Xvec));
Lvec = zeros(1,length(Xvec));
figure; hold on; %axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N); s = Lcat*(t-0.5);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); gam = y(5,:); ar = y(6,:);
    mu = y(7,1);  eta = y(9,1); g1 = gam(1); psi1 = psi(1); L = y(8,1);
    muvec(j) = y(7,1); L = y(8,1);
    H = (psi_s+sin(psi)./r)/2;
    K = psi_s.*sin(psi)./r;
    Ham = r/2.*(psi_s.^2-(sin(psi)./r).^2) - mu*r + gam.*cos(psi) - eta*sin(psi);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Lvec(j) = y(8,1);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r)+2*pi*kbar*(0-2*cos(psi(1)));
    hold on;  grid on;

    eps = abs(L-Lcat)/Lcat;    
    delta = 1/abs(mu);
    n_int = L/pi/sqrt(delta)
    zero_test = cos(L/2/sqrt(delta))

    %plot(z,r); hold on;
    %plot(t-0.5,r1);
    Hh = kbar*bcat/2/(1-Lcat/2*asinh(L/2/bcat))*(1-s./(s.^2+bcat^2).*asinh(s./bcat));
    Hwkb = kbar*bcat/2*cos(s*sqrt(-mu))/cos((sqrt(-mu))*Lcat/2);
    
    Kcat = -bcat^2./(s.^2+bcat^2).^2;
    rcat = sqrt(s.^2 + bcat^2);
    %plot(t,H,'b','linewidth',2)
    %plot(t,kbar*bcat/2*cosh(s*sqrt(mu))/cosh(L*sqrt(mu)/2),'r','linewidth',2)
    
    %r1 = -kbar*bcat/mu*(cosh(s*sqrt(mu))/cosh(Lcat*sqrt(mu)/2));
    Avec(j) = trapz(t,r-rcat)*L;
    

    if X == Xvec(end)
        mu
        delta = 1/abs(mu)
        plot(t-0.5,H,'b','linewidth',2); hold on; grid on;
        %plot(t-0.5,kbar*bcat/2*cosh(s*sqrt(mu))./cosh(Lcat/2*sqrt(mu)),'linewidth',2);
        plot(t-0.5,kbar*bcat/2*cosh(s*sqrt(mu))./cosh(Lcat/2*sqrt(mu))./(s.^2+bcat^2).^(1/4),'r','linewidth',2);
        %plot(t-0.5,s.*asinh(s/bcat)./sqrt(s.^2+bcat^2)*CC + 2*cosh(s*sqrt(mu)));
        %plot(t-0.5, kbar*bcat/2*exp((s-Lcat/2)*sqrt(abs(mu))),'linewidth',2);
        u = s/sqrt(delta);
        %plot(t-0.5,1.8*cos(u)./rcat.^(1/2),'linewidth',2)
        kbar*bcat/2/0.0467
        CC = kbar*bcat/2/cos(Lcat/2/sqrt(delta))
        DD = kbar*bcat/2/cos(Lcat/2/sqrt(delta) + sqrt(delta)/8*(5*Lcat/4+7/2/bcat*atan(Lcat/2/bcat)))
%         plot(t-0.5,DD./sqrt(rcat).*cos(s/sqrt(delta) + sqrt(delta)/8*(5*s/2./(s.^2+bcat^2) + 7*atan(s/bcat)/2/bcat)),'linewidth',2);
        %plot(t-0.5, 2*cos(u)+2*delta*kbar*bcat/2*(u.*sin(u) + 2*(1+bcat^2)*cos(u)),'linewidth',2)
        legend('numerical sol.','WKB approx.'); box on; grid on;
%         title('Tension vs. extension')
        xlabel('$t$','interpreter','latex'); ylabel('$Ha$','interpreter','latex');  grid on; box on;
        set(gca,'FontSize',18);
        set(gcf,'color','w');

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% positive tension case (upper branch)
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      plot(t-0.5,r-rcat,'b','linewidth',2); %perturbation to catenoid
%      C = kbar*bcat^2/mu + eps*Lcat^2/4; %matching constant
%      router = (-eps*s.^2 + C)./rcat;
%      rinner = -kbar*bcat^2/mu*exp(-sqrt(mu)*(Lcat/2-s))+kbar*bcat^2/mu;
%      rinner = [fliplr(rinner(52:end)) rinner(51:end)]; %use even symmetry
%      rcomp = router + rinner - kbar*bcat^2/mu;
%      plot(t-0.5,router,'r--','linewidth',2);
%      plot(t-0.5,rinner,'g--','linewidth',2);
%      plot(t-0.5,rcomp,'k','linewidth',2);
%      box on; grid on;
%      legend('num. solution','outer approx.','inner approx.','composite')
     
     
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% negative tension case (lower branch)
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot(t-0.5,-kbar*bcat/mu*(cosh(s*sqrt(mu))/cosh(Lcat*sqrt(mu)/2) - 1/(1-Lcat/2*asinh(Lcat/2/bcat))*(1-s./(s.^2+bcat^2).*asinh(s./bcat))),'r','linewidth',2);
    %plot(t-0.5,-kbar*bcat/mu*(cosh(s*sqrt(mu))/cosh(Lcat*sqrt(mu)/2) - 1./rcat),'r','linewidth',2);

    figure();
    plot(t-0.5,(r-rcat)/eps,'b','linewidth',2); hold on; grid on;
    D = sec(Lcat/2*sqrt(-mu));
    C = kbar*bcat;
    
    r1approx = -1/mu*C*bcat./rcat./sqrt(rcat).*cos(s*sqrt(-mu))*D; %r_osc, note sqrt(rcat) [higher order approx]
    M = -r1approx(1);
    r1approx_hom = M+2*eps-2*eps*rcat-2*sqrt(1-bcat^2)*eps*atanh(sqrt(1-bcat^2))+2*eps*s.*atanh(s./rcat);
    plot(t-0.5,r1approx/eps,'r--','linewidth',2)
    plot(t-0.5,r1approx_hom/eps,'g--','linewidth',2)
    plot(t-0.5,(r1approx + r1approx_hom)/eps,'ko','linewidth',2);
    legend('numerical sol.','oscillating sol.','remaining sol.','total approx.')

%       plot(t-0.5,H,'b','linewidth',2)
%       plot(t-0.5,Hwkb,'r','linewidth',2);
%       sqrt_mu_over_twopi =  sqrt(abs(mu))/2/pi

      %plot(z,r,'k')
    
    
   end
   

    y_guess = y;
end
%plot(t-0.5,-bcat^2./((Lcat*(t-0.5)).^2+bcat^2).^2,'ko-')


%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec2 = [];
%Xvec2 = [ 1.054:-0.01:1 0.95:-0.05:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    muvec2(j) = y(7,1); L = y(8,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    H = (psi_s+sin(psi)./r)/2;
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r)+2*pi*kbar*(0-2*cos(psi(1)));
    %check = 2*pi*L*trapz(t,psi_s.*sin(psi)) - 2*pi*(0+2*cos(psi(1)))
    hold on;  grid on;
    plot(z,H,'r','LineWidth',2);
%     z2= z./(r.*r+z.*z)*(1+z(end)^2);
%     r2 = r./(r.*r+z.*z)*(1+z(end)^2);
%     plot(z2,r2,'r--');
    %Psi = acos(((r.*r+z.*z).*cos(psi) - 2*r.*(r.*cos(psi) - z.*sin(psi)))./(r.*r+z.*z).^2);
    %plot(t,Psi);
   % y_guess = y;
end

% %% prepare a solution as initial guess
% solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
% sol = bvp4c(@odesystem,@bcs,solinit);
% t = linspace(0,1,N);
% y_guess = deval(sol,t);
%  
% Xvec3 = [1.05 0.527*2]/2;
% muvec3 = zeros(1,length(Xvec3));
% Fvec3 = zeros(1,length(Xvec3));
% for j = 1:length(Xvec3)
%     X = Xvec3(j);
%     solinit = bvpinit(linspace(0,1,11),@newguess);
%     sol = bvp4c(@odesystem,@bcs,solinit);
%     t = linspace(0,1,N);
%     y = deval(sol,t);
%     %r = y(3,:); z = y(4,:); psi_s = y(2,:);
%     muvec3(j) = y(7,1);
%     %b = r((N+1)/2); k1 = psi_s((N+1)/2);
%     %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
%     Fvec3(j) = -2*pi*y(9,1);
%     hold on;  grid on;
%    % plot(z,r,'LineWidth',2);
%    % y_guess = y;
% end

set(gca,'FontSize',18);
set(gcf,'color','w');
 
%% catenoid solution
bcat = 0.825517; v = linspace(-0.527697,0.527697);
%plot(v,bcat*cosh(v/bcat),'k','LineWidth',3)
%% plot(v./(v.*v+bcat^2*cosh(v/bcat).^2)*(v(1)^2+bcat^2*cosh(v(1)/bcat)^2),...
%%     bcat*cosh(v/bcat)./(v.*v+bcat^2*cosh(v/bcat).^2)*(v(1)^2+bcat^2*cosh(v(1)/bcat)^2),'k--');
xlabel('$t$','interpreter','latex'); ylabel('$r_1/a$','interpreter','latex');
% title(['WKB approximation of r_1 (kbar/\kappa = ' num2str(kbar/kappa) ')']);

figure();
% subplot(1,3,1)
% plot(Xvec,Fvec,'bo-','LineWidth',2); hold on;
% plot(Xvec2,Fvec2,'ro-','LineWidth',2); hold on;
% %plot(Xvec3(end),Fvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
% plot([0.527697 0.527697],[-250 250],'k--');
% title('Force vs. extension')
% xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
% set(gca,'FontSize',18);
%  
% subplot(1,3,2)
% plot(Xvec,muvec,'bo-','LineWidth',2); hold on;
% plot(Xvec2,muvec2,'ro-','LineWidth',2); hold on;
% %plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
% plot([0.527697 0.527697],[-50 200],'k--');
% title('Tension vs. extension')
% xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% 
% subplot(1,3,3)
% Lvec = [Lvec(1:end-3) fliplr(Lvec(end-2:end))];
% muvec = [muvec(1:end-3) fliplr(muvec(end-2:end))];
%plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
%plot(Xvec2,Evec2,'ro-','LineWidth',2); hold on;
%plot(0.527697,-4*pi*kbar*cos(atan(sinh(0.527697/0.825517))+pi/2),'kp','MarkerSize',18,'MarkerFaceColor','k')
% plot(log(abs(Xcat-Xvec)/Xcat)/log(10),log(abs(muvec))/log(10),'bo-','linewidth',2); hold on;
plot(log(abs(Lcat-Lvec)/Lcat)/log(10),log(abs(muvec))/log(10),'bo-','linewidth',2); hold on;
% plot(Xvec,Lvec,'o-'); hold on;
% plot(Xcat,Lcat,'p')
%plot([-1.5 2.576],[-2 ],'linewidth',2)
% title('Tension vs. |h-h_m_a_x| (loglog)');
xlabel('$\log(|L-L_0|/L_0)$','interpreter','latex'); ylabel('$\log(\mu a^2/\kappa)$','interpreter','latex'); grid on; box on;
plot([-3.7 -3.3],[3.3 2.9],'--','linewidth',2)
%%plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2)% factor of 2
%%plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2) %from extension/2
%%plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
%%plot([0.527697 0.527697],[-50 200],'k--');
%%title('Energy vs. extension')
%%xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

end

function [] = kbarcomparison()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
 
%% prepare a solution as initial guess
%solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %wacky
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi 0 2 0]); %also wacky, A = 1.11*2*pi
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]);

% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.200001 0.1 0 pi 0 2 0]); %m = 5b
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %m = 5
% solinit = bvpinit(linspace(0,1,5),@guess4); %m = 4
%solinit = bvpinit(linspace(0,1,11),@guess3); %m = 3
%solinit = bvpinit(linspace(0,1,9),@guess3b); %m = 3b
%solinit = bvpinit(linspace(0,1,11),@guess2); %m = 2
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
%X = 0.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.6 0]); %5th order, lower branch
% X = 0.25; solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 0.08 0.1 0 pi 0 1.2 -80]); 
% X = 0.25; solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 0.09 0.1 0 pi 0 1.2 -80]); 
% X = 0.25; solinit = bvpinit(linspace(0,1,11),[pi/2 0.1 0.30002 0.1 0 pi 5 1.2 15]); 
% X = 0.25; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.13 0.1 0 pi 50 1 15]); 
% X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 -20 0.11 0.1 0 pi 70 1.2 70]); 
% X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 -20 0.21 0.1 0 pi 70 1.2 70]); 
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
gam = y(5,:);
mu = y(7,1)
L = y(8,1)
eta = y(9,1)
% plot(z,r); return;

% psi_ss = 1./r.*(-cos(psi).*psi_s + sin(psi).*cos(psi)./r + gam.*sin(psi) + eta*cos(psi));
% C = psi_ss.^2 + psi_s.^4/4 - mu*psi_s.^2;
% plot(t,C)
% %plot(z,r);
% return;

%%%%%%%%%%%%%%% 
%lower branches
%%%%%%%%%%%%%%%

%% calculate shape for various widths
Xvec = [1.0553 1.05 1.04:-0.01:0.14 0.135 0.133]/2; %(2:-0.1:0.1)/2; 
% Xvec = [0.4:-0.01:0.1];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
evvec = zeros(1,length(Xvec));
hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 7
    plot(z,r,'k','LineWidth',2);
    end
    y_guess = y;
    
        %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'k-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'k-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'k-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'k-.','linewidth',2); hold on;
        end
    end
end

kbar = 1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec2 = [1.04:-0.01:0.14 0.135 0.133]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
evvec2 = zeros(1,length(Xvec2));
% hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
    plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec2(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'b-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'b-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-.','linewidth',2); hold on;
        end
    end
end

kbar = 1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec2b = [1.04:0.001:1.047 1.0473 1.0476 1.048]/2;
muvec2b = zeros(1,length(Xvec2b));
Fvec2b = zeros(1,length(Xvec2b));
Evec2b = zeros(1,length(Xvec2b));
evvec2b = zeros(1,length(Xvec2b));
% hold on; axis equal; box on;
for j = 1:length(Xvec2b)
    X = Xvec2b(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec2b(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2b(j) = -2*pi*y(9,1);
    Evec2b(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
%     plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
    
        %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec2b(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'b-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'b-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-.','linewidth',2); hold on;
        end
    end
end

kbar = -1; X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi -10 2 0]); %m = 1
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec3 = [1.04:-0.01:0.14 0.135 0.133]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
evvec3 = zeros(1,length(Xvec3));
% hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
    plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec3(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end

kbar = -1; X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi -10 2 0]); %m = 1
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec3b = [1.04:0.001:1.05 1.051 1.052 1.053 1.054]/2;
muvec3b = zeros(1,length(Xvec3b));
Fvec3b = zeros(1,length(Xvec3b));
Evec3b = zeros(1,length(Xvec3b));
evvec3b = zeros(1,length(Xvec3b));
% hold on; axis equal; box on;
for j = 1:length(Xvec3b)
    X = Xvec3b(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3b(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3b(j) = -2*pi*y(9,1);
    Evec3b(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
%     plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec3b(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end

%%%%%%%%%%%%%%% 
%upper branches
%%%%%%%%%%%%%%%
kbar = 0; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec4 = [1.0553 1.05 1.04:-0.01:0.0]/2;
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
evvec4 = zeros(1,length(Xvec4));
% hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 7
    P1 = plot(z,r,'k','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec4(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'k-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'k-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'k-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'k-.','linewidth',2); hold on;
        end
    end
end

kbar = 1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec5 = [1.04:-0.01:0.0]/2;
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
Evec5 = zeros(1,length(Xvec5));
evvec5 = zeros(1,length(Xvec5));
% hold on; axis equal; box on;
for j = 1:length(Xvec5)
    X = Xvec5(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec5(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec5(j) = -2*pi*y(9,1);
    Evec5(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
    P2 = plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec5(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'b-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'b-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-.','linewidth',2); hold on;
        end
    end
end

kbar = 1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
% Xvec5b = [1.04:0.001:1.055 1.0551 1.0552 1.0553]/2;
Xvec5b = [1.04:0.001:1.053 ]/2;
muvec5b = zeros(1,length(Xvec5b));
Fvec5b = zeros(1,length(Xvec5b));
Evec5b = zeros(1,length(Xvec5b));
evvec5b = zeros(1,length(Xvec5b));
% hold on; axis equal; box on;
for j = 1:length(Xvec5b)
    X = Xvec5b(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec5b(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec5b(j) = -2*pi*y(9,1);
    Evec5b(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
%     plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec5b(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'b-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'b-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'b-.','linewidth',2); hold on;
        end
    end
end

kbar = -1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec6 = [1.04:-0.01:0.0]/2;
muvec6 = zeros(1,length(Xvec6));
Fvec6 = zeros(1,length(Xvec6));
Evec6 = zeros(1,length(Xvec6));
evvec6 = zeros(1,length(Xvec6));
% hold on; axis equal; box on;
for j = 1:length(Xvec6)
    X = Xvec6(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec6(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec6(j) = -2*pi*y(9,1);
    Evec6(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
    P3 = plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec6(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end

kbar = -1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec6b = [1.04:0.001:1.048 1.0485 1.049 1.05]/2;
muvec6b = zeros(1,length(Xvec6b));
Fvec6b = zeros(1,length(Xvec6b));
Evec6b = zeros(1,length(Xvec6b));
evvec6b = zeros(1,length(Xvec6b));
% hold on; axis equal; box on;
for j = 1:length(Xvec6b)
    X = Xvec6b(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec6b(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec6b(j) = -2*pi*y(9,1);
    Evec6b(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
%     plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
    
    %% stability analysis
    [eval,efun]= evmatrix(sol);
    evvec6b(j) = eval;
    if mod(j,4) == 1
        if eval > 0
        subplot(1,2,1)
        plot(z,r,'r-','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-','linewidth',2); hold on;
        else
        subplot(1,2,1)
        plot(z,r,'r-.','linewidth',2); hold on;
        subplot(1,2,2)
        plot(t,efun,'r-.','linewidth',2); hold on;
        end
    end
end

plot([0,0],[0 1],'r','linewidth',2)
bcat = 0.825517; v = linspace(-0.527697,0.527697);
plot(v,bcat*cosh(v/bcat),'k:','LineWidth',2)
legend([P2 P1 P3], {'$\bar\kappa/\kappa = 1$','$\bar\kappa/\kappa = 0$','$\bar\kappa/\kappa = -1$'},'interpreter','latex')
xlabel('$z/a$','interpreter','latex'); ylabel('$r/a$','interpreter','latex');
set(gca,'fontsize',18);
set(gcf,'color','w');

figure();
plot([2*0.527697 2*0.527697],[-250 100],'--','linewidth',2,'color',[0.6 0.6 0.6]); hold on;
g1 = text(1.1,-200,'$h = h^*$','interpreter','latex','color',[0.6 0.6 0.6],'fontsize',18); set(g1,'rotation',90);
p1 = plot(2*Xvec,Fvec,'k-','linewidth',3); hold on; grid on;
p2 = plot(2*Xvec2,Fvec2,'b-','linewidth',3);
plot(2*Xvec2b,Fvec2b,'b-','linewidth',3);
p3 = plot(2*Xvec3,Fvec3,'r-','linewidth',3);
plot(2*Xvec3b,Fvec3b,'r-','linewidth',3);
p4 = plot(2*Xvec4,Fvec4,'k-','linewidth',3);
p5 = plot(2*Xvec5,Fvec5,'b-','linewidth',3);
plot(2*Xvec5b,Fvec5b,'b-','linewidth',3);
p6 = plot(2*Xvec6,Fvec6,'r-','linewidth',3);
plot(2*Xvec6b,Fvec6b,'r-','linewidth',3);
plot([1.055 1.055 ],[-28.94 -30.72],'k','linewidth',3);
plot([0 0.2],[0 0],'r-','linewidth',3)
xlabel('$h/a$','interpreter','latex'); 
ylabel('$Fa/\kappa$','interpreter','latex'); 
legend([p2 p1 p3], {'$\bar\kappa/\kappa = -1$','$\bar\kappa/\kappa = 0$','$\bar\kappa/\kappa = 1$'},'interpreter','latex')
set(gca,'fontsize',18);
set(gcf,'color','w');
axis([0 1.13 -225 75])

figure();
plot([2*0.527697 2*0.527697],[-250 1000],'--','linewidth',2,'color',[0.6 0.6 0.6]); hold on;
g1 = text(1.1,-110,'$h = h^*$','interpreter','latex','color',[0.6 0.6 0.6],'fontsize',18); set(g1,'rotation',90);
q1 = plot([2*Xvec4(1) 2*Xvec],[evvec4(1) evvec],'k-','linewidth',3); hold on; grid on;
q2 = plot(2*Xvec2,evvec2,'b-','linewidth',3);
plot(2*Xvec2b,evvec2b,'b-','linewidth',3);
q3 = plot(2*Xvec3,evvec3,'r-','linewidth',3);
plot(2*Xvec3b,evvec3b,'r-','linewidth',3);
plot(2*Xvec4,evvec4,'k-','linewidth',3);
plot(2*Xvec5,evvec5,'b-','linewidth',3);
plot(2*Xvec5b,evvec5b,'b-','linewidth',3);
plot(2*Xvec6,evvec6,'r-','linewidth',3);
plot(2*Xvec6b,evvec6b,'r-','linewidth',3);
set(gca,'fontsize',18);
set(gcf,'color','w');
axis([0 1.13 -200 1000]);
legend([q2 q1 q3], {'$\bar\kappa/\kappa = -1$','$\bar\kappa/\kappa = 0$','$\bar\kappa/\kappa = 1$'},'interpreter','latex')
xlabel('$h/a$','interpreter','latex');
ylabel('$\lambda_{\mathrm{min}}a^4/\kappa$','interpreter','latex');

end

function [] = kbarcomparison2()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
 
%% prepare a solution as initial guess
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.200001 0.1 0 pi 0 2 0]); %m = 5b
% solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %m = 5
% solinit = bvpinit(linspace(0,1,5),@guess4); %m = 4
%solinit = bvpinit(linspace(0,1,11),@guess3); %m = 3
%solinit = bvpinit(linspace(0,1,9),@guess3b); %m = 3b
%solinit = bvpinit(linspace(0,1,11),@guess2); %m = 2

% kbar = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.27 0.1 0 pi -10 2 0]); %m = 1
% kbar = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.24 0.1 0 pi -10 2 0]); %m = 1
X = 1; kbar = 1; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi -10 2 0]); %m = 1
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
gam = y(5,:);
mu = y(7,1)
L = y(8,1)
eta = y(9,1)
% plot(z,r); return;



%% catenoids
figure(); hold on; axis equal; box on;
v1 = linspace(-0.590668192994804,0.590668192994804);
v2 = linspace(-0.521378393180893,0.521378393180893);
b1 = 0.757966837819559;
b2 = 0.255417422349254;
plot(v1,b1*cosh(v1/b1),'g--','linewidth',3); hold on;
plot(v2,b2*cosh(v2/b2),'g--','linewidth',3);

%%%%%%%%%%%%%%% 
%lower branches
%%%%%%%%%%%%%%%

%% calculate shape for various widths
Xvec = [1:-0.02:0.53 0.525:-0.0001:0.5214 0.52138 0.5213784];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
%     if mod(j,20) == 7
    plot(z,r,'k','LineWidth',2);
%     end
    y_guess = y;
end

figure();
plot(Xvec,Fvec,'linewidth',2)
return

% kbar = 1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
% sol = bvp4c(@odesystem,@bcs,solinit,options);
% t = linspace(0,1,N);
% y_guess = deval(sol,t);
% 
% %% calculate shape for various widths
% Xvec2 = [1.04:-0.01:0.2]/2;
% muvec2 = zeros(1,length(Xvec2));
% Fvec2 = zeros(1,length(Xvec2));
% Evec2 = zeros(1,length(Xvec2));
% % hold on; axis equal; box on;
% for j = 1:length(Xvec2)
%     X = Xvec2(j)
%     solinit = bvpinit(linspace(0,1,11),@newguess);
%     sol = bvp4c(@odesystem,@bcs,solinit,options);
%     t = linspace(0,1,N);
%     y = deval(sol,t);
%     r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
%     mu = y(7,1);
%     muvec2(j) = y(7,1);
%     %b = r((N+1)/2); k1 = psi_s((N+1)/2);
%     %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
%     Fvec2(j) = -2*pi*y(9,1);
%     Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
%     hold on;  grid on;
%     if mod(j,20) == 5
%     plot(z,r,'b','LineWidth',2);
%     end
%     y_guess = y;
% end

kbar = 1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.24 0.1 0 pi -10 2 0]); %m = 1
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shape for various widths
% Xvec2b = [1.04:0.001:1.047 1.0473 1.0476 1.048]/2;
Xvec2b = [1.04:0.01:1.17 1.175 1.18 ]/2;
muvec2b = zeros(1,length(Xvec2b));
Fvec2b = zeros(1,length(Xvec2b));
Evec2b = zeros(1,length(Xvec2b));
% hold on; axis equal; box on;
for j = 1:length(Xvec2b)
    X = Xvec2b(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec2b(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2b(j) = -2*pi*y(9,1);
    Evec2b(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    %if mod(j,20) == 5
    plot(z,r,'b','LineWidth',2);
    %end
    y_guess = y;
end



kbar = -1; X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi -10 2 0]); %m = 1
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec3 = [1.04:-0.01:0.2]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
% hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
    plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
end

kbar = -1; X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.9 0.1 0 pi -10 2 0]); %m = 1
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec3b = [1.04:0.01:1.17  1.171:0.001:1.173]/2;
muvec3b = zeros(1,length(Xvec3b));
Fvec3b = zeros(1,length(Xvec3b));
Evec3b = zeros(1,length(Xvec3b));
% hold on; axis equal; box on;
for j = 1:length(Xvec3b)
    X = Xvec3b(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3b(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3b(j) = -2*pi*y(9,1);
    Evec3b(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
%     plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
end

figure()
plot(Xvec,Fvec,'k','linewidth',2); hold on; grid on;
% plot(Xvec2,Fvec2,'b','linewidth',2);
plot(Xvec2b,Fvec2b,'b','linewidth',2);
plot(Xvec3b,Fvec3b,'r','linewidth',2);

return

%%%%%%%%%%%%%%% 
%upper branches
%%%%%%%%%%%%%%%
kbar = 0; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec4 = [1.0553 1.05 1.04:-0.01:0.0]/2;
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
% hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec4(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 7
    P1 = plot(z,r,'k','LineWidth',2);
    end
    y_guess = y;
end

kbar = 1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec5 = [1.04:-0.01:0.0]/2;
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
Evec5 = zeros(1,length(Xvec5));
% hold on; axis equal; box on;
for j = 1:length(Xvec5)
    X = Xvec5(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec5(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec5(j) = -2*pi*y(9,1);
    Evec5(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
    P2 = plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
end

kbar = 1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec5b = [1.04:0.001:1.055 1.0551 1.0552 1.0553]/2;
muvec5b = zeros(1,length(Xvec5b));
Fvec5b = zeros(1,length(Xvec5b));
Evec5b = zeros(1,length(Xvec5b));
% hold on; axis equal; box on;
for j = 1:length(Xvec5b)
    X = Xvec5b(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec5b(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec5b(j) = -2*pi*y(9,1);
    Evec5b(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
%     plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
end

kbar = -1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec6 = [1.04:-0.01:0.0]/2;
muvec6 = zeros(1,length(Xvec6));
Fvec6 = zeros(1,length(Xvec6));
Evec6 = zeros(1,length(Xvec6));
% hold on; axis equal; box on;
for j = 1:length(Xvec6)
    X = Xvec6(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec6(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec6(j) = -2*pi*y(9,1);
    Evec6(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
    P3 = plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
end

kbar = -1; X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
%% calculate shape for various widths
Xvec6b = [1.04:0.001:1.048 1.0485 1.049]/2;
muvec6b = zeros(1,length(Xvec6b));
Fvec6b = zeros(1,length(Xvec6b));
Evec6b = zeros(1,length(Xvec6b));
% hold on; axis equal; box on;
for j = 1:length(Xvec6b)
    X = Xvec6b(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec6b(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec6b(j) = -2*pi*y(9,1);
    Evec6b(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,20) == 5
%     plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
end

plot([0,0],[0 1],'r','linewidth',2)
bcat = 0.825517; v = linspace(-0.527697,0.527697);
plot(v,bcat*cosh(v/bcat),'k:','LineWidth',2)
legend([P2 P1 P3], {'$\bar\kappa/\kappa = 1$','$\bar\kappa/\kappa = 0$','$\bar\kappa/\kappa = -1$'},'interpreter','latex')
xlabel('$z/a$','interpreter','latex'); ylabel('$r/a$','interpreter','latex');
set(gca,'fontsize',18);
set(gcf,'color','w');

figure();
p1 = plot(2*Xvec,Fvec,'k-','linewidth',3); hold on; grid on;
p2 = plot(2*Xvec2,Fvec2,'b-','linewidth',3);
plot(2*Xvec2b,Fvec2b,'b-','linewidth',3);
p3 = plot(2*Xvec3,Fvec3,'r-','linewidth',3);
plot(2*Xvec3b,Fvec3b,'r-','linewidth',3);
p4 = plot(2*Xvec4,Fvec4,'k-','linewidth',3);
p5 = plot(2*Xvec5,Fvec5,'b-','linewidth',3);
plot(2*Xvec5b,Fvec5b,'b-','linewidth',3);
p6 = plot(2*Xvec6,Fvec6,'r-','linewidth',3);
plot(2*Xvec6b,Fvec6b,'r-','linewidth',3);
plot([2*0.527697 2*0.527697],[-250 100],'k:','linewidth',2)
plot([0 0.2],[0 0],'r-','linewidth',3)
xlabel('$h/a$','interpreter','latex'); 
ylabel('$Fa/\kappa$','interpreter','latex'); 
legend([p2 p1 p3], {'$\bar\kappa/\kappa = 1$','$\bar\kappa/\kappa = 0$','$\bar\kappa/\kappa = -1$'},'interpreter','latex')
set(gca,'fontsize',18);
set(gcf,'color','w');
axis([0 1.1 -225 75])

end

function [] = asymmetricshapes()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints

X = 0.40;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi 0 2 0]); %tether
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether10); %upper loop, increase mu by 1 for lower
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric); %lower loop
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric2); %double loop
% X = 0.40;  solinit = bvpinit(linspace(0,1,6),@guess_asymmetric3); %double loop


options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    K = psi_s.*sin(psi)./r;
    y(7,1)
    y(8,1)
    y(9,1)
    indices = find([0 diff(sign(H))]~=0);

% plot(z,r); hold on; plot(z(indices),r(indices),'o'); return;
% plot(t,H); hold on; plot(t(indices),H(indices),'o'); return;

Xvec = [0.45:-0.05:-0.2 ];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(); hold on; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    H = (psi_s + sin(psi)./r)/2;
    Hp = diff(H)./diff(t);
    Hpp = diff(Hp)./diff(t(1:end-1));
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,1) == 0 || j ==length(Xvec)
%         plot(t(1:end-2),Hpp,'LineWidth',2);
    end
    y_guess = y;
   
end



X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether10); %upper kink, increase mu by 1 for lower
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec2 = [0.45:-0.05:0.25 0.24:-0.01:0 -0.05:-0.05:-0.45];
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
for j = 1:length(Xvec2)
    X = Xvec2(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec2(j) = y(7,1);
    H = (psi_s + sin(psi)./r)/2;
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,5) == 1 || j==length(Xvec2)
%     plot(z,r,'b','LineWidth',2);
    end
    y_guess = y;
end

X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric); %lower kink
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec3 = [0.45 0.39:-0.01:0.1 0.095 0.093 0.092];
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
for j = 1:length(Xvec3)
    X = Xvec3(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec3(j) = y(7,1);
    H = (psi_s + sin(psi)./r)/2;
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,5) == 1 || j == length(Xvec3)
%      plot(z,r,'r','LineWidth',2);
    end
    y_guess = y;
end


X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric2); %double kink
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec4 = [0.45:-0.05:0.2 0.1];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
for j = 1:length(Xvec4)
    X = Xvec4(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec4(j) = y(7,1);
    H = (psi_s + sin(psi)./r)/2;
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,2) == 1 || j == length(Xvec4)
    plot(z,r,'g','LineWidth',2);
    end
    y_guess = y;
end

X = 0.40;  solinit = bvpinit(linspace(0,1,6),@guess_asymmetric3); %double kink
options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec5 = [0.4 0.32 0.24 0.16];
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
Evec5 = zeros(1,length(Xvec5));
for j = 1:length(Xvec5)
    X = Xvec5(j)
    solinit = bvpinit(linspace(0,1,7),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec5(j) = y(7,1);
    H = (psi_s + sin(psi)./r)/2;
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec5(j) = -2*pi*y(9,1);
    Evec5(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,1) == 0 || j == length(Xvec5)
    plot(z,r,'m','LineWidth',2);
    end
    y_guess = y;
end

figure();
subplot(1,2,1)
plot(Xvec,muvec,'ko-');hold on; grid on;
plot(Xvec2,muvec2,'bo-'); 
plot(Xvec3,muvec3,'ro-');
plot(Xvec4,muvec4,'go-');
plot(Xvec5,muvec5,'mo-');
set(gca,'fontsize',18);

subplot(1,2,2)
plot(Xvec,Fvec,'ko-');hold on; grid on;
plot(Xvec2,Fvec2,'bo-'); 
plot(Xvec3,Fvec3,'ro-');
plot(Xvec4,Fvec4,'go-');
plot(Xvec5,Fvec5,'mo-');
set(gca,'fontsize',18);
set(gcf,'color','w');



end

function [] = asymmetricshapes2()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints

% X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5); %n=2
% X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric5); %n=2
% X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric6); %n=2
% X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric7); %n=2, mu=-33 for other side
% X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric8); %n=3
%  X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric9); %n=3 ddd
%  X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric10); %n=3ddd,mu = -84
%  X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric11); %n=3, udd
%  X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric12); %n=1
%  X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric13); %n=2
%  X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric14); %n=3, dud. mu = 38
 X = 0.1; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric15); %n=
% X = 0.40;  solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 1 0.1 0 pi 0 2 0]); %tether
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether10); %upper loop, increase mu by 1 for lower
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric); %lower loop
% X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess_asymmetric2); %double loop
% X = 0.40;  solinit = bvpinit(linspace(0,1,6),@guess_asymmetric3); %double loop


options = bvpset('RelTol',1e-6);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    K = psi_s.*sin(psi)./r;
    y(7,1)
    y(8,1)
    y(9,1)
    indices = find([0 diff(sign(H))]~=0);

plot(z,r,'linewidth',2); hold on; grid on; box on; %plot(z(indices),r(indices),'o'); 
set(gcf,'color','w'); set(gca,'fontsize',18); axis equal; return;
% plot(t,H); hold on; plot(t(indices),H(indices),'o'); return;

Xvec = [0.09 0.08 0.07 0.06];
% Xvec = [0.11 0.16 0.20 0.25 0.28 0.30 0.33 0.35 0.40];
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure(); hold on; box on; axis equal;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,8),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1); gam = y(5,:);
    muvec(j) = y(7,1);
    H = (psi_s + sin(psi)./r)/2;
    Hp = diff(H)./diff(t);
    Hpp = diff(Hp)./diff(t(1:end-1));
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1)
    Falt = L*trapz(t,kappa*2*H.*psi_s + kbar*K + gam.*cos(psi) - y(9,1).*sin(psi))
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,1) == 0 || j ==length(Xvec)
        plot(z,r,'LineWidth',2);
    end
    y_guess = y;
   
end
figure();
plot(Xvec,Fvec,'o-')

end

%% stability eigenvalues by discretization and by bvp4c
function [eval,efun] = evmatrix(eqsol)
global kappa; %kappa = 1; %<---- dont change
global kbar; %kbar = 0;
% global A; A = 2*pi*1; %surface area
% global X; X = 0.5276; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
% global y_guess;
global N; %N = 101; %number of gridpoints
 
%% given an equilibrium solution y_guess, do stability analysis
t = linspace(0,1,N);
y = deval(eqsol,t);
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
gam = y(5,:);
mu = y(7,1);
L = y(8,1);
eta = y(9,1);
psi_ss = cos(psi).*sin(psi)./r.^2 - psi_s.*cos(psi)./r + gam./r.*sin(psi) + eta./r.*cos(psi);
psi_sss = -eta*cos(psi).^2./r.^2 - 2*cos(psi).^2.*sin(psi)./r.^3 + mu*sin(psi)./r - sin(psi).^3./r.^3 ...
          -cos(psi).*sin(psi).*gam./r.^2 + 2*cos(psi).^2.*psi_s./r.^2 - eta*sin(psi).*psi_s./r ...
          -sin(psi).^2.*psi_s./r.^2 + cos(psi).*gam.*psi_s./r + 3*sin(psi).*psi_s.^2/2./r - cos(psi).*psi_ss./r;
H = -(psi_s + sin(psi)./r)/2;
K = psi_s.*sin(psi)./r;
psi0 = psi(1);
psiL = psi(end);

% %work with s as opposed to t
s = L*t; 
h = L/(N-1);
w = r*h; w(1) = w(1)/2; w(end) = w(end)/2; %trapezoid rule weights (measure dA)
Pmat = eye(N) - (2*H)'*(w.*(2*H))/trapz(s,r.*(2*H).^2); %projection onto (2H)^perp

AA = 1/2;
%BB = -(1./(2*r).*(cos(psi).^2./r + sin(psi).*psi_s) + 2*K + mu/2 + H.*(2*psi_s-H));
BB = -mu/2 - cos(psi).^2/2./r.^2 + sin(psi).^2/4./r.^2 - sin(psi).*psi_s./r + 5/4*psi_s.^2;
% old G (Eqn. 129)
% G = cos(psi).^2.*sin(psi).^2./r.^4 - cos(psi).^2.*sin(psi).*psi_s./r.^3 + sin(psi).^3.*psi_s./(2*r.^3) ...
%     -sin(psi).*psi_s.^3./(2*r) + 2*cos(psi).*psi_s.*psi_ss./r + 1.5*psi_ss.^2 - sin(psi).*psi_sss./(2*r) + 1.5*psi_s.*psi_sss;
%CC = (mu*K + 2*(H.^2-K).*(4*H.^2-K) + G);
CC = cos(psi).^2.*sin(psi).^2/2./r.^4 + sin(psi).^4/2./r.^4 - 3/4*cos(psi).^2.*sin(psi).*psi_s./r.^3 ...
     + mu*sin(psi).*psi_s./r - sin(psi).^3.*psi_s/4./r.^3 + cos(psi).^2.*psi_s.^2/4./r.^2 - sin(psi).^2.*psi_s.^2/4./r.^2 ...
     - sin(psi).*psi_s.^3./r + psi_s.^4/2 + cos(psi).*sin(psi).*psi_ss/4./r.^2 + 2*cos(psi).*psi_s.*psi_ss./r ...
     + 3/2*psi_ss.^2 - sin(psi).*psi_sss/2./r + 3/2*psi_s.*psi_sss;

%% (rAu'')'' term, pentadiagonal (which includes BC contributions)
L2 = 4*diag(r(2:end-1)) + diag(r(1:end-2)) + diag(r(3:end)) ...
    -2*diag(r(2:end-2),-1) -2*diag(r(3:end-1),-1) -2*diag(r(2:end-2),1) -2*diag(r(3:end-1),1) ...
    +diag(r(3:end-2),2) + diag(r(3:end-2),-2);
C1 = -(1 + (1+kbar)*h/2*cos(psi(1)))/(1 - (1+kbar)*h/2*cos(psi(1)));   %u_{-1} = C*u_1
Cend = -(1 - (1+kbar)*h/2*cos(psi(end)))/(1 + (1+kbar)*h/2*cos(psi(end))); %u_{N+1} = C*u_{N-1}
L2(1,1) = L2(1,1) + C1*1;
L2(end,end) = L2(end,end) + Cend*1;
L2 = AA*L2/h^4;

%% (rBu')' term, tridiagonal
rB = r.*BB;
rBavg = (rB + circshift(rB,1))/2; rBavg(1) = [];
L1 = diag(rBavg(2:end-1),1) + diag(rBavg(2:end-1),-1) - diag(rB(2:end-1)) - 0.5*diag(rB(1:end-2)) - 0.5*diag(rB(3:end));
L1 = L1/h^2;

%% rCu term, diagonal
L0 = (r(2:end-1).*CC(2:end-1))'.*eye(N-2);

%% compute P'*L*P, which is equivalent to P*L
Lmat = (1./r(2:end-1)').*(L2 + L1 + L0);
M = Pmat(2:end-1,2:end-1)'*Lmat*Pmat(2:end-1,2:end-1);

[vecs,vals] = eig(M);
EV = diag(vals);
[eval,pos] = min(EV);
if abs(eval) < 8
    EV(pos)=[]; vecs(:,pos) = [];
    [eval,pos] = min(EV);
end
% [eval,pos] = min(EV);
% EV(pos)=[]; vecs(:,pos) = [];
% [eval,pos] = min(EV);
% if abs(eval) < 1e-6
%     EV(pos)=[]; vecs(:,pos) = [];
%     [eval,pos] = min(EV);
% end
% eval
efun = [0 vecs(:,pos)' 0];
duds1 = efun(2)/h/L;
efun = efun/duds1; %normalize so slope is 1 at edge
% plot(s,efun); grid on;
p_approx = 0;


return;

%% using approximate eigenfunction, solve ODEs
%% centered finite difference matrices 
D = (diag(ones(1,N-1),1) -diag(ones(1,N-1),-1))/(2*h);
D2 =  (diag(-2*ones(1,N)) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1))/h^2;
D3 = (-0.5*diag(ones(1,N-2),-2) + diag(ones(1,N-1),-1) - diag(ones(1,N-1),1) + 0.5*diag(ones(1,N-2),2))/h^3;
% D4 = (diag(ones(1,N-2),-2) -4*diag(ones(1,N-1),-1) +6*diag(ones(1,N)) -4*diag(ones(1,N-1),1) + diag(ones(1,N-2),2))/h^4;

%% boundary conditions
%  psi = pi/2*ones(1,N); psi_s = zeros(1,N); psi_ss = zeros(1,N); psi_sss = zeros(1,N);
C1 = -(1 + h/2*cos(psi(1)))/(1 - h/2*cos(psi(1)));   %u_{-1} = C*u_1
Cend = -(1 - h/2*cos(psi(end)))/(1 + h/2*cos(psi(end))); %u_{N+1} = C*u_{N-1}
B1 = (1 + h*cos(psi(1)))/(-1 + h*cos(psi(1))); %u_{-2} = B*u_2
Bend = (-1 + h*cos(psi(end)))/(1 + h*cos(psi(end))); %u_{N+2} = B*u_{N-2}
% no torque conditions:
D(1,2) = (1 - C1)/(2*h); D(end,end-1) = (-1 + Cend)/(2*h);
D2(1,2) = (1 + C1)/h^2; D2(end,end-1) = (1 + Cend)/h^2;
D3(1,2) = (-1 + C1)/h^3; D3(end,end-1) = (1 - Cend)/h^3;
D3(1,3) = (0.5 - 0.5*B1)/h^3; D3(end,end-2) = (-0.5 + 0.5*Bend)/h^3;
D3(2,2) = -0.5*C1/h^3; D3(end-1,end-1) = 0.5*Cend/h^3;

lambda = eval; %initial guess for eigenvalue
solinit = bvpinit(linspace(0,1,11),@(t) stabguess(t,efun,(D*efun')',(D2*efun')',(D3*efun')',p_approx,cumtrapz(w,(2*H).*efun),1/(N-1)),lambda);
options = bvpset('RelTol',1e-6);
sol = bvp4c(@(t,y,lambda) ams_ode_new(t,y,lambda,eqsol,mu,L),@(ya,yb,lambda) ams_bc_new(ya,yb,lambda,psi0,psiL),solinit,options);
yperturb = deval(sol,t); 
u = yperturb(1,:); %perturbation in normal direction
eval = sol.parameters;
efun = u;
%plot(s,u);





end

%% u(1) = u, u(2) = u_s, u(3) = Lap(u), u(4) = (Lap(u))_s, u(5) = p, u(6) = integral(r*2H*u)
function out = stabguess(t,u1,u2,u3,u4,u5,u6,h)
out = [interp1(0:h:1,u1,t);
     interp1(0:h:1,u2,t);
     interp1(0:h:1,u3,t);
     interp1(0:h:1,u4,t);
     u5
     interp1(0:h:1,u6,t);
    ];
end

%% area = 2*pi/3
function [] = otherareas1()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi/3; %surface area
global X; X = 0.334/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints

xcat=fsolve(@catarea,[0.1,1])
bcat = xcat(1);
hmax = xcat(2);
x = linspace(-hmax,hmax);
plot(x,bcat*cosh(x/bcat),'k'); hold on; axis equal;

%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 0.3 0]);
%solinit = bvpinit(linspace(0,1,7),@guess2_oa1); %ALL GUESSES USE 7 POINTS
%solinit = bvpinit(linspace(0,1,7),@guess2_oa1);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
mu = y(7,1)
plot(z,r); hold on;

%% calculate shape for various widths
Xvec = [0.3 0.167*2 0.1674*2 0.16745*2 0.167457*2 0.16745716*2]/2; %(2:-0.1:0.1)/2; 
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure; hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec(j) = y(7,1)
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    plot(z,r,'b','LineWidth',2);
    y_guess = y;
end

figure();
subplot(1,3,1)
plot([Xvec],[Fvec],'bo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18);
 
subplot(1,3,2)
plot([Xvec],[ muvec],'bo-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
%plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
%plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
%plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
%plot([0.527697 0.527697],[-50 200],'k--');
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

end

%% initial guess (m = 2)
function out = guess2_oa1(x)
global X;

    %m = 3;
    b = 0.98; %b = 0.825517;
    mu = -361; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 5)
function out = guess5_oa1(x)
global X;

    %m = 3;
    b = 0.98; %b = 0.825517;
    mu = -710; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_oa1(x)
global X;

    %m = 3;
    b = 0.98; %b = 0.825517;
    mu = -350; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 6)
function out = guess6_oa1(x)
global X;

    %m = 3;
    b = 0.98; %b = 0.825517;
    mu = -999; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.001*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 3)
function out = guess3_oa1(x)
global X;

    %m = 3;
    b = 0.98; %b = 0.825517;
    mu = -1010; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.001*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% area = 2*pi*2/3
function [] = otherareas2()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*2/3; %surface area
global X; X = 0.678/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints

xcat=fsolve(@catarea,[0.1,1])
bcat = xcat(1);
hmax = xcat(2);
x = linspace(-hmax,hmax);
plot(x,bcat*cosh(x/bcat),'k'); hold on; axis equal;

%% prepare a solution as initial guess
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 0.6 0]);
%solinit = bvpinit(linspace(0,1,7),@guess2_oa1); %ALL GUESSES USE 7 POINTS
solinit = bvpinit(linspace(0,1,7),@guess3_oa2);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
mu = y(7,1)
plot(z,r); hold on;

%% calculate shape for various widths
Xvec = [0.65 0.68 0.6802]/2; %(2:-0.1:0.1)/2; 
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure; hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec(j) = y(7,1)
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    plot(z,r,'b','LineWidth',2);
    y_guess = y;
end

figure();
subplot(1,3,1)
plot([Xvec],[Fvec],'bo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18);
 
subplot(1,3,2)
plot([Xvec],[ muvec],'bo-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
%plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
%plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
%plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
%plot([0.527697 0.527697],[-50 200],'k--');
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

end

%% initial guess (m = 2)
function out = guess2_oa2(x)
global X;

    %m = 3;
    b = 0.937; %b = 0.825517;
    mu = -185; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 3)
function out = guess3_oa2(x)
global X;

    %m = 3;
    b = 0.9; %b = 0.825517;
    mu = -190; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(7*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 5)
function out = guess5_oa2(x)
global X;

    %m = 3;
    b = 0.9; %b = 0.825517;
    mu = -500; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(5*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_oa2(x)
global X;

    %m = 3;
    b = 0.9; %b = 0.825517;
    mu = -500.01; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(5*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% area = 2*pi*1.2
function [] = otherareas3()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.199; %surface area
global X; X = 0.66;
%global X; X = 0.662;
global y_guess;
global N; N = 101; %number of gridpoints

xcat=fsolve(@catarea,[0.1,1])
bcat = xcat(1);
hmax = xcat(2);
x = linspace(-hmax,hmax);
plot(x,bcat*cosh(x/bcat),'k'); hold on; axis equal;

%% prepare a solution as initial guess
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.3 0.1 0 pi 0 1.3 0])
solinit = bvpinit(linspace(0,1,7),@guess4_oa3); %1 < m < 5, need X = 0.66
%solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.1 0.1 0 pi 0 1.3 0]); %m =5
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
mu = y(7,1)
plot(z,r); hold on;

%% calculate shape for various widths
Xvec = [0.65 0.66 0.662 0.66214]; %(2:-0.1:0.1)/2; 
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
figure; hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec(j) = y(7,1)
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    plot(z,r,'b','LineWidth',2);
    y_guess = y;
end

figure();
subplot(1,3,1)
plot([Xvec],[Fvec],'bo-','LineWidth',2); hold on;
title('Force vs. extension')
xlabel('Half-width z'); ylabel('Force F');  grid on; box on;
set(gca,'FontSize',18);
 
subplot(1,3,2)
plot([Xvec],[ muvec],'bo-','LineWidth',2); hold on;
title('Tension vs. extension')
xlabel('Half-width z'); ylabel('Tension \mu');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

subplot(1,3,3)
plot(Xvec,Evec,'bo-','LineWidth',2); hold on;
%plot(Xvec(1:end-1),diff(Evec)./diff(Xvec)/2); hold on; %factor of 2
%plot(Xvec2(1:end-1),diff(Evec2)./diff(Xvec2)/2); %due to extension/2
%plot(Xvec3(end),muvec3(end),'kp','MarkerSize',18,'MarkerFaceColor','k')
%plot([0.527697 0.527697],[-50 200],'k--');
title('Energy vs. extension')
xlabel('Half-width z'); ylabel('Energy E');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');

end

%% initial guess (m = 2)
function out = guess2_oa3(x)
global X;

    %m = 3;
    b = 0.5; %b = 0.825517;
    mu = -12.5; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(3*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 3)
function out = guess3_oa3(x)
global X;

    %m = 3;
    b = 0.51; %b = 0.825517;
    mu = -29.9; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(7*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_oa3(x)
global X;

    %m = 3;
    b = 0.5; %b = 0.825517;
    mu = -94; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(4*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 5)
function out = guess5_oa3(x)
global X;

    %m = 3;
    b = 0.55; %b = 0.825517;
    mu = -125; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(5*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2)
function out = guess2_tether(x)
global X;

    %m = 3;
    b = 0.6; %b = 0.825517;
    mu = 20; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2)
function out = guess2_tether2(x)
global X;

    %m = 3;
    b = 0.6; %b = 0.825517;
    mu = -20; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether(x)
global X;

    %m = 3;
    b = 0.60; %b = 0.825517;
    mu = -69; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.02*sin(4*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2)
function out = guess2_tether3(x)
global X;

    %m = 3;
    b = 0.9; %b = 0.825517;
    mu = -509; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(4*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2)
function out = guess2_tether4(x)
global X;

    %m = 3;
    b = 0.9; %b = 0.825517;
    mu = -6; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(4*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2)
function out = guess2_tether5(x)
global X;

    %m = 3;
    b = 0.7; %b = 0.825517;
    mu = -5; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 1*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2)
function out = guess2_tether6(x)
global X;

    %m = 3;
    b = 0.8; %b = 0.825517;
    mu = 15; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether2(x)
global X;



    b = 0.9; %b = 0.825517;
    mu = -343; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x/2);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
    
        
%     b = 0.2; %b = 0.825517;
%     mu = -25; %mu = -67.9574;
%     %xs = sinh(X/b) - x/b;
%     out = [pi - acot(sinh(X/b)-x/b);
%         -1/b/(1+(sinh(X/b) - x/b)^2);
%         b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.09*sin(14*pi*x/2/sqrt(1-b^2));
%         b*asinh(x-b*sinh(X/b))/(-b);
%         mu*x;
%         0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
%         mu;
%         2*sqrt(1-b^2);
%         -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether3(x)
global X;

    %m = 3;
    b = 0.61; %b = 0.825517;
    mu = -62; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.02*sin(4*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether4(x)
global X;

    b = 0.91; %b = 0.825517;
    mu = -28; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2)
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x/2);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether5(x)
global X;

    b = 0.9; %b = 0.825517;
    mu = -519; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.1*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether6(x)
global X;

    b = 0.9; %b = 0.825517;
    mu = -547; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether7(x)
global X;

    b = 0.9; %b = 0.825517;
    mu = -556; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether8(x)
global X;

    b = 0.9;
    mu = -505; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(2*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
    %% sixth mode
% b = 0.9; %b = 0.825517;
%     mu = -615; %mu = -67.9574;
%     %xs = sinh(X/b) - x/b;
%     out = [pi - acot(sinh(X/b)-x/b);
%         -1/b/(1+(sinh(X/b) - x/b)^2);
%         b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
%         b*asinh(x-b*sinh(X/b))/(-b);
%         mu*x;
%         0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
%         mu;
%         2*sqrt(1-b^2);
%         -mu*b];


end

%% initial guess (m = 4)
function out = guess2_tether9(x)
global X;

    b = 0.52; %b = 0.825517;
    mu = -20; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess2_tether10(x)
global X;

    b = 0.52; %b = 0.825517;
    mu = -22; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess2_tether11(x)
global X;

    b = 0.45; %b = 0.825517;
    mu = -21; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.01*sin(2*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess2_tether12(x)
global X;

    b = 0.45; %b = 0.825517;
    mu = 28.999999; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.4*sin(2*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether9(x)
global X;

    b = 0.45; %b = 0.825517;
    mu = -38.4; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.1*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether10(x)
global X;

    b = 0.4;
    mu = -42.1985; %-42.1;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.1*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess2_tether13(x)
global X;

    b = 0.4;
    mu = -40; %-42.1;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess2_tether14(x)
global X;

    b = 0.4;
    mu = -40.2; %-42.1;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether11(x)
global X;

    b = 0.45;
    mu = -22; %-42.1;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether12(x)
global X;

    b = 0.4;
    mu = -42.4; %-42.1;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether13(x)
global X;

    b = 0.4;
    mu = -86; %-42.1;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.15*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess2_tether15(x)
global X;

    b = 0.4;
    mu = -21; %-42.1;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) - 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 2)
function out = guess2_tether_kbar(x)
global X;

    %m = 3;
    b = 0.9; %b = 0.825517;
    mu = -526; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

    %kbar = 0.5;
%     b = 0.9; %b = 0.825517;
%     mu = -522; %mu = -67.9574;
%     %xs = sinh(X/b) - x/b;
%     out = [pi - acot(sinh(X/b)-x/b);
%         -1/b/(1+(sinh(X/b) - x/b)^2);
%         b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
%         b*asinh(x-b*sinh(X/b))/(-b);
%         mu*x;
%         0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
%         mu;
%         2*sqrt(1-b^2);
%         -mu*b];
end

%% initial guess (m = 2)
function out = guess2_tether_kbar2(x)
global X;

    %m = 3;
    b = 0.9; %b = 0.825517;
    mu = -558; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether_kbar(x)
global X;

    b = 0.9; %b = 0.825517;
    mu = -525; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

    %kbar= 1
%     b = 0.9; %b = 0.825517;
%     mu = -525; %mu = -67.9574;
%     %xs = sinh(X/b) - x/b;
%     out = [pi - acot(sinh(X/b)-x/b);
%         -1/b/(1+(sinh(X/b) - x/b)^2);
%         b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
%         b*asinh(x-b*sinh(X/b))/(-b);
%         mu*x;
%         0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
%         mu;
%         2*sqrt(1-b^2);
%         -mu*b];

    %kbar = 0.5
%     b = 0.9; %b = 0.825517;
%     mu = -521; %mu = -67.9574;
%     %xs = sinh(X/b) - x/b;
%     out = [pi - acot(sinh(X/b)-x/b);
%         -1/b/(1+(sinh(X/b) - x/b)^2);
%         b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
%         b*asinh(x-b*sinh(X/b))/(-b);
%         mu*x;
%         0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
%         mu;
%         2*sqrt(1-b^2);
%         -mu*b];
end

%% initial guess (m = 4)
function out = guess4_tether_kbar2(x)
global X;

    %m = 3;
    b = 0.9; %b = 0.825517;
    mu = -526; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

    %kbar = 0.5;
%     b = 0.9; %b = 0.825517;
%     mu = -526; %mu = -67.9574;
%     %xs = sinh(X/b) - x/b;
%     out = [pi - acot(sinh(X/b)-x/b);
%         -1/b/(1+(sinh(X/b) - x/b)^2);
%         b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.01*sin(2*pi*x/2/sqrt(1-b^2));
%         b*asinh(x-b*sinh(X/b))/(-b);
%         mu*x;
%         0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
%         mu;
%         2*sqrt(1-b^2);
%         -mu*b];
end

function out = guess_asymmetric(x)

global X;

    b = 0.52; %b = 0.825517;
    mu = -21; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric2(x)

global X;

    b = 0.52; %b = 0.825517;
    mu = -34; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric3(x)

global X;

    b = 0.52; %b = 0.825517;
    mu = -241; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric4(x)

global X;

    b = 0.52; %b = 0.825517;
    mu = -247; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric5(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =-2; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric6(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =-5; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric7(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =-10; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric8(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =-24; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric9(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =-35; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric10(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =-41; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric11(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =-97; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric12(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =1.5; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric13(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =3; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric14(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =28; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric15(x)

global X;

    b = 0.52; 
    mu = 3; %86; %-46; %mu = -12;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

function out = guess_asymmetric16(x)

global X;

    b = 0.52; %b = 0.825517;
    mu =99; %mu = -67.9574;
    %xs = sinh(X/b) - x/b;
    out = [pi - acot(sinh(X/b)-x/b);
        -1/b/(1+(sinh(X/b) - x/b)^2);
        b*sqrt(1+(sinh(X/b) - x/b)^2) + 0.2*sin(4*pi*x);
        b*asinh(x-b*sinh(X/b))/(-b);
        mu*x;
        0.5*(b*asinh(x/b-sinh(X/b))+(x-b*sinh(X/b))*sqrt(1+(sinh(X/b)-x/b)^2))-0.5*(b*asinh(0/b-sinh(X/b))+(0-b*sinh(X/b))*sqrt(1+(sinh(X/b)-0/b)^2));
        mu;
        2*sqrt(1-b^2);
        -mu*b];

end

%% force vs extension for catenoid soap film
function fvhfig()

%% force vs extension
hvec = linspace(0,1.32543,501);
bvec1 = zeros(size(hvec));
bvec2 = zeros(size(hvec));
for j = 1:length(hvec)
    if j == 1
    bvec1(j) = fzero(@(b) b*cosh(hvec(j)/2/b) - 1,0.5);
    bvec2(j) = fzero(@(b) b*cosh(hvec(end)/2/b) - 1,0.5);
    else
    bvec1(j) = fzero(@(b) b*cosh(hvec(j)/2/b) - 1,bvec1(j-1));
    bvec2(j) = fzero(@(b) b*cosh(hvec(end-j+1)/2/b) - 1,bvec2(j-1));

    end
    
end
bvec2(end) = 0;
bvec2 = fliplr(bvec2);
bmax = (bvec1(end) + bvec2(end) )/2;


figure();
plot([hvec 1.32548],[bvec1 bmax],'b-','linewidth',3); hold on; grid on; box on;
plot([hvec 1.32548],[bvec2 bmax],'g--','linewidth',3);
set(gcf,'color','w')
set(gca,'fontsize',18')
xlabel('$h/a$','interpreter','latex')
ylabel('$F/(2\pi a \mu)$','interpreter','latex')

%% area vs extension
Avec1 = pi*bvec1.*(hvec+bvec1.*sinh(hvec./bvec1))/2/pi;
Avec2 = pi*bvec2.*(hvec+bvec2.*sinh(hvec./bvec2))/2/pi;


figure();
plot(hvec,2*sqrt(1-bvec1.^2),'b-','linewidth',3); hold on; grid on;
plot(hvec,2*sqrt(1-bvec2.^2),'g--','linewidth',3);
set(gcf,'color','w')
set(gca,'fontsize',18')
xlabel('$h/a$','interpreter','latex')
ylabel('$L/a$','interpreter','latex')

figure();
plot([0 1.4],[1 1],'k:','linewidth',3); hold on; grid on; box on;
plot(hvec,Avec1,'b-','linewidth',3);
plot(hvec,Avec2,'g--','linewidth',3); 
set(gcf,'color','w')
set(gca,'fontsize',18')
xlabel('$h/a$','interpreter','latex')
ylabel('$A/(2\pi a^2)$','interpreter','latex')

dAdb1 = diff(Avec1)./diff(bvec1);
dAdb2 = diff(Avec2)./diff(bvec2);
dhdb = (hvec./bvec1 -2./sqrt(1-bvec1.^2));
d2hdb2 = -2*bvec1./(1-bvec1.^2).^1.5 - hvec./bvec1.^2 + dhdb./bvec1;
Lvec1 = 2*sqrt(1-bvec1.^2);
temp = (16*acsch(2*bvec1./Lvec1) -16*bvec1.^4.*atanh(Lvec1/2) - Lvec1.^2.*(2*Lvec1+4*bvec1.^2.*log(1+Lvec1.*(Lvec1+2)./(2*bvec1.^2))) )./Lvec1.^4;

% figure()
% plot(hvec(2:end-1),diff(dAdb1)./diff(bvec1(1:end-1)),'b-','linewidth',3); hold on; grid on;
% % plot(hvec,-4*bvec1./Lvec1,'linewidth',3); hold on; grid on;
% plot(hvec,temp,'g--','linewidth',3)
% set(gcf,'color','w')
% set(gca,'fontsize',18')
% xlabel('$h/a$','interpreter','latex')
% ylabel('$L/a$','interpreter','latex')

end

%% makes figures for A = 2*pi*1
function [] = onecatfig()
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints
 
%% prepare a solution as initial guess
solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]); %m = 1
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); 
gam = y(5,:);
% mu = y(7,1)
% L = y(8,1)
% eta = y(9,1)
% plot(z,r); return;

%% ==== colors ====
col1 = [0, 0.4470, 0.7410];
col2 = [0.8500, 0.3250, 0.0980];
col3 = [0.9290, 0.6940, 0.1250];
col4 = [0.4940, 0.1840, 0.5560];
col5 = [0.4660, 0.6740, 0.1880];

disp('n = 1')
%% calculate shape for various widths
Xvec = [1.054 1.05:-0.01:0.14]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%         plot(z,r,'color',col1,'LineWidth',2);
    end
    if 2*X == 0.7
%         plot(z,r,'color',[col1 0.75],'LineWidth',2);
%         shape1C = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1C = centercrop(shape1C,0.2,0.2); imwrite(shape1C,'shape1C.png');
    elseif 2*X == 0.35
%         plot(z,r,'color',[col1 0.5],'LineWidth',2);
%         shape1B = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1B = centercrop(shape1B,0.2,0.2); imwrite(shape1B,'shape1B.png');
    end
    y_guess = y;
end
 
%% prepare a solution as initial guess
X = 0.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]); %m = 1b
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec2 = [0.2:0.01:1.04 1.045:0.005:1.05 1.054]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec2(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) ==1
%          plot(z,r,'color',col1,'LineWidth',2);
    end
    if 2*X == 0.7
%         plot(z,r,'color',[col1 0.75],'LineWidth',2);
%         shape1E = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1E = centercrop(shape1E,0.2,0.2); imwrite(shape1E,'shape1E.png');
    elseif 2*X == 0.35
%         plot(z,r,'color',[col1 0.5],'LineWidth',2);
%         shape1F = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1F = centercrop(shape1F,0.2,0.2); imwrite(shape1F,'shape1F.png');
    end

   y_guess = y;
end

%% calculate shape for various widths
Xvec3 = [0.2:-0.01:0]/2;
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec3(j) = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) ==1
%          plot(z,r,'color',col1,'LineWidth',2);
    end
    if 2*X == 0.0
%         plot(z,r,'color',[col1 0.4],'LineWidth',2);
%         shape1G = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1G = centercrop(shape1G,0.2,0.2); imwrite(shape1G,'shape1G.png');
    end

   y_guess = y;
end

disp('n = 2');
%% prepare a solution as initial guess
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),@guess2); %m = 2
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec4 = [1.05:-0.01:0]/2;
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec4(j) = y(7,1);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col2,'LineWidth',2);
    end
    if 2*X == 0.0
%         plot(z,r,'color',[col2 0.25],'LineWidth',2);
%         shape2A = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2A = centercrop(shape2A,0.2,0.2); imwrite(shape2A,'shape2A.png');
%         shape2Ap = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Ap = centercrop(shape2Ap,0.2,0.2); imwrite(shape2Ap,'shape2Ap.png');
    elseif abs(2*X- 0.35) < 0.001
%         plot(z,r,'color',[col2 0.5],'LineWidth',2);
%         shape2B = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2B = centercrop(shape2B,0.2,0.2); imwrite(shape2B,'shape2B.png');    
%         shape2Bp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Bp = centercrop(shape2Bp,0.2,0.2); imwrite(shape2Bp,'shape2Bp.png'); 
    elseif 2*X == 0.7
%         plot(z,r,'color',[col2 0.75],'LineWidth',2);
%         shape2C = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2C = centercrop(shape2C,0.2,0.2); imwrite(shape2C,'shape2C.png');
%         shape2Cp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Cp = centercrop(shape2Cp,0.2,0.2); imwrite(shape2Cp,'shape2Cp.png');
    end
    
   y_guess = y;
end


disp('n = 3');
%% prepare a solution as initial guess
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),@guess3); %m = 3
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec5 = [1.05:-0.01:0]/2;
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
Evec5 = zeros(1,length(Xvec5));
for j = 1:length(Xvec5)
    X = Xvec5(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec5(j) = y(7,1);
    Fvec5(j) = -2*pi*y(9,1);
    Evec5(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col3,'LineWidth',2);
    end
    if 2*X == 0.7
%         plot(z,r,'color',[col3 0.75],'LineWidth',2);
        shape3C = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
        shape3C = centercrop(shape3C,0.2,0.2); imwrite(shape3C,'shape3C.tif','tif','resolution',600,'compression','none');       
    elseif abs(2*X- 0.35) < 0.001
%         plot(z,r,'color',[col3 0.5],'LineWidth',2);
%         shape3B = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3B = centercrop(shape3B,0.2,0.2); imwrite(shape3B,'shape3B.png');          
    elseif 2*X == 0
%         plot(z,r,'color',[col3 0.25],'LineWidth',2);
%         shape3A = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3A = centercrop(shape3A,0.2,0.2); imwrite(shape3A,'shape3A.png');         
    end

   y_guess = y;
end

return

%% prepare a solution as initial guess
X = 1.04/2; solinit = bvpinit(linspace(0,1,9),@guess3b); %m = 3b
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec6 = [1.05:-0.01:0]/2;
muvec6 = zeros(1,length(Xvec6));
Fvec6 = zeros(1,length(Xvec6));
Evec6 = zeros(1,length(Xvec6));
for j = 1:length(Xvec6)
    X = Xvec6(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec6(j) = y(7,1);
    Fvec6(j) = -2*pi*y(9,1);
    Evec6(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col3,'LineWidth',2);
    end
    if 2*X == 0.7
        plot(z,r,'color',[col3 0.75],'LineWidth',2);
%         shape3E = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3E = centercrop(shape3E,0.2,0.2); imwrite(shape3E,'shape3E.png');       
    elseif abs(2*X - 0.35) < 0.001
%         plot(z,r,'color',[col3 0.5],'LineWidth',2);
%         shape3F = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3F = centercrop(shape3F,0.2,0.2); imwrite(shape3F,'shape3F.png');         
    elseif 2*X == 0
%         plot(z,r,'color',[col3 0.25],'LineWidth',2);
%         shape3G = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3G = centercrop(shape3G,0.2,0.2); imwrite(shape3G,'shape3G.png');         
    end
   y_guess = y;
end


disp('n = 4')
%% prepare a solution as initial guess
X = 1.04/2; solinit = bvpinit(linspace(0,1,5),@guess4); %m = 4
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec7 = [1.05:-0.01:0]/2;
muvec7 = zeros(1,length(Xvec7));
Fvec7 = zeros(1,length(Xvec7));
Evec7 = zeros(1,length(Xvec7));
for j = 1:length(Xvec7)
    X = Xvec7(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec7(j) = y(7,1);
    Fvec7(j) = -2*pi*y(9,1);
    Evec7(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col4,'LineWidth',2);
    end
    if 2*X == 0.7
%         plot(z,r,'color',[col4 0.75],'LineWidth',2);
%         shape4C = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4C = centercrop(shape4C,0.2,0.2); imwrite(shape4C,'shape4C.png');   
%         shape4Cp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Cp = centercrop(shape4Cp,0.2,0.2); imwrite(shape4Cp,'shape4Cp.png');  
    elseif abs(2*X - 0.35) < 0.001
%         plot(z,r,'color',[col4 0.5],'LineWidth',2);
%         shape4B = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4B = centercrop(shape4B,0.2,0.2); imwrite(shape4B,'shape4B.png');  
%         shape4Bp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Bp = centercrop(shape4Bp,0.2,0.2); imwrite(shape4Bp,'shape4Bp.png');  
    elseif 2*X == 0
%         plot(z,r,'color',[col4 0.25],'LineWidth',2);
%         shape4A = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4A = centercrop(shape4A,0.2,0.2); imwrite(shape4A,'shape4A.png');   
%         shape4Ap = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Ap = centercrop(shape4Ap,0.2,0.2); imwrite(shape4Ap,'shape4Ap.png');  
    end

   y_guess = y;
end


disp('n = 5')
%% prepare a solution as initial guess
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 2 0]); %m = 5
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec8 = [1.05:-0.01:0]/2;
muvec8 = zeros(1,length(Xvec8));
Fvec8 = zeros(1,length(Xvec8));
Evec8 = zeros(1,length(Xvec8));
for j = 1:length(Xvec8)
    X = Xvec8(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec8(j) = y(7,1);
    Fvec8(j) = -2*pi*y(9,1);
    Evec8(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col5,'LineWidth',2);
    end
    if 2*X == 0.7
%         plot(z,r,'color',[col5 0.75],'LineWidth',2);
%         shape5C = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5C = centercrop(shape5C,0.2,0.2); imwrite(shape5C,'shape5C.png');       
    elseif abs(2*X - 0.35) < 0.001
%         plot(z,r,'color',[col5 0.5],'LineWidth',2);
%         shape5B = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5B = centercrop(shape5B,0.2,0.2); imwrite(shape5B,'shape5B.png');         
    elseif 2*X == 0
%         plot(z,r,'color',[col5 0.25],'LineWidth',2);
%         shape5A = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5A = centercrop(shape5A,0.2,0.2); imwrite(shape5A,'shape5A.png');         
    end

   y_guess = y;
end

%% prepare a solution as initial guess
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.200001 0.1 0 pi 0 2 0]); %m = 5b
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec9 = [1.05:-0.01:0]/2;
muvec9 = zeros(1,length(Xvec9));
Fvec9 = zeros(1,length(Xvec9));
Evec9 = zeros(1,length(Xvec9));
for j = 1:length(Xvec9)
    X = Xvec9(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec9(j) = y(7,1);
    Fvec9(j) = -2*pi*y(9,1);
    Evec9(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col5,'LineWidth',2);
    end
    if 2*X == 0.7
        plot(z,r,'color',[col5 0.75],'LineWidth',2);
%         shape5E = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5E = centercrop(shape5E,0.2,0.2); imwrite(shape5E,'shape5E.png');       
    elseif abs(2*X - 0.35) < 0.001
        plot(z,r,'color',[col5 0.5],'LineWidth',2);
%         shape5F = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5F = centercrop(shape5F,0.2,0.2); imwrite(shape5F,'shape5F.png');         
    elseif 2*X == 0
        plot(z,r,'color',[col5 0.25],'LineWidth',2);
%         shape5G = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5G = centercrop(shape5G,0.2,0.2); imwrite(shape5G,'shape5G.png');         
    end

   y_guess = y;
end


%% catenoid solution
bcat = 0.825517; v = linspace(-0.527697,0.527697,101);
r = bcat*cosh(v/bcat);
shape1D = snapshot(0.6,1.35,zeros(size(v)),zeros(size(v)),v,r,v(end),kappa,N);
shape1D = centercrop(shape1D,0.2,0.2); imwrite(shape1D,'shape1D.tif','tif','resolution',600);  

return
% plot(v,bcat*cosh(v/bcat),'k','LineWidth',3)
% % plot(v./(v.*v+bcat^2*cosh(v/bcat).^2)*(v(1)^2+bcat^2*cosh(v(1)/bcat)^2),...
% %     bcat*cosh(v/bcat)./(v.*v+bcat^2*cosh(v/bcat).^2)*(v(1)^2+bcat^2*cosh(v(1)/bcat)^2),'k--');
% xlabel('z'); ylabel('r');
% title('Shapes');

%% disk solution
fig = figure;
[T,P] = meshgrid([0,1],linspace(0,2*pi,60));
r = [zeros(60,1) ones(60,1)]; z = zeros(60,2); C = zeros(60,2);
s=surf(z,r.*cos(P),r.*sin(P),C); hold on;
%         s=surf(z(T),r(T).*cos(P),r(T).*sin(P),log(C)/log(10)); hold on;
% s=surf(z(T),r(T).*cos(P),r(T).*sin(P),C); hold on;

%         cb = colorbar('southoutside'); cb.Label.String = 'log_{10} (\kappa/2)(2H)^2';
%         caxis([-2.0 2.0])
colormap('winter')
caxis([-20 20])
s.EdgeColor = 'none';
X = 0;
plot3(ones([1,100])*X,cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',4);
plot3(-ones([1,100])*X,cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',4);
axis equal;
axis([-0.6 0.6 -1.35   1.35 -1.35 1.35]); %box on;
%         axis([-width width -height height -height height]); %box on;
camorbit(15,0);
axis off;
set(gcf,'color','w');set(gca,'FontSize',18);
for q = 1:2:60
    plot3([0 0],[0,1*cos(2*pi*q/60)],[0,1*sin(2*pi*q/60)],'k')
end

frame = getframe(fig);
shape1A = frame2im(frame);
shape1A = centercrop(shape1A,0.2,0.2); imwrite(shape1A,'shape1A.png');


%% force and tension plots

set(gca,'FontSize',18);
set(gcf,'color','w');
close all;

figure();
lw = 3;
% subplot(1,2,1)
subplot('position',[0.1 0.1 0.25 0.85])
Xvecu = Xvec(Xvec < 0.08); Xvecs = Xvec(Xvec >= 0.08);
Fvecu = Fvec(Xvec < 0.08); Fvecs = Fvec(Xvec >= 0.08);
% plot([2*Xvec 0],[Fvec 0],'-','linewidth',lw,'color',col1); hold on;
plot([2*Xvecu 0],[Fvecu 0],':','linewidth',lw,'color',col1); hold on;
plot([2*Xvecs ],[Fvecs ],'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec2,Fvec2,'-','linewidth',lw,'color',col1); hold on;
% plot(2*Xvec3,Fvec3,'-','linewidth',lw,'color',col1); hold on;
Xvec3u = Xvec3(Xvec3 < 0.07); Xvec3s = Xvec3(Xvec3 >= 0.07);
Fvec3u = Fvec3(Xvec3 < 0.07); Fvec3s = Fvec3(Xvec3 >= 0.07);
plot(2*Xvec3s,Fvec3s,'-','linewidth',lw,'color',col1); 
plot(2*Xvec3u,Fvec3u,':','linewidth',lw,'color',col1); 
plot(2*Xvec4,Fvec4,':','linewidth',lw,'color',col2); hold on;
plot(2*Xvec5,Fvec5,':','linewidth',lw,'color',col3); hold on;
plot(2*Xvec6,Fvec6,':','linewidth',lw,'color',col3); hold on;
plot(2*Xvec7,Fvec7,':','linewidth',lw,'color',col4); hold on;
plot(2*Xvec8,Fvec8,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec9,Fvec9,':','linewidth',lw,'color',col5); hold on;
plot([1.055 1.055],[-3200 0],'--','color',[0.6 0.6 0.6],'linewidth',2)
hh = text(1.027,-2950,'$h=h^*$','color',[0.6 0.6 0.6],'interpreter','latex','fontsize',18); set(hh,'rotation',90);
% title('Force vs. extension')
xlabel('$h/a$','interpreter','latex'); ylabel('$Fa/\kappa$','interpreter','latex');  grid on; box on;
set(gca,'FontSize',18);
xlim([0 0.55*2])
ylim([-3200 0]) 
% legend('n = 1','n = 2','n = 3','n = 4','n = 5');
plotpoint(0,0,'1A',col1,pi/2,0.08,50);
plotpoint(0.35,-0.01,'1B',col1,pi/2,0.08,50);
plotpoint(0.70,-1.57,'1C',col1,pi/2,0.08,50);
plotpoint(1.055,-29.8,'1D',col1,0,0.04,50);
plotpoint(0.70,-115.3,'1E',col1,-pi/2,0.08,50);
plotpoint(0.35,-161.1,'1F',col1,-pi/2,0.08,50);
plotpoint(0,-189.5,'1G',col1,-pi/2,0.08,50);
plotpoint(0,-405.7,'2A',col2,-pi/2,0.08,50);
plotpoint(0.35,-386.9,'2B',col2,-pi/2,0.08,50);
plotpoint(0.7,-292.6,'2C',col2,-pi/2,0.08,50);
plotpoint(1.055,-151.4,'2D',col2,0,0.04,50);
plotpoint(0,-1130,'3A',col3,pi/2,0.08,50);
plotpoint(0.35,-828.9,'3B',col3,pi/2,0.08,50);
plotpoint(0.7,-601.1,'3C',col3,pi/2,0.08,50);
plotpoint(1.055,-352.4,'3D',col3,0,0.04,50);
plotpoint(0.7,-698.4,'3E',col3,-pi/2,0.08,50);
plotpoint(0.35,-930.7,'3F',col3,-pi/2,0.08,50);
plotpoint(0,-1155,'3G',col3,-pi/2,0.08,50);
plotpoint(0,-1964,'4A',col4,-pi/3,0.08,50);
plotpoint(0.35,-1554,'4B',col4,-pi/2.5,0.08,50);
plotpoint(0.7,-1153,'4C',col4,-pi/2.5,0.08,50);
plotpoint(1.055,-633.8,'4D',col4,0,0.04,50);
plotpoint(0,-3153,'5A',col5,pi/2,0.08,50);
plotpoint(0.35,-2376,'5B',col5,pi/2,0.08,50);
plotpoint(0.7,-1748,'5C',col5,pi/2,0.08,50);
plotpoint(1.055,-995.8,'5D',col5,0,0.04,50);
plotpoint(0.7,-1846,'5E',col5,-pi/2,0.08,50);
plotpoint(0.35,-2477,'5F',col5,-pi/2,0.08,50);
plotpoint(0,-3198,'5G',col5,pi/12,0.06,50);


% subplot(1,2,2)
subplot('position',[0.5 0.1 0.25 0.85])
plot(2*Xvec,muvec,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec2,muvec2,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec3,muvec3,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec4,muvec4,'-','linewidth',lw,'color',col2); hold on;
plot(2*Xvec5,muvec5,'-','linewidth',lw,'color',col3); hold on;
plot(2*Xvec6,muvec6,'-','linewidth',lw,'color',col3); hold on;
plot(2*Xvec7,muvec7,'-','linewidth',lw,'color',col4); hold on;
plot(2*Xvec8,muvec8,'-','linewidth',lw,'color',col5); hold on;
plot(2*Xvec9,muvec9,'-','linewidth',lw,'color',col5); hold on;
% title('Tension vs. extension')
xlabel('$h/a$','interpreter','latex'); ylabel('$\mu a^2/\kappa$','interpreter','latex');  grid on; box on;
set(gca,'FontSize',28);
set(gcf,'color','w');
xlim([0 0.53*2])
ylim([-200 150])
% plotpoint(0,0,'1A',col1,pi/2);
% plotpoint2(0.35,35.04,'1B',col1,-pi/2);
plotpoint2(0.70,7.373,'1C',col1,-pi/2);
plotpoint2(1.055,-5.75,'1D',col1,-pi/2);
plotpoint2(0.70,-9.474,'1E',col1,-pi/4);
% plotpoint2(0.35,-0.747,'1F',col1,-pi/2);
plotpoint2(0,14.63,'1G',col1,-pi/2);
plotpoint2(0,59.26,'2A',col2,-pi/2);
% plotpoint2(0.35,9.149,'2B',col2,-5*pi/6);
plotpoint2(0.7,-19.79,'2C',col2,-pi/2);
plotpoint2(1.055,-29.3,'2D',col2,-pi/2);
plotpoint2(0,120.6,'3A',col3,-pi/2);
% plotpoint2(0.35,22.2,'3B',col3,5*pi/6);
plotpoint2(0.7,-38.1,'3C',col3,pi/6);
plotpoint2(1.055,-68.5,'3D',col3,-pi/2);
plotpoint2(0.7,-48.93,'3E',col3,-pi/2);
% plotpoint2(0.35,18.04,'3F',col3,-5*pi/6);
plotpoint2(0,129.3,'3G',col3,pi/2);
% plotpoint2(0,222.5,'4A',col4,-pi/2);
% plotpoint2(0.35,36.76,'4B',col4,pi/2);
plotpoint2(0.7,-77.12,'4C',col4,-pi/2);
plotpoint2(1.055,-122.3,'4D',col4,-pi/2);
% plotpoint2(0,345,'5A',col5,-pi/2);
% plotpoint2(0.35,59.81,'5B',col5,pi/6);
plotpoint2(0.7,-114.5,'5C',col5,pi/4);
plotpoint2(1.055,-192.7,'5D',col5,-pi/2);
plotpoint2(0.7,-125.6,'5E',col5,-pi/(4/3));
% plotpoint2(0.35,55.58,'5F',col5,-pi/6);
plotpoint2(0,352.3,'5G',col5,-pi/2);

plot([0.33, 0.37, 0.37, 0.33, 0.33],[-1,-1,66,66,-1],'k--','linewidth',2)


% zoomed in portion
axes('position',[.52 .15 .05 .3]); box on;
indexOfInterest = (2*Xvec >= 0.34) & (2*Xvec <= 0.36);
indexOfInterest2 = (2*Xvec2 >= 0.34) & (2*Xvec2 <= 0.36);
% indexOfInterest3 = (2*Xvec3 >= 0.34) & (2*Xvec3 <= 0.36);
indexOfInterest4 = (2*Xvec4 >= 0.34) & (2*Xvec4 <= 0.36);
indexOfInterest5 = (2*Xvec5 >= 0.34) & (2*Xvec5 <= 0.36);
indexOfInterest6 = (2*Xvec6 >= 0.34) & (2*Xvec6 <= 0.36);
indexOfInterest7 = (2*Xvec7 >= 0.34) & (2*Xvec7 <= 0.36);
indexOfInterest8 = (2*Xvec8 >= 0.34) & (2*Xvec8 <= 0.36);
indexOfInterest9 = (2*Xvec9 >= 0.34) & (2*Xvec9 <= 0.36);
plot(2*Xvec(indexOfInterest),muvec(indexOfInterest),'color',col1,'linewidth',lw); hold on;
plot(2*Xvec2(indexOfInterest2),muvec2(indexOfInterest2),'color',col1,'linewidth',lw);
% plot(2*Xvec3(indexOfInterest3),muvec3(indexOfInterest3),'color',col1,'linewidth',lw);
plot(2*Xvec4(indexOfInterest4),muvec4(indexOfInterest4),'color',col2,'linewidth',lw);
plot(2*Xvec5(indexOfInterest5),muvec5(indexOfInterest5),'color',col3,'linewidth',lw);
plot(2*Xvec6(indexOfInterest6),muvec6(indexOfInterest6),'color',col3,'linewidth',lw);
plot(2*Xvec7(indexOfInterest7),muvec7(indexOfInterest7),'color',col4,'linewidth',lw);
plot(2*Xvec8(indexOfInterest8),muvec8(indexOfInterest8),'color',col5,'linewidth',lw);
plot(2*Xvec9(indexOfInterest9),muvec9(indexOfInterest9),'color',col5,'linewidth',lw);
plotpoint(0.35,35.04,'1B',col1,-pi/2,0.08,3);
plotpoint(0.35,-0.747,'1F',col1,pi/2,0.08,3);
plotpoint(0.35,9.149,'2B',col2,-pi/2,0.08,3);
plotpoint(0.35,22.2,'3B',col3,pi/2,0.08,3);
plotpoint(0.35,18.04,'3F',col3,-pi/2,0.08,3);
plotpoint(0.35,36.76,'4B',col4,pi/2,0.08,3);
plotpoint(0.35,59.81,'5B',col5,pi/2,0.08,3);
plotpoint(0.35,55.58,'5F',col5,-pi/2,0.08,3);

axis tight;

end

%% makes figures for A = 2*pi*1.1
function [] = twocatfig()
global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.1; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 101; %number of gridpoints

%% ==== colors ====
col1 = [0, 0.4470, 0.7410];
col2 = [0.8500, 0.3250, 0.0980];
col3 = [0.9290, 0.6940, 0.1250];
col4 = [0.4940, 0.1840, 0.5560];
col5 = [0.4660, 0.6740, 0.1880];

%% prepare a solution as initial guess (thick catenoid branch)
X = 0.2/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t); y = y_guess;
r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); gam = y(5,:);
% plot(z,r); return;

disp('n = 1')
%% calculate shape for various widths
Xvec = [1:0.01:0.57*2 0.585*2 0.59*2 0.5901*2]/2;
muvec = zeros(1,length(Xvec));
Fvec = zeros(1,length(Xvec));
Evec = zeros(1,length(Xvec));
% hold on; axis equal; box on;
for j = 1:length(Xvec)
    X = Xvec(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec(j) = y(7,1);
    Fvec(j) = -2*pi*y(9,1);
    Evec(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%         plot(z,r,'color',col1,'LineWidth',2);
    end
    if abs(2*X - 0.5901*2) < 0.00001
%         plot(z,r,'color',[col1 4/7],'LineWidth',2);
%         shape1K = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1K = centercrop(shape1K,0.2,0.2); imwrite(shape1K,'shape1K.png');
    end

    y_guess = y;
end

%% prepare a solution as initial guess (thick catenoid, reverse)
X = 0.2/2; solinit = bvpinit(linspace(0,1,11),[pi/2 -0.01 1 0.1 0 pi 0 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec2 = [1:-0.01:0]/2;
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
% hold on; axis equal; box on;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    mu = y(7,1);
    muvec2(j) = y(7,1);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%         plot(z,r,'color',col1,'LineWidth',2);
    end

    if abs(2*X - 0) < 0.001
%         plot(z,r,'color',[col1 1/7],'LineWidth',2);
%         shape1H = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1H = centercrop(shape1H,0.2,0.2); imwrite(shape1H,'shape1H.png');
    end
    if abs(X - 0.35/2) < 0.001
%         plot(z,r,'color',[col1 2/7],'LineWidth',2);
%         shape1I = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1I = centercrop(shape1I,0.2,0.2); imwrite(shape1I,'shape1I.png');
    end
    if abs(X - 0.7/2) < 0.001
%         plot(z,r,'color',[col1 3/7],'LineWidth',2);
%         shape1J = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape1J = centercrop(shape1J,0.2,0.2); imwrite(shape1J,'shape1J.png');
    end
    
    y_guess = y;
end

 
%% prepare a solution as initial guess (thick cat to thin cat)
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.19 0.1 0 pi -10 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);
 
%% calculate shape for various widths
Xvec3 = [0.5901 0.59:-0.005:0.525 0.522];
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec3(j) = y(7,1);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,2) ==1
%         plot(z,r,'color',col1,'LineWidth',2);
    end

    if abs(X - 0.522) < 0.0001
%         plot(z,r,'color',[col1 5/7],'LineWidth',2);
%        shape1l = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape1l = centercrop(shape1l,0.2,0.2); imwrite(shape1l,'shape1lt.png');
    end
   y_guess = y;
end

%% prepare a solution as initial guess (thin cat to tether)
X = 0.4; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.25 0.1 0 pi 100 2 0]);
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec4 = [0.523 0.525:0.005:0.65];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec4(j) = y(7,1);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) ==1
%          plot(z,r,'color',col1,'LineWidth',2);
    end
    if abs(X - 1.3/2) < 0.0001
%         plot(z,r,'color',[col1 6/7],'LineWidth',2);
%        shape1m = snapshot(0.7,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape1m = centercrop(shape1m,0.2,0.2); imwrite(shape1m,'shape1mt.png');
    end

   y_guess = y;
end

%% calculate shape for various widths (tether)
Xvecinset = [0.65:0.05:1.5];
muvecinset = zeros(1,length(Xvecinset));
Fvecinset = zeros(1,length(Xvecinset));
Evecinset = zeros(1,length(Xvecinset));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvecinset)
    X = Xvecinset(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvecinset(j) = y(7,1);
    Fvecinset(j) = -2*pi*y(9,1);
    Evecinset(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) ==1
%          plot(z,r,'color',col1,'LineWidth',2);
    end

     if abs(X - 1.35) < 0.001
%     plot(z,r,'color',[col1 1],'LineWidth',2);
%        shape1n = snapshot(1.4,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape1n = centercrop(shape1n,0.2,0.2); imwrite(shape1n,'shape1nt.tif','tif','resolution',600);
    end
    
   y_guess = y;
end

disp('n = 2')
%% prepare a solution as initial guess
X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether9); %n = 2, thick
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec5 = [0.45:0.005:0.59];
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
Evec5 = zeros(1,length(Xvec5));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec5)
    X = Xvec5(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec5(j) = y(7,1);
    Fvec5(j) = -2*pi*y(9,1);
    Evec5(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) ==1
%          plot(z,r,'color',col2,'LineWidth',2);
    end
    if abs(X - 0.59) < 0.001
%         plot(z,r,'color',col2,'LineWidth',2);
%         shape2K = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2K = centercrop(shape2K,0.2,0.2); imwrite(shape2K,'shape2K.png');
    end

   y_guess = y;
end

%% prepare a solution as initial guess
X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether9); %n = 2, thick
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec6 = [0.45:-0.005:0];
muvec6 = zeros(1,length(Xvec6));
Fvec6 = zeros(1,length(Xvec6));
Evec6 = zeros(1,length(Xvec6));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec6)
    X = Xvec6(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec6(j) = y(7,1);
    Fvec6(j) = -2*pi*y(9,1);
    Evec6(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col2,'LineWidth',2);
    end
    if abs(X - 0.7/2) < 0.001
        plot(z,r,'color',[col2 0.75],'LineWidth',2);
%         shape2J = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2J = centercrop(shape2J,0.2,0.2); imwrite(shape2J,'shape2J.png');
%         shape2Jp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Jp = centercrop(shape2Jp,0.2,0.2); imwrite(shape2Jp,'shape2Jp.png');
    elseif abs(X - 0.35/2) < 0.001
        plot(z,r,'color',[col2 0.5],'LineWidth',2);
%         shape2I = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2I = centercrop(shape2I,0.2,0.2); imwrite(shape2I,'shape2I.png');
%         shape2Ip = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Ip = centercrop(shape2Ip,0.2,0.2); imwrite(shape2Ip,'shape2Ip.png');
    elseif abs(X - 0) < 0.001
        plot(z,r,'color',[col2 0.25],'LineWidth',2);
%         shape2H = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2H = centercrop(shape2H,0.2,0.2); imwrite(shape2H,'shape2H.png');
%         shape2Hp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Hp = centercrop(shape2Hp,0.2,0.2); imwrite(shape2Hp,'shape2Hp.png');
    end

   y_guess = y;
end

return

%% prepare a solution as initial guess
X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess2_tether12); %n = 2, thin
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec7 = [0.45:-0.005:0.33];
muvec7 = zeros(1,length(Xvec7));
Fvec7 = zeros(1,length(Xvec7));
Evec7 = zeros(1,length(Xvec7));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec7)
    X = Xvec7(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec7(j) = y(7,1);
    Fvec7(j) = -2*pi*y(9,1);
    Evec7(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col2,'LineWidth',2);
    end
    if abs(X - 0.85/2) < 0.001
%         plot(z,r,'color',[col2 2/3],'LineWidth',2);
%         shape2j = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2j = centercrop(shape2j,0.2,0.2); imwrite(shape2j,'shape2jt.png');
%         shape2jp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2jp = centercrop(shape2jp,0.2,0.2); imwrite(shape2jp,'shape2jpt.png');
    elseif abs(X - 0.66/2) < 0.001
%         plot(z,r,'color',[col2 1/3],'LineWidth',2);
%         shape2i = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2i = centercrop(shape2i,0.2,0.2); imwrite(shape2i,'shape2it.png');
%         shape2ip = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2ip = centercrop(shape2ip,0.2,0.2); imwrite(shape2ip,'shape2ipt.png');
    end

   y_guess = y;
end

%% calculate shape for various widths
Xvec8 = [0.45:0.005:0.52]; %n = 2, thin
muvec8 = zeros(1,length(Xvec8));
Fvec8 = zeros(1,length(Xvec8));
Evec8 = zeros(1,length(Xvec8));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec8)
    X = Xvec8(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec8(j) = y(7,1);
    Fvec8(j) = -2*pi*y(9,1);
    Evec8(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col2,'LineWidth',2);
    end
    if abs(X - 0.52/2) < 0.001
%         plot(z,r,'color',col2,'LineWidth',2);
%         shape2k = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2k = centercrop(shape2k,0.2,0.2); imwrite(shape2k,'shape2kt.png');
%         shape2kp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2kp = centercrop(shape2kp,0.2,0.2); imwrite(shape2kp,'shape2kpt.png');
    end
    
   y_guess = y;
end

disp('n = 3')
%% prepare a solution as initial guess
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.21 0.1 0 pi -10 2 0]); %n = 3, upper
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec9 = [0.59:-0.005:0];
muvec9 = zeros(1,length(Xvec9));
Fvec9 = zeros(1,length(Xvec9));
Evec9 = zeros(1,length(Xvec9));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec9)
    X = Xvec9(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec9(j) = y(7,1);
    Fvec9(j) = -2*pi*y(9,1);
    Evec9(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col3,'LineWidth',2);
    end

    if abs(X - 0.59) < 0.001
       plot(z,r,'color',col3,'LineWidth',2);
%        shape3K = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape3K = centercrop(shape3K,0.2,0.2); imwrite(shape3K,'shape3K.png');
    elseif abs(X - 0.7/2) < 0.001
        plot(z,r,'color',[col3 0.75],'LineWidth',2);
%        shape3L = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape3L = centercrop(shape3L,0.2,0.2); imwrite(shape3L,'shape3L.png');
    elseif abs(X - 0.35/2) < 0.001
        plot(z,r,'color',[col3 0.5],'LineWidth',2);
%        shape3M = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape3M = centercrop(shape3M,0.2,0.2); imwrite(shape3M,'shape3M.png');
    elseif abs(X - 0) < 0.001
        plot(z,r,'color',[col3 0.25],'LineWidth',2);
%        shape3N = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape3N = centercrop(shape3N,0.2,0.2); imwrite(shape3N,'shape3N.png'); 
    end
   y_guess = y;
end

return
%% prepare a solution as initial guess
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.215 0.1 0 pi -10 2 0]);%n = 3, upper
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec10 = [0.52:-0.005:0.395];
muvec10 = zeros(1,length(Xvec10));
Fvec10 = zeros(1,length(Xvec10));
Evec10 = zeros(1,length(Xvec10));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec10)
    X = Xvec10(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec10(j) = y(7,1);
    Fvec10(j) = -2*pi*y(9,1);
    Evec10(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col3,'LineWidth',2);
    end

    if abs(X - 0.92/2) < 0.0001
%         plot(z,r,'color',[col3 2/3],'LineWidth',2);
%         shape3l = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3l = centercrop(shape3l,0.2,0.2); imwrite(shape3l,'shape3lt.png'); 
    elseif abs(X - 0.395) < 0.0001
%         plot(z,r,'color',[col3 1/3],'LineWidth',2);
%         shape3m = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3m = centercrop(shape3m,0.2,0.2); imwrite(shape3m,'shape3mt.png');        
    end
    
   y_guess = y;
end


%% prepare a solution as initial guess
X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.13 0.1 0 pi -10 2 0]); %n= 3, lower
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec11 = [0.59:-0.005:0];
muvec11 = zeros(1,length(Xvec11));
Fvec11 = zeros(1,length(Xvec11));
Evec11 = zeros(1,length(Xvec11));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec11)
    X = Xvec11(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec11(j) = y(7,1);
    Fvec11(j) = -2*pi*y(9,1);
    Evec11(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col3,'LineWidth',2);
    end
    
    if abs(X - 0.7/2) < 0.001
%         plot(z,r,'color',[col3 3/4],'LineWidth',2);
%        shape3J = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape3J = centercrop(shape3J,0.2,0.2); imwrite(shape3J,'shape3J.png');
    elseif abs(X - 0.35/2) < 0.001
%         plot(z,r,'color',[col3 2/4],'LineWidth',2);
%        shape3I = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape3I = centercrop(shape3I,0.2,0.2); imwrite(shape3I,'shape3I.png');
    elseif abs(X - 0) < 0.001
%         plot(z,r,'color',[col3 1/4],'LineWidth',2);
%        shape3H = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape3H = centercrop(shape3H,0.2,0.2); imwrite(shape3H,'shape3H.png'); 
    end

   y_guess = y;
end


%% prepare a solution as initial guess
X = 0.45; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.113 0.1 0 pi -10 2 0]); %n = 3
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec12 = [0.52:-0.005:0];
muvec12 = zeros(1,length(Xvec12));
Fvec12 = zeros(1,length(Xvec12));
Evec12 = zeros(1,length(Xvec12));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec12)
    X = Xvec12(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec12(j) = y(7,1);
    Fvec12(j) = -2*pi*y(9,1);
    Evec12(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col3,'LineWidth',2);
    end
    
    if abs(X - 0.52) < 0.0001
%         plot(z,r,'color',col3,'LineWidth',2);
%         shape3k = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3k = centercrop(shape3k,0.2,0.2); imwrite(shape3k,'shape3kt.png'); 
    elseif abs(X - 0.7/2) < 0.0001
%         plot(z,r,'color',[col3 3/4],'LineWidth',2);
%         shape3j = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3j = centercrop(shape3j,0.2,0.2); imwrite(shape3j,'shape3jt.png'); 
    elseif abs(X - 0.35/2) < 0.0001
%         plot(z,r,'color',[col3 2/4],'LineWidth',2);
%         shape3i = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3i = centercrop(shape3i,0.2,0.2); imwrite(shape3i,'shape3it.png'); 
    elseif abs(X - 0) < 0.0001   
%         plot(z,r,'color',[col3 1/4],'LineWidth',2);
%         shape3h = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3h = centercrop(shape3h,0.2,0.2); imwrite(shape3h,'shape3ht.png');         
    end

   y_guess = y;
end

disp('n = 4')
%% prepare a solution as initial guess
X = 0.45; solinit = bvpinit(linspace(0,1,6),@guess4_tether9); %n =4
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec13 = [0.59:-0.005:0];
muvec13 = zeros(1,length(Xvec13));
Fvec13 = zeros(1,length(Xvec13));
Evec13 = zeros(1,length(Xvec13));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec13)
    X = Xvec13(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec13(j) = y(7,1);
    Fvec13(j) = -2*pi*y(9,1);
    Evec13(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col4,'LineWidth',2);
    end
    
    if abs(X - 0.59) < 0.001
%         plot(z,r,'color',col4,'LineWidth',2);
%         shape4K = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4K = centercrop(shape4K,0.2,0.2); imwrite(shape4K,'shape4K.png');
    elseif abs(X - 0.7/2) < 0.001
%         plot(z,r,'color',[col4 3/4],'LineWidth',2);
%         shape4J = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4J = centercrop(shape4J,0.2,0.2); imwrite(shape4J,'shape4J.png');
%         shape4Jp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Jp = centercrop(shape4Jp,0.2,0.2); imwrite(shape4Jp,'shape4Jp.png');
    elseif abs(X - 0.35/2) < 0.001
%         plot(z,r,'color',[col4 2/4],'LineWidth',2);
%         shape4I = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4I = centercrop(shape4I,0.2,0.2); imwrite(shape4I,'shape4I.png');
%         shape4Ip = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Ip = centercrop(shape4Ip,0.2,0.2); imwrite(shape4Ip,'shape4Ip.png');
    elseif abs(X - 0) < 0.001
%         plot(z,r,'color',[col4 1/4],'LineWidth',2);
%         shape4H = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4H = centercrop(shape4H,0.2,0.2); imwrite(shape4H,'shape4H.png');
%         shape4Hp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Hp = centercrop(shape4Hp,0.2,0.2); imwrite(shape4Hp,'shape4Hp.png');        
    end

   y_guess = y;
end

%% prepare a solution as initial guess
N = 1001; X = 0.5; solinit = bvpinit(linspace(0,1,6),@guess4_tether10); %n =4
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec14= [0.52:-0.005:0.295];
muvec14 = zeros(1,length(Xvec14));
Fvec14 = zeros(1,length(Xvec14));
Evec14 = zeros(1,length(Xvec14));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec14)
    X = Xvec14(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec14(j) = y(7,1);
    Fvec14(j) = -2*pi*y(9,1);
    Evec14(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col4,'LineWidth',2);
    end
    if abs(X - 0.52) < 0.001
%         plot(z,r,'color',col4,'LineWidth',2);
%         shape4k = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4k = centercrop(shape4k,0.2,0.2); imwrite(shape4k,'shape4kt.png');
%         shape4kp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4kp = centercrop(shape4kp,0.2,0.2); imwrite(shape4kp,'shape4kpt.png');
    elseif abs(X - 0.91/2) < 0.001
%         plot(z,r,'color',[col4 2/3],'LineWidth',2);
%         shape4j = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4j = centercrop(shape4j,0.2,0.2); imwrite(shape4j,'shape4jt.png');
%         shape4jp = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4jp = centercrop(shape4jp,0.2,0.2); imwrite(shape4jp,'shape4jpt.png');
    elseif abs(X - 0.295) < 0.001
%         plot(z,r,'color',[col4 1/3],'LineWidth',2);
%         shape4i = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4i = centercrop(shape4i,0.2,0.2); imwrite(shape4i,'shape4it.png');
%         shape4ip = snapshot(0.6,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4ip = centercrop(shape4ip,0.2,0.2); imwrite(shape4ip,'shape4ipt.png');
    end
   y_guess = y;
end

disp('n = 5')
%% prepare a solution as initial guess
N = 1001; X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n= 5, lower
%  X = 0.466; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.514; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.352 0.1 0 pi -10 2 0]); %n = 5, upper
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.142 0.1 0 pi -150 2 70]); %n = 5 lower
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec15= [0.59:-0.005:0];
muvec15 = zeros(1,length(Xvec15));
Fvec15 = zeros(1,length(Xvec15));
Evec15 = zeros(1,length(Xvec15));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec15)
    X = Xvec15(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec15(j) = y(7,1);
    Fvec15(j) = -2*pi*y(9,1);
    Evec15(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col5,'LineWidth',2);
    end
    if abs(X - 0.59) < 0.001
%         plot(z,r,'color',col5,'LineWidth',2);
%        shape5K = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5K = centercrop(shape5K,0.2,0.2); imwrite(shape5K,'shape5K.png');
    elseif abs(X - 0.7/2) < 0.001
%         plot(z,r,'color',[col5 3/4],'LineWidth',2);
%        shape5L = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5L = centercrop(shape5L,0.2,0.2); imwrite(shape5L,'shape5L.png');
    elseif abs(X - 0.35/2) < 0.001
%         plot(z,r,'color',[col5 2/4],'LineWidth',2);
%        shape5M = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5M = centercrop(shape5M,0.2,0.2); imwrite(shape5M,'shape5M.png');
    elseif abs(X - 0) < 0.001
%         plot(z,r,'color',[col5 1/4],'LineWidth',2);
%        shape5N = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5N = centercrop(shape5N,0.2,0.2); imwrite(shape5N,'shape5N.png'); 
    end
   y_guess = y;
end

%% prepare a solution as initial guess
N = 1001; X = 0.466; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.514; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.44 0.1 0 pi -10 2 0]); %n = 5, upper
%  X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.352 0.1 0 pi -10 2 0]); %n = 5, upper
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.142 0.1 0 pi -150 2 70]); %n = 5 lower
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec16= [0.52:-0.005:0];
muvec16 = zeros(1,length(Xvec16));
Fvec16 = zeros(1,length(Xvec16));
Evec16 = zeros(1,length(Xvec16));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec16)
    X = Xvec16(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec16(j) = y(7,1);
    Fvec16(j) = -2*pi*y(9,1);
    Evec16(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col5,'LineWidth',2);
    end
    if abs(X - 0.52) < 0.001
%         plot(z,r,'color',col5,'LineWidth',2);
%        shape5k = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5k = centercrop(shape5k,0.2,0.2); imwrite(shape5k,'shape5kt.png');
    elseif abs(X - 0.7/2) < 0.001
%         plot(z,r,'color',[col5 3/4],'LineWidth',2);
%        shape5j = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5j = centercrop(shape5j,0.2,0.2); imwrite(shape5j,'shape5jt.png');
    elseif abs(X - 0.35/2) < 0.001
%         plot(z,r,'color',[col5 2/4],'LineWidth',2);
%        shape5i = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5i = centercrop(shape5i,0.2,0.2); imwrite(shape5i,'shape5it.png');
    elseif abs(X - 0) < 0.001
%         plot(z,r,'color',[col5 1/4],'LineWidth',2);
%        shape5h = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5h = centercrop(shape5h,0.2,0.2); imwrite(shape5h,'shape5ht.png'); 
    end
    
   y_guess = y;
   
end

%% prepare a solution as initial guess
N = 1001; X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.352 0.1 0 pi -10 2 0]); %n = 5, upper
% X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.142 0.1 0 pi -150 2 70]); %n = 5 lower
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec17= [0.52:-0.005:0.73/2];
muvec17 = zeros(1,length(Xvec17));
Fvec17 = zeros(1,length(Xvec17));
Evec17 = zeros(1,length(Xvec17));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec17)
    X = Xvec17(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec17(j) = y(7,1);
    Fvec17(j) = -2*pi*y(9,1);
    Evec17(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col5,'LineWidth',2);
    end
    if abs(X - 0.89/2) < 0.001
%         plot(z,r,'color',[col5 2/3],'LineWidth',2);
%        shape5l = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5l = centercrop(shape5l,0.2,0.2); imwrite(shape5l,'shape5lt.png');
    elseif abs(X - 0.73/2) < 0.001
%         plot(z,r,'color',[col5 1/3],'LineWidth',2);
%        shape5m = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5m = centercrop(shape5m,0.2,0.2); imwrite(shape5m,'shape5mt.png');
    end
   y_guess = y;
end


%% prepare a solution as initial guess
N = 1001; %X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.352 0.1 0 pi -10 2 0]); %n = 5, upper
X = 0.47; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.142 0.1 0 pi -150 2 70]); %n = 5 lower
sol = bvp4c(@odesystem,@bcs,solinit); t = linspace(0,1,N); y_guess = deval(sol,t);

%% calculate shape for various widths
Xvec18= [0.59:-0.005:0];
muvec18 = zeros(1,length(Xvec18));
Fvec18 = zeros(1,length(Xvec18));
Evec18 = zeros(1,length(Xvec18));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec18)
    X = Xvec18(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:); L = y(8,1);
    muvec18(j) = y(7,1);
    Fvec18(j) = -2*pi*y(9,1);
    Evec18(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    if mod(j,10) == 1
%          plot(z,r,'color',col5,'LineWidth',2);
    end
    if abs(X - 0.7/2) < 0.001
        plot(z,r,'color',[col5 3/4],'LineWidth',2);
%        shape5J = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5J = centercrop(shape5J,0.2,0.2); imwrite(shape5J,'shape5J.png');
    elseif abs(X - 0.35/2) < 0.001
        plot(z,r,'color',[col5 2/4],'LineWidth',2);
%        shape5I = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5I = centercrop(shape5I,0.2,0.2); imwrite(shape5I,'shape5I.png');
    elseif abs(X - 0) < 0.001
        plot(z,r,'color',[col5 1/4],'LineWidth',2);
%        shape5H = snapshot(0.6,1.35,psi,psi_s,z,r,X,kappa,N);
%        shape5H = centercrop(shape5H,0.2,0.2); imwrite(shape5H,'shape5H.png'); 
    end
    
   y_guess = y;
end



%% ==== catenoids ====
v1 = linspace(-0.590668192994804,0.590668192994804,101);
v2 = linspace(-0.521378393180893,0.521378393180893,101);
b1 = 0.757966837819559;
b2 = 0.255417422349254;
shape1K = snapshot(0.6,1.35,zeros(size(v1)),zeros(size(v1)),v1,b1*cosh(v1/b1),v1(end),kappa,101);
shape1K = centercrop(shape1K,0.2,0.2); imwrite(shape1K,'shape1K.png');
shape1l = snapshot(0.6,1.35,zeros(size(v2)),zeros(size(v2)),v2,b2*cosh(v2/b2),v2(end),kappa,101);
shape1l = centercrop(shape1l,0.2,0.2); imwrite(shape1l,'shape1lt.png');
% plot(v1,b1*cosh(v1/b1),'k--','linewidth',3); hold on;
% plot(v2,b2*cosh(v2/b2),'k--','linewidth',3);

% bcat = 0.825517; v = linspace(-0.527697,0.527697,101);
% r = bcat*cosh(v/bcat);
% shape1D = snapshot(0.6,1.35,zeros(size(v)),zeros(size(v)),v,r,v(end),kappa,N);
% shape1D = centercrop(shape1D,0.2,0.2); imwrite(shape1D,'shape1D.png');  

set(gca,'FontSize',18);
set(gcf,'color','w');
close all;

figure();
lw = 3;
% subplot(2,2,[1 3])
subplot('position',[0.1 0.1 0.25 0.85])
plot(2*Xvec,Fvec,'-','linewidth',lw,'color',col1); hold on;
% plot(2*Xvec2,Fvec2,'-','linewidth',lw,'color',col1); hold on;
Xvec2s = Xvec2(Xvec2>0.09); Xvec2u = Xvec2(Xvec2<=0.09);
Fvec2s = Fvec2(Xvec2>0.09); Fvec2u = Fvec2(Xvec2<=0.09);
plot(2*Xvec2s,Fvec2s,'-','linewidth',lw,'color',col1);
plot(2*Xvec2u,Fvec2u,':','linewidth',lw,'color',col1);
plot(2*Xvec3,Fvec3,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec4,Fvec4,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec5,Fvec5,':','linewidth',lw,'color',col2); hold on;
plot(2*Xvec6,Fvec6,':','linewidth',lw,'color',col2); hold on;
plot(2*Xvec7,Fvec7,'-','linewidth',lw,'color',col2); hold on;
plot(2*Xvec8,Fvec8,'-','linewidth',lw,'color',col2); hold on;
plot(2*Xvec9,Fvec9,':','linewidth',lw,'color',col3); hold on;
% plot(2*Xvec10,Fvec10,'-','linewidth',lw,'color',col3); hold on;
Xvec10s = Xvec10(Xvec10<0.5) 
Xvec10u = [Xvec10(Xvec10>=0.5) Xvec10s(1)]
Fvec10s = Fvec10(Xvec10<0.5) 
Fvec10u = [Fvec10(Xvec10>=0.5) Fvec10s(1)]
plot(2*Xvec10s,Fvec10s,'-','linewidth',lw,'color',col3);
plot(2*Xvec10u,Fvec10u,':','linewidth',lw,'color',col3);
plot(2*Xvec11,Fvec11,':','linewidth',lw,'color',col3); hold on;
% plot(2*Xvec12,Fvec12,'-','linewidth',lw,'color',col3); hold on;
Xvec12s = Xvec12(Xvec12<0.02); Xvec12u = Xvec12(Xvec12>=0.02);
Fvec12s = Fvec12(Xvec12<0.02); Fvec12u = Fvec12(Xvec12>=0.02);
plot(2*Xvec12s,Fvec12s,'-','linewidth',lw,'color',col3);
plot(2*Xvec12u,Fvec12u,':','linewidth',lw,'color',col3);
plot(2*Xvec13,Fvec13,':','linewidth',lw,'color',col4); hold on;
plot(2*Xvec14,Fvec14,':','linewidth',lw,'color',col4); hold on;
plot(2*Xvec15,Fvec15,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec16,Fvec16,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec17,Fvec17,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec18,Fvec18,':','linewidth',lw,'color',col5); hold on;
xlabel('$h/a$','interpreter','latex'); ylabel('$Fa/\kappa$','interpreter','latex');  grid on; box on;
set(gca,'FontSize',18);
% plot([1.18 1.18],[-2700, 100],'--','color',[0.6 0.6 0.6],'linewidth',2);
% plot([1.04 1.04],[-2700, 100],'--','color',[0.6 0.6 0.6],'linewidth',2);
% hh = text(1.15,-1900,'$h=h^*$','interpreter','latex','color',[0.6 0.6 0.6],'fontsize',18); set(hh,'rotation',90);
% vv = text(1.01,-1900,'$h=h^*$','interpreter','latex','color',[0.6 0.6 0.6],'fontsize',18); set(vv,'rotation',90);
xlim([0 0.65*2])
ylim([-2700 100]) 
% thick catenoid branches
plotpoint(0,-160.1,'1H',col1,-pi/2,0.08,50);
plotpoint(0.35,-141.5,'1I',col1,-pi/2,0.08,50);
plotpoint(0.70,-106.6,'1J',col1,-pi/2,0.08,50);
plotpoint(1.18,-18,'1K',col1,-pi/2,0.08,50);
plotpoint(1.045,15,'1\lambda',col1,pi/2,0.08,50);
plotpoint(1.3,71.25,'1m',col1,-pi/2,0.08,50);
% plotpoint(0,-189.5,'1G',col1,-pi/2,0.08,50);
plotpoint(0,-321,'2H',col2,-pi/2,0.08,50);
plotpoint(0.35,-322.9,'2I',col2,-pi/2,0.08,50);
plotpoint(0.7,-259.7,'2J',col2,-pi/2,0.08,50);
plotpoint(1.18,-102.2,'2K',col2,-pi/2,0.08,50);
plotpoint(0,-932.6,'3H',col3,pi/2,0.08,50);
plotpoint(0.35,-708.1,'3I',col3,pi/2,0.08,50);
plotpoint(0.7,-528,'3J',col3,pi/2,0.08,50);
plotpoint(1.18,-240,'3K',col3,-pi/2,0.08,50);
plotpoint(0.7,-621.4,'3L',col3,-pi/2,0.08,50);
plotpoint(0.35,-790.8,'3M',col3,-pi/2,0.08,50);
plotpoint(0,-935,'3N',col3,-pi/2,0.08,50);
plotpoint(0,-1612,'4H',col4,-pi/3,0.08,50);
plotpoint(0.35,-1321,'4I',col4,-pi/2.5,0.08,50);
plotpoint(0.7,-1018,'4J',col4,-pi/2.5,0.08,50);
plotpoint(1.18,-435.5,'4K',col4,-pi/2,0.08,50);
plotpoint(0,-2581,'5H',col5,pi/2,0.08,50);
plotpoint(0.35,-2023,'5I',col5,3*pi/4,0.08,50);
plotpoint(0.7,-1539,'5J',col5,3*pi/4,0.08,50);
plotpoint(1.18,-684.5,'5K',col5,pi/2,0.08,50);
plotpoint(0.7,-1633,'5L',col5,-pi/4,0.08,50);
plotpoint(0.35,-2105,'5M',col5,-pi/4,0.08,50);
plotpoint(0,-2611,'5N',col5,0,0.08,50);

% zoomed in portion
plot([0 1.04 1.04 0 0],[-95 -95 17 17 -95],'k--','linewidth',2);

% inset
axes('position',[.22 .13 .12 .2]); box on; 
plot([2*Xvec4 2*Xvecinset],[Fvec4 Fvecinset],'-','linewidth',lw,'color',col1); hold on; grid on;
indexOfInterest = (Fvec3 >=0) & (Xvec3 < 3);
plot(2*Xvec3(indexOfInterest),Fvec3(indexOfInterest),'-','linewidth',lw,'color',col1);
plotpoint(1.045,15,'1\lambda',col1,0,0.12,15);
plotpoint(1.3,71.25,'1m',col1,-pi/2,0.08,15);
plotpoint(2.7,164.8,'1n',col1,-pi/2,0.08,15);


% thin branches
% subplot(2,2,2)
subplot('position',[0.5 0.3 0.25 0.3])
plot(2*Xvec7,Fvec7,'-','linewidth',lw,'color',col2); hold on;
plot(2*Xvec8,Fvec8,'-','linewidth',lw,'color',col2); hold on;
% plot(2*Xvec10,Fvec10,'-','linewidth',lw,'color',col3); hold on;
plot(2*Xvec10s,Fvec10s,'-','linewidth',lw,'color',col3);
plot(2*Xvec10u,Fvec10u,':','linewidth',lw,'color',col3);
% plot(2*Xvec12,Fvec12,'-','linewidth',lw,'color',col3); hold on;
plot(2*Xvec12s,Fvec12s,'-','linewidth',lw,'color',col3);
plot(2*Xvec12u,Fvec12u,':','linewidth',lw,'color',col3);
plot(2*Xvec14,Fvec14,':','linewidth',lw,'color',col4); hold on;
plot(2*Xvec16,Fvec16,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec17,Fvec17,':','linewidth',lw,'color',col5); hold on;
xlabel('$h/a$','interpreter','latex'); ylabel('$F a/\kappa$','interpreter','latex');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');
% plot([1.04 1.04],[-100, 20],'--','color',[0.6 0.6 0.6],'linewidth',2);
% gg = text(1.015,-3.5,'$h=h^*$','interpreter','latex','color',[0.6 0.6 0.6],'fontsize',18); set(gg,'rotation',90);
xlim([0 1.05])
ylim([-100 15])
plotpoint(0.66,-3.34,'2i',col2,pi/2,0.08,5);
plotpoint(0.85,-5.84,'2j',col2,pi/2,0.08,5);
plotpoint(1.04,-11.57,'2k',col2,pi/2,0.08,5);
plotpoint(0,4.69,'3h',col3,-pi/2,0.08,5);
plotpoint(0.35,0.614,'3i',col3,pi/2,0.08,5);
plotpoint(0.7,-8.82,'3j',col3,-pi/2,0.08,5);
plotpoint(1.04,-29.7,'3k',col3,-pi/2,0.08,5);
plotpoint(0.92,-10.12,'3\lambda',col3,-pi/2,0.08,5);
plotpoint(0.79,-4.67,'3m',col3,-pi/2,0.08,5);
plotpoint(0.59,-10.57,'4i',col4,pi/2,0.08,5);
plotpoint(0.91,-34.41,'4j',col4,pi/2,0.08,5);
plotpoint(1.04,-59.5,'4k',col4,pi/2,0.08,5);
plotpoint(0,12.07,'5h',col5,pi/2,0.08,5);
plotpoint(0.35,-4.68,'5i',col5,-pi/2,0.08,5);
plotpoint(0.7,-30.96,'5j',col5,-pi/2,0.08,5);
plotpoint(1.04,-97,'5k',col5,pi,0.04,5);
plotpoint(0.89,-38.2,'5\lambda',col5,-pi/2,0.08,5);
plotpoint(0.73,-17.66,'5m',col5,-pi/2,0.08,5);

% return; 
% figure()
% subplot('position',[0.1 0.1 0.25 0.85])
% plot(2*Xvec,Evec,'-','linewidth',lw,'color',col1); hold on;
% plot(2*Xvec2,Evec2,'-','linewidth',lw,'color',col1); hold on;
% plot(2*Xvec3,Evec3,'-','linewidth',lw,'color',col1); hold on;
% plot(2*Xvec4,Evec4,'-','linewidth',lw,'color',col1); hold on;
% plot(2*Xvec5,Evec5,'-','linewidth',lw,'color',col2); hold on;
% plot(2*Xvec6,Evec6,'-','linewidth',lw,'color',col2); hold on;
% plot(2*Xvec7,Evec7,'-','linewidth',lw,'color',col2); hold on;
% plot(2*Xvec8,Evec8,'-','linewidth',lw,'color',col2); hold on;
% plot(2*Xvec9,Evec9,'-','linewidth',lw,'color',col3); hold on;
% plot(2*Xvec10,Evec10,'-','linewidth',lw,'color',col3); hold on;
% plot(2*Xvec11,Evec11,'-','linewidth',lw,'color',col3); hold on;
% plot(2*Xvec12,Evec12,'-','linewidth',lw,'color',col3); hold on;
% plot(2*Xvec13,Evec13,'-','linewidth',lw,'color',col4); hold on;
% plot(2*Xvec14,Evec14,'-','linewidth',lw,'color',col4); hold on;
% plot(2*Xvec15,Evec15,'-','linewidth',lw,'color',col5); hold on;
% plot(2*Xvec16,Evec16,'-','linewidth',lw,'color',col5); hold on;
% plot(2*Xvec17,Evec17,'-','linewidth',lw,'color',col5); hold on;
% plot(2*Xvec18,Evec18,'-','linewidth',lw,'color',col5); hold on;
% xlabel('$h/a$','interpreter','latex'); ylabel('$E/\kappa$','interpreter','latex');  grid on; box on;
% set(gca,'FontSize',18);
% set(gcf,'color','w');

figure();
subplot('position',[0.1 0.1 0.25 0.85])
plot(2*Xvec,muvec,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec2,muvec2,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec3,muvec3,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec4,muvec4,'-','linewidth',lw,'color',col1); hold on;
plot(2*Xvec5,muvec5,'-','linewidth',lw,'color',col2); hold on;
plot(2*Xvec6,muvec6,'-','linewidth',lw,'color',col2); hold on;
plot(2*Xvec7,muvec7,'-','linewidth',lw,'color',col2); hold on;
plot(2*Xvec8,muvec8,'-','linewidth',lw,'color',col2); hold on;
plot(2*Xvec9,muvec9,'-','linewidth',lw,'color',col3); hold on;
plot(2*Xvec10,muvec10,'-','linewidth',lw,'color',col3); hold on;
plot(2*Xvec11,muvec11,'-','linewidth',lw,'color',col3); hold on;
plot(2*Xvec12,muvec12,'-','linewidth',lw,'color',col3); hold on;
plot(2*Xvec13,muvec13,'-','linewidth',lw,'color',col4); hold on;
plot(2*Xvec14,muvec14,'-','linewidth',lw,'color',col4); hold on;
plot(2*Xvec15,muvec15,'-','linewidth',lw,'color',col5); hold on;
plot(2*Xvec16,muvec16,'-','linewidth',lw,'color',col5); hold on;
plot(2*Xvec17,muvec17,'-','linewidth',lw,'color',col5); hold on;
plot(2*Xvec18,muvec18,'-','linewidth',lw,'color',col5); hold on;
xlabel('$h/a$','interpreter','latex'); ylabel('$\mu a^2/\kappa$','interpreter','latex');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');


end

%% makes figures for A = 2*pi*1.3
function [] = nocatfig()

global kappa; kappa = 1;
global kbar; kbar = 0;
global A; A = 2*pi*1.3; %surface area
global X; X = 1.04/2; %half width; %0.24/2 for upper branches, 0.2/2 for thin ones
global y_guess;
global N; N = 1001; %number of gridpoints

%% ==== colors ====
col1 = [0, 0.4470, 0.7410];
col2 = [0.8500, 0.3250, 0.0980];
col3 = [0.9290, 0.6940, 0.1250];
col4 = [0.4940, 0.1840, 0.5560];
col5 = [0.4660, 0.6740, 0.1880];

disp('n = 1')
%% ------------------ 1st branch --------------------
%% prepare a solution as initial guess
t = linspace(0,1,N);
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.5 0.1 0 pi 0 1.7 0]);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
y_guess = deval(sol,t);
Xvec2 = [ 0:0.005:1.5];
muvec2 = zeros(1,length(Xvec2));
Fvec2 = zeros(1,length(Xvec2));
Evec2 = zeros(1,length(Xvec2));
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec2)
    X = Xvec2(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec2(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec2(j) = -2*pi*y(9,1);
    Evec2(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%         plot(z,r,'color',col1,'LineWidth',2);
    end
    if abs(2*X - 0) < 0.00001
        plot(z,r,'color',[col1 1/5],'LineWidth',2);
%         shape1O = snapshot(0.6,1.5,psi,psi_s,z,r,X,kappa,N);
%         shape1O = centercrop(shape1O,0.2,0.2); imwrite(shape1O,'shape1O.png');
    elseif abs(2*X - 0.5) < 0.00001
        plot(z,r,'color',[col1 2/5],'LineWidth',2);
%         shape1P = snapshot(0.6,1.5,psi,psi_s,z,r,X,kappa,N);
%         shape1P = centercrop(shape1P,0.2,0.2); imwrite(shape1P,'shape1P.png');
    elseif abs(2*X - 1) < 0.0001
        plot(z,r,'color',[col1 3/5],'LineWidth',2);
%         shape1Q = snapshot(0.6,1.5,psi,psi_s,z,r,X,kappa,N);
%         shape1Q = centercrop(shape1Q,0.2,0.2); imwrite(shape1Q,'shape1Q.png');
    elseif abs(2*X - 1.5) < 0.0001
        plot(z,r,'color',[col1 4/5],'LineWidth',2);
%         shape1R = snapshot(0.8,1.5,psi,psi_s,z,r,X,kappa,N);
%         shape1R = centercrop(shape1R,0.2,0.2); imwrite(shape1R,'shape1R.png');
    elseif abs(2*X - 2.5) < 0.0001
        plot(z,r,'color',col1,'LineWidth',2);
%         shape1s = snapshot(1.6,1.5,psi,psi_s,z,r,X,kappa,N);
%         shape1s = centercrop(shape1s,0.2,0.2); imwrite(shape1s,'shape1st.png'); 
    end

    y_guess = y;
end

return
disp('n = 2')
%% ------------------ 1st branch --------------------
X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
y_guess = deval(sol,t);
Xvec3 = [0.6:0.01:0.7 0.701:0.001:0.712];
muvec3 = zeros(1,length(Xvec3));
Fvec3 = zeros(1,length(Xvec3));
Evec3 = zeros(1,length(Xvec3));
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec3)
    X = Xvec3(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec3(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec3(j) = -2*pi*y(9,1);
    Evec3(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%         plot(-z,r,'color',col2,'LineWidth',2);
    end
    if abs(X - 0.712) < 0.00001
%         plot(-z,r,'color',[col2 4/5],'linewidth',2);
%         shape2R = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2R = centercrop(shape2R,0.2,0.2); imwrite(shape2R,'shape2R.png');
%         shape2Rp = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Rp = centercrop(shape2Rp,0.2,0.2); imwrite(shape2Rp,'shape2Rp.png');
    end

    y_guess = y;
end

X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether6);
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
y_guess = deval(sol,t);
Xvec4 = [0.6:-0.01:0.47];
muvec4 = zeros(1,length(Xvec4));
Fvec4 = zeros(1,length(Xvec4));
Evec4 = zeros(1,length(Xvec4));
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec4)
    X = Xvec4(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec4(j) = y(7,1); mu = y(7,1);
    Fvec4(j) = -2*pi*y(9,1);
    Evec4(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%         plot(-z,r,'color',col2,'LineWidth',2);
    end
    if abs(X - 0.47) < 0.00001
%         plot(-z,r,'color',col2,'linewidth',2);
%         shape2s = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2s = centercrop(shape2s,0.2,0.2); imwrite(shape2s,'shape2st.png');
%         shape2sp = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2sp = centercrop(shape2sp,0.2,0.2); imwrite(shape2sp,'shape2spt.png');
    end

    y_guess = y;
end

%% ------------------ 2nd branch --------------------
X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5); %n=2
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
figure(1)
Xvec5 = [0.6:0.01:0.7 0.701:0.001:0.712];
muvec5 = zeros(1,length(Xvec5));
Fvec5 = zeros(1,length(Xvec5));
Evec5 = zeros(1,length(Xvec5));
for j = 1:length(Xvec5)
    X = Xvec5(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec5(j) = y(7,1); mu = y(7,1);
    Fvec5(j) = -2*pi*y(9,1);
    Evec5(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col2,'linewidth',2);
    end

    
    y_guess = y;
end

X = 0.6; solinit = bvpinit(linspace(0,1,5),@guess2_tether5); %n=2
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);
figure(1)
Xvec6 = [0.6:-0.01:0];
muvec6 = zeros(1,length(Xvec6));
Fvec6 = zeros(1,length(Xvec6));
Evec6 = zeros(1,length(Xvec6));
for j = 1:length(Xvec6)
    X = Xvec6(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec6(j) = y(7,1); mu = y(7,1);
    Fvec6(j) = -2*pi*y(9,1);
    Evec6(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col2,'linewidth',2);
    end
    if abs(X - 0.5) < 0.00001
%          plot(z,r,'color',[col2, 3/5],'linewidth',2);
%         shape2Q = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2Q = centercrop(shape2Q,0.2,0.2); imwrite(shape2Q,'shape2Q.png');
%         shape2Qp = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Qp = centercrop(shape2Qp,0.2,0.2); imwrite(shape2Qp,'shape2Qp.png');
    elseif abs(X - 0.25) < 0.00001
%         plot(z,r,'color',[col2 2/5],'linewidth',2);
%         shape2P = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2P = centercrop(shape2P,0.2,0.2); imwrite(shape2P,'shape2P.png');
%         shape2Pp = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Pp = centercrop(shape2Pp,0.2,0.2); imwrite(shape2Pp,'shape2Pp.png');
    elseif abs(X - 0) < 0.0001
%         plot(z,r,'color',[col2 1/5],'linewidth',2);
%         shape2O = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape2O = centercrop(shape2O,0.2,0.2); imwrite(shape2O,'shape2O.png');
%         shape2Op = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape2Op = centercrop(shape2Op,0.2,0.2); imwrite(shape2Op,'shape2Op.png');
    end
    y_guess = y;
end


disp('n = 3')
%% ------------------ 1st branch --------------------
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
y_guess = deval(sol,t);
Xvec7 = [1.04:0.01:1.4 1.401:0.001:1.431 1.4312]/2;
muvec7 = zeros(1,length(Xvec7));
Fvec7 = zeros(1,length(Xvec7));
Evec7 = zeros(1,length(Xvec7));
for j = 1:length(Xvec7)
    X = Xvec7(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec7(j) = y(7,1); mu = y(7,1);
    Fvec7(j) = -2*pi*y(9,1);
    Evec7(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%         plot(z,r,'color',col3,'LineWidth',2);
    end
    if abs(2*X - 1.4312) < 0.000001
%         plot(z,r,'color',[col3 4/5],'linewidth',2);
%         shape3R = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3R = centercrop(shape3R,0.2,0.2); imwrite(shape3R,'shape3R.png');
    end
    y_guess = y;
end

disp('8')
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.2 0.1 0 pi 0 1.7 0]); %n=3U
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
y_guess = deval(sol,t);
Xvec8 = [1.04:-0.01:0]/2;
muvec8 = zeros(1,length(Xvec8));
Fvec8 = zeros(1,length(Xvec8));
Evec8 = zeros(1,length(Xvec8));
for j = 1:length(Xvec8)
    X = Xvec8(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec8(j) = y(7,1); mu = y(7,1);
    Fvec8(j) = -2*pi*y(9,1);
    Evec8(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%         plot(z,r,'color',col3,'LineWidth',2);
    end
    if abs(2*X - 1) < 0.0001
        disp('3Q')
%         plot(z,r,'color',[col3 3/5],'linewidth',2);
%         shape3Q = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3Q = centercrop(shape3Q,0.2,0.2); imwrite(shape3Q,'shape3Q.png');
    elseif abs(2*X - 0.5) < 0.0001
%         plot(z,r,'color',[col3 2/5],'linewidth',2);
%         shape3P = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3P = centercrop(shape3P,0.2,0.2); imwrite(shape3P,'shape3P.png'); 
    elseif abs(2*X - 0) < 0.0001
%         plot(z,r,'color',[col3 1/5],'linewidth',2);
%         shape3O = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3O = centercrop(shape3O,0.2,0.2); imwrite(shape3O,'shape3O.png'); 
    end
    
    y_guess = y;
end

%% ------------------ 2nd branch --------------------
X = 1.24/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.1 0.11 0 pi -3 1.9 1.4]); %n=3T
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

disp('9')
Xvec9 = [1.04:0.01:1.41 1.411 1.412 1.413 1.414 1.4155]/2;
muvec9 = zeros(1,length(Xvec9));
Fvec9 = zeros(1,length(Xvec9));
Evec9 = zeros(1,length(Xvec9));
for j = 1:length(Xvec9)
    X = Xvec9(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec9(j) = y(7,1); mu = y(7,1);
    Fvec9(j) = -2*pi*y(9,1);
    Evec9(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col3,'linewidth',2);
    end
    if abs(2*X - 1.4155) < 0.00001
%         plot(z,r,'color',[col3 4/7],'linewidth',2);
%         shape3W = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3W = centercrop(shape3W,0.2,0.2); imwrite(shape3W,'shape3W.png');
    end

    y_guess = y;
end

disp('10')
X = 1.24/2; solinit = bvpinit(linspace(0,1,10),[pi/2 1 0.1 0.11 0 pi -3 1.9 1.4]); %n=3T
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec10 = [1.04:-0.01:0]/2;
muvec10 = zeros(1,length(Xvec10));
Fvec10 = zeros(1,length(Xvec10));
Evec10 = zeros(1,length(Xvec10));
for j = 1:length(Xvec10)
    X = Xvec10(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec10(j) = y(7,1); mu = y(7,1);
    Fvec10(j) = -2*pi*y(9,1);
    Evec10(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col3,'linewidth',2);
    end

    if abs(2*X - 1) < 0.0001
%         plot(z,r,'color',[col3 5/7],'linewidth',2);
%         shape3x = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3x = centercrop(shape3x,0.2,0.2); imwrite(shape3x,'shape3xt.png');
    elseif abs(2*X - 0.5) < 0.0001
%         plot(z,r,'color',[col3 6/7],'linewidth',2);
%         shape3y = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3y = centercrop(shape3y,0.2,0.2); imwrite(shape3y,'shape3yt.png');       
    elseif abs(2*X - 0) < 0.0001
%         plot(z,r,'color',col3,'linewidth',2);
%         shape3z = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3z = centercrop(shape3z,0.2,0.2); imwrite(shape3z,'shape3zt.png'); 
    end
    
    y_guess = y;
end
 
%% ------------------ 3rd branch --------------------
X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.14 0.1 0 pi 0 1.7 0]); %n=3L
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec11 = [1.04:0.01:1.4 1.401:0.001:1.415 1.4155]/2;
muvec11 = zeros(1,length(Xvec11));
Fvec11 = zeros(1,length(Xvec11));
Evec11 = zeros(1,length(Xvec11));
for j = 1:length(Xvec11)
    X = Xvec11(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec11(j) = y(7,1); mu = y(7,1);
    Fvec11(j) = -2*pi*y(9,1);
    Evec11(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col3,'linewidth',2);
    end
    y_guess = y;
end

X = 1.24/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.14 0.1 0 pi 0 1.7 0]); %n=3L
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec12 = [1.04:-0.01:0]/2;
muvec12 = zeros(1,length(Xvec12));
Fvec12 = zeros(1,length(Xvec12));
Evec12 = zeros(1,length(Xvec12));
for j = 1:length(Xvec12)
    X = Xvec12(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec12(j) = y(7,1); mu = y(7,1);
    Fvec12(j) = -2*pi*y(9,1);
    Evec12(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col3,'linewidth',2);
    end
    if abs(2*X - 1) < 0.0001
        disp('3V')
%         plot(z,r,'color',[col3 3/7],'linewidth',2);
%         shape3V = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3V = centercrop(shape3V,0.2,0.2); imwrite(shape3V,'shape3V.png'); 
    elseif abs(2*X - 0.5) < 0.0001
%         plot(z,r,'color',[col3 2/7],'linewidth',2);
%         shape3U = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3U = centercrop(shape3U,0.2,0.2); imwrite(shape3U,'shape3U.png'); 
    elseif abs(2*X - 0) < 0.0001
%         plot(z,r,'color',[col3 1/7],'linewidth',2);
%         shape3T = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3T = centercrop(shape3T,0.2,0.2); imwrite(shape3T,'shape3T.png'); 
    end
    y_guess = y;
end


%% ------------------ 4th branch --------------------
X = 1.43/2; solinit = bvpinit(linspace(0,1,10),[pi/2 2 0.16 0.11 0 pi -1 1.7 1]); %n=3S
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec13 = [1.4311 1.43:-0.001:1.401 1.4:-0.01:1.28]/2;
muvec13 = zeros(1,length(Xvec13));
Fvec13 = zeros(1,length(Xvec13));
Evec13 = zeros(1,length(Xvec13));
for j = 1:length(Xvec13)
    X = Xvec13(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec13(j) = y(7,1); mu = y(7,1);
    Fvec13(j) = -2*pi*y(9,1);
    Evec13(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col3,'linewidth',2);
    end
    if abs(2*X - 1.28) < 0.0001
%         plot(z,r,'color',col3,'linewidth',2);
%         shape3s = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape3s = centercrop(shape3s,0.2,0.2); imwrite(shape3s,'shape3st.png');
    end
    y_guess = y;
end
Xvec13 = [Xvec7(end) Xvec13];
Fvec13 = [Fvec7(end) Fvec13];
muvec13 = [muvec7(end) muvec13];



disp('n = 4')
%% ------------------ 1st branch --------------------
X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether6); %n =4
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    H = (psi_s + sin(psi)./r)/2;
    muvec14 = y(7,1);
    L = y(8,1);
    Fvec14 = -2*pi*y(9,1);
    indices = find([0 diff(sign(H))]~=0);
    plot(z,r,'b','linewidth',2);
    Xvec14 = 0.59;
    Evec14 = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));

options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
%Xvec = [0.58 0.56 0.54];
Xvec14 = [Xvec14 0.61 0.65 0.68 0.7 0.705 0.706 0.707 0.708 0.709];
% Xvec = [0.59:0.01:0.7 0.705];
% Xvec = [1.04:-0.01:0]/2;
muvec14 = [muvec14 zeros(1,length(Xvec14)-1)];
Fvec14 = [Fvec14 zeros(1,length(Xvec14)-1)];
Evec14 = [Evec14 zeros(1,length(Xvec14)-1)];
figure(1); hold on;  box on; grid on; axis equal;
for j = 2:length(Xvec14)
    X = Xvec14(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec14(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec14(j) = -2*pi*y(9,1);
    Evec14(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,1) == 0
%         plot(z,r,'color',col4,'LineWidth',2);
    end
    
    if abs(X - 0.709) < 0.0001
%         plot(z,r,'color',[col4 4/5],'LineWidth',2);
%         shape4R = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4R = centercrop(shape4R,0.2,0.2); imwrite(shape4R,'shape4R.png'); 
%         shape4Rp = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Rp = centercrop(shape4Rp,0.2,0.2); imwrite(shape4Rp,'shape4Rp.png');     
    end
    y_guess = y;
end

X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether6); %n =4
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec15 = [0.58 0.57 0.56];
muvec15 = [zeros(1,length(Xvec15))];
Fvec15 = [zeros(1,length(Xvec15))];
Evec15 = [ zeros(1,length(Xvec15))];
figure(1); hold on;  box on; grid on; axis equal;
for j = 1:length(Xvec15)
    X = Xvec15(j)
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec15(j) = y(7,1); mu = y(7,1);
    Fvec15(j) = -2*pi*y(9,1);
    Evec15(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,1) == 0
%         plot(z,r,'color',col4,'LineWidth',2);
    end
    if abs(X - 0.56) < 0.00001
%         plot(z,r,'color',col4,'LineWidth',2);
%         shape4s = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4s = centercrop(shape4s,0.2,0.2); imwrite(shape4s,'shape4st.png'); 
%         shape4sp = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4sp = centercrop(shape4sp,0.2,0.2); imwrite(shape4sp,'shape4spt.png'); 
    end
    y_guess = y;
end

Xvec15 = [Xvec14(1) Xvec15];
Fvec15 = [Fvec14(1) Fvec15];
muvec15 = [muvec14(1) muvec15];

%% ------------------ 2nd branch --------------------
disp('16')
X = 0.59; solinit = bvpinit(linspace(0,1,6),@guess4_tether5); %n=4
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec16 = [0.59:0.01:0.7 0.705 0.706 0.707 0.708 0.7085 0.7088];
muvec16 = zeros(1,length(Xvec16));
Fvec16 = zeros(1,length(Xvec16));
Evec16 = zeros(1,length(Xvec16));
for j = 1:length(Xvec16)
    X = Xvec16(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec16(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec16(j) = -2*pi*y(9,1);
    Evec16(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col4,'linewidth',2);
    end
    
    y_guess = y;
end
Xvec16 = [Xvec16 Xvec14(end)];
Fvec16 = [Fvec16 Fvec14(end)];
muvec16 = [muvec16 muvec14(end)];


Xvec17 = [0.59:-0.01:0];
muvec17 = zeros(1,length(Xvec17));
Fvec17 = zeros(1,length(Xvec17));
Evec17 = zeros(1,length(Xvec17));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec17)
    X = Xvec17(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec17(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec17(j) = -2*pi*y(9,1);
    Evec17(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col4,'linewidth',2);
    end
    
    if abs(X - 0.5) < 0.0001
%         plot(z,r,'color',[col4 3/5],'LineWidth',2);
%         shape4Q = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4Q = centercrop(shape4Q,0.2,0.2); imwrite(shape4Q,'shape4Q.png');  
%         shape4Qp = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Qp = centercrop(shape4Qp,0.2,0.2); imwrite(shape4Qp,'shape4Qp.png');
    elseif abs(X - 0.25) < 0.0001
%         plot(z,r,'color',[col4 2/5],'LineWidth',2);
%         shape4P = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4P = centercrop(shape4P,0.2,0.2); imwrite(shape4P,'shape4P.png');  
%         shape4Pp = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Pp = centercrop(shape4Pp,0.2,0.2); imwrite(shape4Pp,'shape4Pp.png');  
    elseif abs(X - 0) < 0.0001
%         plot(z,r,'color',[col4 1/5],'LineWidth',2);
%         shape4O = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape4O = centercrop(shape4O,0.2,0.2); imwrite(shape4O,'shape4O.png');  
%         shape4Op = snapshot(0.75,1.35,psi,psi_s,-z,r,X,kappa,N);
%         shape4Op = centercrop(shape4Op,0.2,0.2); imwrite(shape4Op,'shape4Op.png');  
    end
    y_guess = y;
end


disp('n = 5')
%% ------------------ 1st branch --------------------
disp('18')
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
Xvec18 = [1.05:0.01:1.40 1.401:0.001:1.415 1.4161]/2;
%  Xvec = [1.05:-0.01:0]/2;
muvec18 = zeros(1,length(Xvec18));
Fvec18 = zeros(1,length(Xvec18));
Evec18 = zeros(1,length(Xvec18));
for j = 1:length(Xvec18)
    X = Xvec18(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec18(j) = y(7,1); mu = y(7,1);
    Fvec18(j) = -2*pi*y(9,1);
    Evec18(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%         plot(z,r,'color',col5,'LineWidth',2);
    end
    if abs(X - 1.4161/2) < 0.00001
%         plot(z,r,'color',[col5 4/5],'LineWidth',2);
%         shape5R = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5R = centercrop(shape5R,0.2,0.2); imwrite(shape5R,'shape5R.png');  
    end

    y_guess = y;
end

disp('19')
X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.09 0.1 0 pi 0 1.7 0]); %n=5U
options = bvpset('RelTol',1e-3);
sol = bvp4c(@odesystem,@bcs,solinit,options);
t = linspace(0,1,N);
y_guess = deval(sol,t);
Xvec19 = [1.05:-0.01:0]/2;
muvec19 = zeros(1,length(Xvec19));
Fvec19 = zeros(1,length(Xvec19));
Evec19 = zeros(1,length(Xvec19));
for j = 1:length(Xvec19)
    X = Xvec19(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec19(j) = y(7,1); mu = y(7,1);
    Fvec19(j) = -2*pi*y(9,1);
    Evec19(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%         plot(z,r,'color',col5,'LineWidth',2);
    end
    if abs(2*X - 1) < 0.0001
%         plot(z,r,'color',[col5 3/5],'LineWidth',2);
%         shape5Q = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5Q = centercrop(shape5Q,0.2,0.2); imwrite(shape5Q,'shape5Q.png');  
    elseif abs(2*X - 0.5) < 0.0001
%         plot(z,r,'color',[col5 2/5],'LineWidth',2);
%         shape5P = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5P = centercrop(shape5P,0.2,0.2); imwrite(shape5P,'shape5P.png');
    elseif abs(2*X - 0) < 0.0001
%         plot(z,r,'color',[col5 1/5],'LineWidth',2);
%         shape5O = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5O = centercrop(shape5O,0.2,0.2); imwrite(shape5O,'shape5O.png');      
    end

    y_guess = y;
end

%% ------------------ 2nd branch --------------------
disp('20')
 X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec20 = [1.04:0.01:1.4 1.401:0.001:1.417 1.4172 1.4174 1.4176 1.4178]/2;
% Xvec2 = [1.04:-0.01:0]/2;
muvec20 = zeros(1,length(Xvec20));
Fvec20 = zeros(1,length(Xvec20));
Evec20 = zeros(1,length(Xvec20));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec20)
    X = Xvec20(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec20(j) = y(7,1); mu = y(7,1);
    Fvec20(j) = -2*pi*y(9,1);
    Evec20(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'r','linewidth',2);
    end
    if abs(X - 1.4178/2) < 0.000001
        plot(z,r,'color',[col5 4/5],'linewidth',2);
%         shape5W = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5W = centercrop(shape5W,0.2,0.2); imwrite(shape5W,'shape5W.png'); 
    end
    y_guess = y;
end

 X = 1.04/2; solinit = bvpinit(linspace(0,1,11),[pi/2 0.01 0.075 0.1 0 pi 0 1.7 0]); %n=5L
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec21 = [1.04:-0.01:0]/2;
muvec21 = zeros(1,length(Xvec21));
Fvec21 = zeros(1,length(Xvec21));
Evec21 = zeros(1,length(Xvec21));
for j = 1:length(Xvec21)
    X = Xvec21(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec21(j) = y(7,1); mu = y(7,1);
    Fvec21(j) = -2*pi*y(9,1);
    Evec21(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,5) == 1
%       plot(z,r,'color',col5,'linewidth',2);
    end
    if abs(X - 0.5) < 0.000001
        plot(z,r,'color',[col5 3/5],'linewidth',2);
%         shape5V = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5V = centercrop(shape5V,0.2,0.2); imwrite(shape5V,'shape5V.png'); 
    elseif abs(X - 0.25) < 0.000001
        plot(z,r,'color',[col5 2/5],'linewidth',2);
%         shape5U = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5U = centercrop(shape5U,0.2,0.2); imwrite(shape5U,'shape5U.png'); 
    elseif abs(X - 0) < 0.000001
        plot(z,r,'color',[col5 1/5],'linewidth',2);
%         shape5T = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5T = centercrop(shape5T,0.2,0.2); imwrite(shape5T,'shape5T.png');   
    end
    y_guess = y;
end
 
%% ------------------ 3rd branch --------------------
disp('22')
X = 1.4/2; solinit = bvpinit(linspace(0,1,11),[pi/2 2 0.13 0.1 0 pi -70 1.3 90]); %n=5S
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec22 = [1.4178 1.4174 1.417 1.415:-0.001:1.401 1.4:-0.01:1.3 1.295:-0.005:1.225 1.2]/2;
% Xvec3 = [1.04:-0.01:0]/2;
muvec22 = zeros(1,length(Xvec22));
Fvec22 = zeros(1,length(Xvec22));
Evec22 = zeros(1,length(Xvec22));
for j = 1:length(Xvec22)
    X = Xvec22(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec22(j) = y(7,1); mu = y(7,1);
    Fvec22(j) = -2*pi*y(9,1);
    Evec22(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,3) == 1
%       plot(z,r,'color',col5,'linewidth',2);
    end
    y_guess = y;
    if abs(X - 1.225/2) < 0.0001
        plot(z,r,'color',col5,'linewidth',2);
%         shape5x = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5x = centercrop(shape5x,0.2,0.2); imwrite(shape5x,'shape5xt.png'); 
    end
        
    
end
Xvec22 = [Xvec20(end) Xvec22];
Fvec22 = [Fvec20(end) Fvec22];
muvec22 = [muvec20(end) muvec22];


%% ------------------ 4th branch --------------------
disp('23')
X = 1.4/2; solinit = bvpinit(linspace(0,1,9),[pi/2 2 0.14 0.1 0 pi -70 1.3 90]); %n=5T
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

Xvec23 = [1.415:-0.001:1.4 1.395:-0.01:1.315 1.3149 1.3148 ]/2;
muvec23 = zeros(1,length(Xvec23));
Fvec23 = zeros(1,length(Xvec23));
Evec23 = zeros(1,length(Xvec23));
for j = 1:length(Xvec23)
    X = Xvec23(j);
    solinit = bvpinit(linspace(0,1,11),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec23(j) = y(7,1); mu = y(7,1);
    Fvec23(j) = -2*pi*y(9,1);
    Evec23(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if mod(j,2) == 1
%        plot(z,r,'color',col5,'linewidth',2);
    end
    y_guess = y;
    if abs(X - 1.3148/2) < 0.00001
%         shape5s = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5s = centercrop(shape5s,0.2,0.2); imwrite(shape5s,'shape5st.png'); 
    end
    
    
end


Xvec23 = [Xvec18(end) Xvec23];
Fvec23 = [Fvec18(end) Fvec23];
muvec23 = [muvec18(end) muvec23];

disp('branch 4 extra')
X = 1.4/2; solinit = bvpinit(linspace(0,1,9),[pi/2 2 0.14 0.1 0 pi -70 1.3 90]); %n=5T
options = bvpset('RelTol',1e-9);
sol = bvp4c(@odesystem,@bcs,solinit);
t = linspace(0,1,N);
y_guess = deval(sol,t);

figure(1);
Xvec24 = [1.415 1.2989]/2;
muvec24 = zeros(1,length(Xvec24));
Fvec24 = zeros(1,length(Xvec24));
Evec24 = zeros(1,length(Xvec24));
%figure; hold on; axis equal; box on;
for j = 1:length(Xvec24)
    X = Xvec24(j);
    solinit = bvpinit(linspace(0,1,9),@newguess);
    sol = bvp4c(@odesystem,@bcs,solinit,options);
    t = linspace(0,1,N);
    y = deval(sol,t);
    r = y(3,:); z = y(4,:); psi_s = y(2,:); psi = y(1,:);
    b = r(51); L = y(8,1);
    muvec24(j) = y(7,1); mu = y(7,1);
    %b = r((N+1)/2); k1 = psi_s((N+1)/2);
    %Fvec(j) = 2*pi*b*y(7,1) + 2*pi*kappa*(1-b^2*k1^2)/(2*b);
    Fvec24(j) = -2*pi*y(9,1);
    Evec24(j) = 2*pi*abs(L)*trapz(t,kappa/2*(psi_s+sin(psi)./r).^2.*r+kbar*psi_s.*sin(psi));
    hold on;  grid on;
    H = (psi_s + sin(psi)./r)/2;
    if  j == length(Xvec24)
%            plot(z,r,'color',col5,'LineWidth',2);
%         shape5s = snapshot(0.75,1.35,psi,psi_s,z,r,X,kappa,N);
%         shape5s = centercrop(shape5s,0.2,0.2); imwrite(shape5s,'shape5st.png'); 
    end
    y_guess = y;
    
    
end

return

set(gca,'FontSize',18);
set(gcf,'color','w');
close all;

figure()
lw = 3;
subplot('position',[0.1 0.1 0.25 0.85])
Xvec2s = Xvec2(Xvec2>= 0.11); Xvec2u = Xvec2(Xvec2 < 0.11);
Fvec2s = Fvec2(Xvec2>= 0.11); Fvec2u = Fvec2(Xvec2 < 0.11);
% plot(2*Xvec2,Fvec2,'linewidth',lw,'color',col1); hold on;
plot(2*Xvec2s,Fvec2s,'linewidth',lw,'color',col1); hold on;
plot(2*Xvec2u,Fvec2u,':','linewidth',lw,'color',col1); hold on;
plot(2*Xvec3,Fvec3,'linewidth',lw,'color',col2); hold on;
plot(2*Xvec4,Fvec4,'linewidth',lw,'color',col2); hold on;
% plot(2*Xvec5,Fvec5,'linewidth',lw,'color',col2); hold on;
Xvec5s = Xvec5(Xvec5>= 0.71); Xvec5u = Xvec5(Xvec5 < 0.71);
Fvec5s = Fvec5(Xvec5>= 0.71); Fvec5u = Fvec5(Xvec5 < 0.71);
plot(2*Xvec5s,Fvec5s,'linewidth',lw,'color',col2); hold on;
plot(2*Xvec5u,Fvec5u,':','linewidth',lw,'color',col2); hold on;
plot(2*Xvec6,Fvec6,':','linewidth',lw,'color',col2); hold on;
plot(2*Xvec7,Fvec7,':','linewidth',lw,'color',col3); hold on;
plot(2*Xvec8,Fvec8,':','linewidth',lw,'color',col3); hold on;
plot(2*Xvec9,Fvec9,':','linewidth',lw,'color',col3); hold on;
% Xvec9s = Xvec9(Xvec9>= 0.71); Xvec9u = Xvec9(Xvec9 < 0.71);
% Fvec9s = Fvec9(Xvec9>= 0.71); Fvec9u = Fvec9(Xvec9 < 0.71);
% plot(2*Xvec9s,Fvec9s,'linewidth',lw,'color',col3); hold on;
% plot(2*Xvec9u,Fvec9u,':','linewidth',lw,'color',col3); hold on;
plot(2*Xvec10,Fvec10,':','linewidth',lw,'color',col3); hold on;
plot(2*Xvec11,Fvec11,':','linewidth',lw,'color',col3); hold on;
plot(2*Xvec12,Fvec12,':','linewidth',lw,'color',col3); hold on;
% plot(2*Xvec13,Fvec13,'linewidth',lw,'color',col3); hold on;
Xvec13s = Xvec13(Xvec13<= 0.71); Xvec13u = Xvec13(Xvec13 > 0.71);
Fvec13s = Fvec13(Xvec13<= 0.71); Fvec13u = Fvec13(Xvec13 > 0.71);
plot(2*Xvec13s,Fvec13s,'linewidth',lw,'color',col3); hold on;
plot(2*Xvec13u,Fvec13u,':','linewidth',lw,'color',col3); hold on;
% plot(2*Xvec14,Fvec14,'linewidth',lw,'color',col4); hold on;
Xvec14s = Xvec14(Xvec14<= 0.61); Xvec14u = Xvec14(Xvec14 > 0.61);
Fvec14s = Fvec14(Xvec14<= 0.61); Fvec14u = Fvec14(Xvec14 > 0.61);
plot(2*Xvec14s,Fvec14s,'linewidth',lw,'color',col4); hold on;
plot(2*Xvec14u,Fvec14u,':','linewidth',lw,'color',col4); hold on;
plot([1.3 1.22],[-25.26 -17.04],':','linewidth',lw, 'color',col4)
plot(2*Xvec15,Fvec15,'linewidth',lw,'color',col4); hold on;
plot(2*Xvec16,Fvec16,':','linewidth',lw,'color',col4); hold on;
plot(2*Xvec17,Fvec17,':','linewidth',lw,'color',col4); hold on;
plot(2*Xvec18,Fvec18,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec19,Fvec19,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec20,Fvec20,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec21,Fvec21,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec22,Fvec22,':','linewidth',lw,'color',col5); hold on;
plot(2*Xvec23,Fvec23,':','linewidth',lw,'color',col5); hold on;
xlabel('$h/a$','interpreter','latex'); ylabel('$Fa/\kappa$','interpreter','latex');  grid on; box on;
set(gca,'FontSize',18);
set(gcf,'color','w');
xlim([0 1.6])
ylim([-1850 50]) 

plotpoint(0,-119.3,'1O',col1,-pi/2,0.08,30);
plotpoint(0.5,-104.5,'1P',col1,-pi/2,0.08,30);
plotpoint(1.0,-64.9,'1Q',col1,-pi/2,0.08,30);
plotpoint(1.5,1.56,'1R',col1,pi/2,0.08,30);
plotpoint(0,-203.3,'2O',col2,-pi/2,0.08,30);
plotpoint(0.5,-229.3,'2P',col2,-pi/2,0.08,30);
plotpoint(1,-155.5,'2Q',col2,-pi/2,0.08,30);
% plotdiamond(1.424,-32.79,'2R',col2,pi/2,0.07,30);
plotpoint(0.94,-3.166,'2s',col2,pi/2,0.08,30);
plotpoint(0,-633.7,'3O',col3,pi/2,0.08,30);
plotpoint(0.5,-556.5,'3P',col3,-pi/2,0.08,30);
plotpoint(1,-375.8,'3Q',col3,-pi/2,0.08,30);
% plotdiamond(1.431,-40.8,'3R',col3,-pi/4,0.08,30);
plotpoint(1.28,-5.845,'3s',col3,pi/2,0.08,30);
plotpoint(0,-669.6,'3T',col3,-pi/2,0.08,30);
plotpoint(0.5,-486.7,'3U',col3,pi/2,0.08,30);
plotpoint(1,-312.2,'3V',col3,pi/2,0.08,30);
% plotdiamond(1.415,-77.73,'3W',col3,-pi/2,0.08,30);
plotpoint(1,-11.42,'3x',col3,-pi/2,0.08,30);
plotpoint(0.5,0.157,'3y',col3,-pi/2,0.08,30);
plotpoint(0,8.55,'3z',col3,-pi/2,0.08,30);
plotpoint(0,-1144,'4O',col4,-pi/2,0.08,30);
plotpoint(0.5,-916.9,'4P',col4,-pi/2,0.08,30);
plotpoint(1,-610.3,'4Q',col4,-pi/2,0.08,30);
% plotdiamond(1.418,-125.3,'4R',col4,-pi/2,0.08,30);
plotpoint(1.12,-12.17,'4s',col4,-pi/2,0.08,30);
plotpoint(0,-1813,'5T',col5,pi/2,0.09,30);
plotpoint(0.5,-1397,'5U',col5,pi/2,0.09,30);
plotpoint(1,-918.4,'5V',col5,pi/2,0.09,30);
% plotdiamond(1.415,-201.4,'5R',col5,pi,0.07,30);
plotpoint(1.315,-51.95,'5s',col5,-pi/2,0.08,30);
plotpoint(0,-1824,'5O',col5,-pi/2,0.09,30);
plotpoint(0.5,-1465,'5P',col5,-pi/2,0.09,30);
plotpoint(1,-983,'5Q',col5,-pi/2,0.09,30);
% plotdiamond(1.417,-231.2,'5W',col5,0,0.07,30);
plotpoint(1.225,-22.2,'5x',col5,-pi/2,0.08,30);

% box
plot([1.38 1.45 1.45 1.38 1.38],[-300 -300 0 0 -300],'k--','linewidth',2)

% inset
axes('position',[.22 .13 .12 .2]); box on;
leftlim = 1.35;
indexOfInterest2 = (2*Xvec2 >= leftlim) & (2*Xvec2 < 1.45) & (Fvec2 >=-300 );
indexOfInterest3 = (2*Xvec3 >= leftlim) & (2*Xvec3 < 1.45) & (Fvec3 >=-300 );
indexOfInterest4 = (2*Xvec4 >= leftlim) & (2*Xvec4 < 1.45) & (Fvec4 >=-300 );
indexOfInterest5s = (2*Xvec5s >= leftlim) & (2*Xvec5s < 1.45) & (Fvec5s >=-300 );
indexOfInterest5u = (2*Xvec5u >= leftlim) & (2*Xvec5u < 1.45) & (Fvec5u >=-300 );
indexOfInterest6 = (2*Xvec6 >= leftlim) & (2*Xvec6 < 1.45) & (Fvec6 >=-300 );
indexOfInterest7 = (2*Xvec7 >= leftlim) & (2*Xvec7 < 1.45) & (Fvec7 >=-300 );
indexOfInterest8 = (2*Xvec8 >= leftlim) & (2*Xvec8 < 1.45) & (Fvec8 >=-300 );
indexOfInterest9 = (2*Xvec9 >= leftlim) & (2*Xvec9 < 1.45) & (Fvec9 >=-300 );
indexOfInterest10 = (2*Xvec10 >= leftlim) & (2*Xvec10 < 1.45) & (Fvec10 >=-300 );
indexOfInterest11 = (2*Xvec11 >= leftlim) & (2*Xvec11 < 1.45) & (Fvec11 >=-300 );
indexOfInterest12 = (2*Xvec12 >= leftlim) & (2*Xvec12 < 1.45) & (Fvec12 >=-300 );
indexOfInterest13s = (2*Xvec13s >= leftlim) & (2*Xvec13s < 1.45) & (Fvec13s >=-300 );
indexOfInterest13u = (2*Xvec13u >= leftlim) & (2*Xvec13u < 1.45) & (Fvec13u >=-300 );
indexOfInterest14 = (2*Xvec14 >= leftlim) & (2*Xvec14 < 1.45) & (Fvec14 >=-300 );
indexOfInterest15 = (2*Xvec15 >= leftlim) & (2*Xvec15 < 1.45) & (Fvec15 >=-350 );
indexOfInterest16 = (2*Xvec16 >= leftlim) & (2*Xvec16 < 1.45) & (Fvec16 >=-350 );
indexOfInterest17 = (2*Xvec17 >= leftlim) & (2*Xvec17 < 1.45) & (Fvec17 >=-350 );
indexOfInterest18 = (2*Xvec18 >= leftlim) & (2*Xvec18 < 1.45) & (Fvec18 >=-350 );
indexOfInterest19 = (2*Xvec19 >= leftlim) & (2*Xvec19 < 1.45) & (Fvec19 >=-350 );
indexOfInterest20 = (2*Xvec20 >= leftlim) & (2*Xvec20 < 1.45) & (Fvec20 >=-350 );
indexOfInterest21 = (2*Xvec21 >= leftlim) & (2*Xvec21 < 1.45) & (Fvec21 >=-350 );
indexOfInterest22 = (2*Xvec22 >= leftlim) & (2*Xvec22 < 1.45) & (Fvec22 >=-350 );
indexOfInterest23 = (2*Xvec23 >= leftlim) & (2*Xvec23 < 1.45) & (Fvec23 >=-350 );
% plot(Xvec2(indexOfInterest2),Fvec2(indexOfInterest2),'linewidth',lw,'color',col1); hold on; grid on;
plot(2*Xvec3(indexOfInterest3),Fvec3(indexOfInterest3),'linewidth',lw,'color',col2); hold on; grid on;
plot(2*Xvec4(indexOfInterest4),Fvec4(indexOfInterest4),'linewidth',lw,'color',col2);
plot(2*Xvec5u(indexOfInterest5u),Fvec5u(indexOfInterest5u),':','linewidth',lw,'color',col2);
plot(2*Xvec5s(indexOfInterest5s),Fvec5s(indexOfInterest5s),'linewidth',lw,'color',col2);
plot([1.418 1.42],[-39.81 -38.11],':','linewidth',lw,'color',col2);
plot(2*Xvec6(indexOfInterest6),Fvec6(indexOfInterest6),':','linewidth',lw,'color',col2);
plot(2*Xvec7(indexOfInterest7),Fvec7(indexOfInterest7),':','linewidth',lw,'color',col3);
plot(2*Xvec8(indexOfInterest8),Fvec8(indexOfInterest8),':','linewidth',lw,'color',col3);
plot(2*Xvec9(indexOfInterest9),Fvec9(indexOfInterest9),':','linewidth',lw,'color',col3);
plot(2*Xvec10(indexOfInterest10),Fvec10(indexOfInterest10),':','linewidth',lw,'color',col3);
plot(2*Xvec11(indexOfInterest11),Fvec11(indexOfInterest11),':','linewidth',lw,'color',col3);
plot(2*Xvec12(indexOfInterest12),Fvec12(indexOfInterest12),':','linewidth',lw,'color',col3);
plot(2*Xvec13s(indexOfInterest13s),Fvec13s(indexOfInterest13s),'linewidth',lw,'color',col3);
plot(2*Xvec13u(indexOfInterest13u),Fvec13u(indexOfInterest13u),':','linewidth',lw,'color',col3);
plot(2*Xvec14(indexOfInterest14),Fvec14(indexOfInterest14),':','linewidth',lw,'color',col4);
plot(2*Xvec15(indexOfInterest15),Fvec15(indexOfInterest15),'linewidth',lw,'color',col4);
plot(2*Xvec16(indexOfInterest16),Fvec16(indexOfInterest16),':','linewidth',lw,'color',col4);
plot(2*Xvec17(indexOfInterest17),Fvec17(indexOfInterest17),':','linewidth',lw,'color',col4);
plot(2*Xvec18(indexOfInterest18),Fvec18(indexOfInterest18),':','linewidth',lw,'color',col5);
plot(2*Xvec19(indexOfInterest19),Fvec19(indexOfInterest19),':','linewidth',lw,'color',col5);
plot(2*Xvec20(indexOfInterest20),Fvec20(indexOfInterest20),':','linewidth',lw,'color',col5);
plot(2*Xvec21(indexOfInterest21),Fvec21(indexOfInterest21),':','linewidth',lw,'color',col5);
plot(2*Xvec22(indexOfInterest22),Fvec22(indexOfInterest22),':','linewidth',lw,'color',col5);
plot(2*Xvec23(indexOfInterest23),Fvec23(indexOfInterest23),':','linewidth',lw,'color',col5);
plotpoint(1.424,-32.79,'2R',col2,0,0.003,20);
plotpoint(1.4312,-40.8,'3R',col3,0,0.003,20);
plotpoint(1.415,-77.73,'3W',col3,0,0.003,20);
plotpoint(1.4182,-130.3,'4R',col4,0,0.003,20);
plotpoint(1.416,-232.7,'5R',col5,pi,0.003,20);
plotpoint(1.4178,-218.2,'5W',col5,0,0.003,20);

axis([1.4 1.435 -300 0]);

% other inset
axes('position',[.28 .35 .06 .1]); box on; 
indexOfInterest = (2*Xvec2 >=1.4 ) & (2*Xvec2 < 3);
plot(2*Xvec2(indexOfInterest),Fvec2(indexOfInterest),'-','linewidth',lw,'color',col1); hold on; grid on;
plotpoint(2.5,46.59,'1s',col1,-pi/2,0.08,9);
plotpoint(1.5,1.56,'1R',col1,-pi/4,0.08,9);



end

function [] = plotpoint(x,y,txt,col,theta,xl,yl)
plot(x,y,'o','markersize',10,'markerfacecolor',col,'markeredgecolor',col);
% plot(x+0.03*cos(theta),y+50*sin(theta),'o','markersize',20,'markerfacecolor','w','markeredgecolor','k');
text(x+xl*cos(theta),y+yl*sin(theta),txt,'fontsize',14,'horizontalalignment','center');

end

function [] = plotsquare(x,y,txt,col,theta,xl,yl)
plot(x,y,'s','markersize',11,'markerfacecolor',col,'markeredgecolor',col);
% plot(x+0.03*cos(theta),y+50*sin(theta),'o','markersize',20,'markerfacecolor','w','markeredgecolor','k');
text(x+xl*cos(theta),y+yl*sin(theta),txt,'fontsize',14,'horizontalalignment','center');

end

function [] = plotdiamond(x,y,txt,col,theta,xl,yl)
plot(x,y,'d','markersize',11,'markerfacecolor',col,'markeredgecolor',col);
% plot(x+0.03*cos(theta),y+50*sin(theta),'o','markersize',20,'markerfacecolor','w','markeredgecolor','k');
text(x+xl*cos(theta),y+yl*sin(theta),txt,'fontsize',14,'horizontalalignment','center');

end

function [] = plotpoint2(x,y,txt,col,theta)
plot(x,y,'o','markersize',10,'markerfacecolor',col,'markeredgecolor',col);
% plot(x+0.03*cos(theta),y+50*sin(theta),'o','markersize',20,'markerfacecolor','w','markeredgecolor','k');
text(x+0.07*cos(theta),y+7*sin(theta),txt,'fontsize',14,'horizontalalignment','center');

end

function [] = mergeshapefigs()

    s1A = imread('shape1A.png'); s1A = insertText(s1A,[10,10],'1A','fontsize',50,'boxcolor','w');
    s1B = imread('shape1B.png'); s1B = insertText(s1B,[10,10],'1B','fontsize',50,'boxcolor','w');
    s1C = imread('shape1C.png'); s1C = insertText(s1C,[10,10],'1C','fontsize',50,'boxcolor','w');
    s1D = imread('shape1D.png'); s1D = insertText(s1D,[10,10],'1D','fontsize',50,'boxcolor','w');
    s1E = imread('shape1E.png'); s1E = insertText(s1E,[10,10],'1E','fontsize',50,'boxcolor','w');
    s1F = imread('shape1F.png'); s1F = insertText(s1F,[10,10],'1F','fontsize',50,'boxcolor','w');
    s1G = imread('shape1G.png'); s1G = insertText(s1G,[10,10],'1G','fontsize',50,'boxcolor','w');
    
    s2A = imread('shape2A.png'); s2A = insertText(s2A,[10,10],'2A','fontsize',50,'boxcolor','w');
    s2B = imread('shape2B.png'); s2B = insertText(s2B,[10,10],'2B','fontsize',50,'boxcolor','w');
    s2C = imread('shape2C.png'); s2C = insertText(s2C,[10,10],'2C','fontsize',50,'boxcolor','w');
    s2D = imread('shape1D.png'); s2D = insertText(s2D,[10,10],'2D','fontsize',50,'boxcolor','w');
    s2Cp = imread('shape2Cp.png'); s2Cp = insertText(s2Cp,[10,10],'2C''','fontsize',50,'boxcolor','w');
    s2Bp = imread('shape2Bp.png'); s2Bp = insertText(s2Bp,[10,10],'2B''','fontsize',50,'boxcolor','w');
    s2Ap = imread('shape2Ap.png'); s2Ap = insertText(s2Ap,[10,10],'2A''','fontsize',50,'boxcolor','w');
    
    s3A = imread('shape3A.png'); s3A = insertText(s3A,[10,10],'3A','fontsize',50,'boxcolor','w');
    s3B = imread('shape3B.png'); s3B = insertText(s3B,[10,10],'3B','fontsize',50,'boxcolor','w');
    s3C = imread('shape3C.png'); s3C = insertText(s3C,[10,10],'3C','fontsize',50,'boxcolor','w');
    s3D = imread('shape1D.png'); s3D = insertText(s3D,[10,10],'3D','fontsize',50,'boxcolor','w');
    s3E = imread('shape3E.png'); s3E = insertText(s3E,[10,10],'3E','fontsize',50,'boxcolor','w');
    s3F = imread('shape3F.png'); s3F = insertText(s3F,[10,10],'3F','fontsize',50,'boxcolor','w');
    s3G = imread('shape3G.png'); s3G = insertText(s3G,[10,10],'3G','fontsize',50,'boxcolor','w');

    s4A = imread('shape4A.png'); s4A = insertText(s4A,[10,10],'4A','fontsize',50,'boxcolor','w');
    s4B = imread('shape4B.png'); s4B = insertText(s4B,[10,10],'4B','fontsize',50,'boxcolor','w');
    s4C = imread('shape4C.png'); s4C = insertText(s4C,[10,10],'4C','fontsize',50,'boxcolor','w');
    s4D = imread('shape1D.png'); s4D = insertText(s4D,[10,10],'4D','fontsize',50,'boxcolor','w');
    s4Cp = imread('shape4Cp.png'); s4Cp = insertText(s4Cp,[10,10],'4C''','fontsize',50,'boxcolor','w');
    s4Bp = imread('shape4Bp.png'); s4Bp = insertText(s4Bp,[10,10],'4B''','fontsize',50,'boxcolor','w');
    s4Ap = imread('shape4Ap.png'); s4Ap = insertText(s4Ap,[10,10],'4A''','fontsize',50,'boxcolor','w');    
    
    s5A = imread('shape5A.png'); s5A = insertText(s5A,[10,10],'5A','fontsize',50,'boxcolor','w');
    s5B = imread('shape5B.png'); s5B = insertText(s5B,[10,10],'5B','fontsize',50,'boxcolor','w');
    s5C = imread('shape5C.png'); s5C = insertText(s5C,[10,10],'5C','fontsize',50,'boxcolor','w');
    s5D = imread('shape1D.png'); s5D = insertText(s5D,[10,10],'5D','fontsize',50,'boxcolor','w');
    s5E = imread('shape5E.png'); s5E = insertText(s5E,[10,10],'5E','fontsize',50,'boxcolor','w');
    s5F = imread('shape5F.png'); s5F = insertText(s5F,[10,10],'5F','fontsize',50,'boxcolor','w');
    s5G = imread('shape5G.png'); s5G = insertText(s5G,[10,10],'5G','fontsize',50,'boxcolor','w');    
    
    bigimg = [s1A s1B s1C s1D s1E s1F s1G;
              s2A s2B s2C s2D s2Cp s2Bp s2Ap;
              s3A s3B s3C s3D s3E s3F s3G;
              s4A s4B s4C s4D s4Cp s4Bp s4Ap;
              s5A s5B s5C s5D s5E s5F s5G;];
    imshow(bigimg)
    imwrite(bigimg,'smallareashapes_winter.png')
               
end

function [] = mergeshapefigs2thick()

    s1H = imread('shape1H.png'); s1H = insertText(s1H,[10,10],'1H','fontsize',50,'boxcolor','w');
    s1I = imread('shape1I.png'); s1I = insertText(s1I,[10,10],'1I','fontsize',50,'boxcolor','w');
    s1J = imread('shape1J.png'); s1J = insertText(s1J,[10,10],'1J','fontsize',50,'boxcolor','w');
    s1K = imread('shape1K.png'); s1K = insertText(s1K,[10,10],'1K','fontsize',50,'boxcolor','w');
    s1l = imread('shape1lt.png'); s1l = insertText(s1l,[10,10],['1' char(955)],'fontsize',50,'boxcolor','w');
    s1m = imread('shape1mt.png'); s1m = insertText(s1m,[10,10],'1m','fontsize',50,'boxcolor','w');
    s1n = imread('shape1nt.png'); s1n = insertText(s1n,[10,10],'1n','fontsize',50,'boxcolor','w');
    
    s2H = imread('shape2H.png'); s2H = insertText(s2H,[10,10],'2H','fontsize',50,'boxcolor','w');
    s2I = imread('shape2I.png'); s2I = insertText(s2I,[10,10],'2I','fontsize',50,'boxcolor','w');
    s2J = imread('shape2J.png'); s2J = insertText(s2J,[10,10],'2J','fontsize',50,'boxcolor','w');
    s2K = imread('shape1K.png'); s2K = insertText(s2K,[10,10],'2K','fontsize',50,'boxcolor','w');
    s2Jp = imread('shape2Jp.png'); s2Jp = insertText(s2Jp,[10,10],'2J''','fontsize',50,'boxcolor','w');
    s2Ip = imread('shape2Ip.png'); s2Ip = insertText(s2Ip,[10,10],'2I''','fontsize',50,'boxcolor','w');
    s2Hp = imread('shape2Hp.png'); s2Hp = insertText(s2Hp,[10,10],'2H''','fontsize',50,'boxcolor','w');
    
    s3H = imread('shape3H.png'); s3H = insertText(s3H,[10,10],'3H','fontsize',50,'boxcolor','w');
    s3I = imread('shape3I.png'); s3I = insertText(s3I,[10,10],'3I','fontsize',50,'boxcolor','w');
    s3J = imread('shape3J.png'); s3J = insertText(s3J,[10,10],'3J','fontsize',50,'boxcolor','w');
    s3K = imread('shape1K.png'); s3K = insertText(s3K,[10,10],'3K','fontsize',50,'boxcolor','w');
    s3L = imread('shape3L.png'); s3L = insertText(s3L,[10,10],'3L','fontsize',50,'boxcolor','w');
    s3M = imread('shape3M.png'); s3M = insertText(s3M,[10,10],'3M','fontsize',50,'boxcolor','w');
    s3N = imread('shape3N.png'); s3N = insertText(s3N,[10,10],'3N','fontsize',50,'boxcolor','w');
 
    s4H = imread('shape4H.png'); s4H = insertText(s4H,[10,10],'4H','fontsize',50,'boxcolor','w');
    s4I = imread('shape4I.png'); s4I = insertText(s4I,[10,10],'4I','fontsize',50,'boxcolor','w');
    s4J = imread('shape4J.png'); s4J = insertText(s4J,[10,10],'4J','fontsize',50,'boxcolor','w');
    s4K = imread('shape1K.png'); s4K = insertText(s4K,[10,10],'4K','fontsize',50,'boxcolor','w');
    s4Jp = imread('shape4Jp.png'); s4Jp = insertText(s4Jp,[10,10],'4J''','fontsize',50,'boxcolor','w');
    s4Ip = imread('shape4Ip.png'); s4Ip = insertText(s4Ip,[10,10],'4I''','fontsize',50,'boxcolor','w');
    s4Hp = imread('shape4Hp.png'); s4Hp = insertText(s4Hp,[10,10],'4H''','fontsize',50,'boxcolor','w');  
     
    s5H = imread('shape5H.png'); s5H = insertText(s5H,[10,10],'5H','fontsize',50,'boxcolor','w');
    s5I = imread('shape5I.png'); s5I = insertText(s5I,[10,10],'5I','fontsize',50,'boxcolor','w');
    s5J = imread('shape5J.png'); s5J = insertText(s5J,[10,10],'5J','fontsize',50,'boxcolor','w');
    s5K = imread('shape1K.png'); s5K = insertText(s5K,[10,10],'5K','fontsize',50,'boxcolor','w');
    s5L = imread('shape5L.png'); s5L = insertText(s5L,[10,10],'5L','fontsize',50,'boxcolor','w');
    s5M = imread('shape5M.png'); s5M = insertText(s5M,[10,10],'5M','fontsize',50,'boxcolor','w');
    s5N = imread('shape5N.png'); s5N = insertText(s5N,[10,10],'5N','fontsize',50,'boxcolor','w');    
    
     
    bigimg = [s1H s1I s1J s1K s1l s1m s1n;
              s2H s2I s2J s2K s2Jp s2Ip s2Hp;
              s3H s3I s3J s3K s3L s3M s3N;
              s4H s4I s4J s4K s4Jp s4Ip s4Hp;
              s5H s5I s5J s5K s5L s5M s5N];
    imshow(bigimg)
    imwrite(bigimg,'mediumareashapes_thick_winter.png')
               
end

function [] = mergeshapefigs2thin()
    
    s2i = imread('shape2it.png'); s2i = insertText(s2i,[10,10],'2i','fontsize',50,'boxcolor','w');
    s2j = imread('shape2jt.png'); s2j = insertText(s2j,[10,10],'2j','fontsize',50,'boxcolor','w');
    s2k = imread('shape1lt.png'); s2k = insertText(s2k,[10,10],'2k','fontsize',50,'boxcolor','w');
    s2jp = imread('shape2jpt.png'); s2jp = insertText(s2jp,[10,10],'2j''','fontsize',50,'boxcolor','w');
    s2ip = imread('shape2ipt.png'); s2ip = insertText(s2ip,[10,10],'2i''','fontsize',50,'boxcolor','w');    
    
    s3h = imread('shape3ht.png'); s3h = insertText(s3h,[10,10],'3h','fontsize',50,'boxcolor','w');
    s3i = imread('shape3it.png'); s3i = insertText(s3i,[10,10],'3i','fontsize',50,'boxcolor','w');
    s3j = imread('shape3jt.png'); s3j = insertText(s3j,[10,10],'3j','fontsize',50,'boxcolor','w');
    s3k = imread('shape1lt.png'); s3k = insertText(s3k,[10,10],'3k','fontsize',50,'boxcolor','w');
    s3l = imread('shape3lt.png'); s3l = insertText(s3l,[10,10],['3' char(955)],'fontsize',50,'boxcolor','w');
    s3m = imread('shape3mt.png'); s3m = insertText(s3m,[10,10],'3m','fontsize',50,'boxcolor','w'); 
    
    s4i = imread('shape4it.png'); s4i = insertText(s4i,[10,10],'4i','fontsize',50,'boxcolor','w');
    s4j = imread('shape4jt.png'); s4j = insertText(s4j,[10,10],'4j','fontsize',50,'boxcolor','w');
    s4k = imread('shape1lt.png'); s4k = insertText(s4k,[10,10],'4k','fontsize',50,'boxcolor','w');
    s4jp = imread('shape4jpt.png'); s4jp = insertText(s4jp,[10,10],'4j''','fontsize',50,'boxcolor','w');
    s4ip = imread('shape4ipt.png'); s4ip = insertText(s4ip,[10,10],'4i''','fontsize',50,'boxcolor','w');  

    s5h = imread('shape5ht.png'); s5h = insertText(s5h,[10,10],'5h','fontsize',50,'boxcolor','w');
    s5i = imread('shape5it.png'); s5i = insertText(s5i,[10,10],'5i','fontsize',50,'boxcolor','w');
    s5j = imread('shape5jt.png'); s5j = insertText(s5j,[10,10],'5j','fontsize',50,'boxcolor','w');
    s5k = imread('shape1lt.png'); s5k = insertText(s5k,[10,10],'5k','fontsize',50,'boxcolor','w');
    s5l = imread('shape5lt.png'); s5l = insertText(s5l,[10,10],['5' char(955)],'fontsize',50,'boxcolor','w');
    s5m = imread('shape5mt.png'); s5m = insertText(s5m,[10,10],'5m','fontsize',50,'boxcolor','w'); 

    wph = 255*ones(size(s2i));
    bigimg =  [wph s2i s2j s2k s2jp s2ip ;
              s3h s3i s3j s3k s3l s3m ;
              wph s4i s4j s4k s4jp s4ip ;
              s5h s5i s5j s5k s5l s5m ];
    imshow(bigimg)
    imwrite(bigimg,'mediumareashapes_thin_winter.png')
end

function [] = mergeshapefigs3()

    s1O = imread('shape1O.png'); s1O = insertText(s1O,[10,10],'1O','fontsize',50,'boxcolor','w');
    s1P = imread('shape1P.png'); s1P = insertText(s1P,[10,10],'1P','fontsize',50,'boxcolor','w');
    s1Q = imread('shape1Q.png'); s1Q = insertText(s1Q,[10,10],'1Q','fontsize',50,'boxcolor','w');
    s1R = imread('shape1R.png'); s1R = insertText(s1R,[10,10],'1R','fontsize',50,'boxcolor','w');
    s1s = imread('shape1st.png'); s1s = insertText(s1s,[10,10],'1s','fontsize',50,'boxcolor','w');
    
    %watch filenames -- some switching needed
    s2O = imread('shape2O.png'); s2O = insertText(s2O,[10,10],'2O','fontsize',50,'boxcolor','w');
    s2P = imread('shape2P.png'); s2P = insertText(s2P,[10,10],'2P','fontsize',50,'boxcolor','w');
    s2Q = imread('shape2Q.png'); s2Q = insertText(s2Q,[10,10],'2Q','fontsize',50,'boxcolor','w');
    s2R = imread('shape2Rp.png'); s2R = insertText(s2R,[10,10],'2R','fontsize',50,'boxcolor','w');
    s2Rp = imread('shape2R.png'); s2Rp = insertText(s2Rp,[10,10],'2R''','fontsize',50,'boxcolor','w');
    s2Qp = imread('shape2Qp.png'); s2Qp = insertText(s2Qp,[10,10],'2Q''','fontsize',50,'boxcolor','w');
    s2Pp = imread('shape2Pp.png'); s2Pp = insertText(s2Pp,[10,10],'2P''','fontsize',50,'boxcolor','w');
    s2Op = imread('shape2Op.png'); s2Op = insertText(s2Op,[10,10],'2O''','fontsize',50,'boxcolor','w');
    s2s = imread('shape2spt.png'); s2s = insertText(s2s,[10,10],'2s','fontsize',50,'boxcolor','w');
    s2sp = imread('shape2st.png'); s2sp = insertText(s2sp,[10,10],'2s''','fontsize',50,'boxcolor','w');
    
    s3O = imread('shape3O.png'); s3O = insertText(s3O,[10,10],'3O','fontsize',50,'boxcolor','w');
    s3P = imread('shape3P.png'); s3P = insertText(s3P,[10,10],'3P','fontsize',50,'boxcolor','w');
    s3Q = imread('shape3Q.png'); s3Q = insertText(s3Q,[10,10],'3Q','fontsize',50,'boxcolor','w');
    s3R = imread('shape3R.png'); s3R = insertText(s3R,[10,10],'3R','fontsize',50,'boxcolor','w');
    s3s = imread('shape3st.png'); s3s = insertText(s3s,[10,10],'3s','fontsize',50,'boxcolor','w');
    
    s3T = imread('shape3T.png'); s3T = insertText(s3T,[10,10],'3T','fontsize',50,'boxcolor','w');
    s3U = imread('shape3U.png'); s3U = insertText(s3U,[10,10],'3U','fontsize',50,'boxcolor','w');
    s3V = imread('shape3V.png'); s3V = insertText(s3V,[10,10],'3V','fontsize',50,'boxcolor','w');
    s3W = imread('shape3W.png'); s3W = insertText(s3W,[10,10],'3W','fontsize',50,'boxcolor','w');
    s3x = imread('shape3xt.png'); s3x = insertText(s3x,[10,10],'3x','fontsize',50,'boxcolor','w');
    s3y = imread('shape3yt.png'); s3y = insertText(s3y,[10,10],'3y','fontsize',50,'boxcolor','w');
    s3z = imread('shape3zt.png'); s3z = insertText(s3z,[10,10],'3z','fontsize',50,'boxcolor','w');
    
    s4O = imread('shape4O.png'); s4O = insertText(s4O,[10,10],'4O','fontsize',50,'boxcolor','w');
    s4P = imread('shape4P.png'); s4P = insertText(s4P,[10,10],'4P','fontsize',50,'boxcolor','w');
    s4Q = imread('shape4Q.png'); s4Q = insertText(s4Q,[10,10],'4Q','fontsize',50,'boxcolor','w');
    s4R = imread('shape4R.png'); s4R = insertText(s4R,[10,10],'4R','fontsize',50,'boxcolor','w');
    s4Rp = imread('shape4Rp.png'); s4Rp = insertText(s4Rp,[10,10],'4R''','fontsize',50,'boxcolor','w');
    s4Qp = imread('shape4Qp.png'); s4Qp = insertText(s4Qp,[10,10],'4Q''','fontsize',50,'boxcolor','w');
    s4Pp = imread('shape4Pp.png'); s4Pp = insertText(s4Pp,[10,10],'4P''','fontsize',50,'boxcolor','w');
    s4Op = imread('shape4Op.png'); s4Op = insertText(s4Op,[10,10],'4O''','fontsize',50,'boxcolor','w');
    s4s = imread('shape4st.png'); s4s = insertText(s4s,[10,10],'4s','fontsize',50,'boxcolor','w');
    s4sp = imread('shape4spt.png'); s4sp = insertText(s4sp,[10,10],'4s''','fontsize',50,'boxcolor','w');

    s5O = imread('shape5O.png'); s5O = insertText(s5O,[10,10],'5O','fontsize',50,'boxcolor','w');
    s5P = imread('shape5P.png'); s5P = insertText(s5P,[10,10],'5P','fontsize',50,'boxcolor','w');
    s5Q = imread('shape5Q.png'); s5Q = insertText(s5Q,[10,10],'5Q','fontsize',50,'boxcolor','w');
    s5R = imread('shape5R.png'); s5R = insertText(s5R,[10,10],'5R','fontsize',50,'boxcolor','w');
    s5s = imread('shape5st.png'); s5s = insertText(s5s,[10,10],'5s','fontsize',50,'boxcolor','w');
    
    s5T = imread('shape5T.png'); s5T = insertText(s5T,[10,10],'5T','fontsize',50,'boxcolor','w');
    s5U = imread('shape5U.png'); s5U = insertText(s5U,[10,10],'5U','fontsize',50,'boxcolor','w');
    s5V = imread('shape5V.png'); s5V = insertText(s5V,[10,10],'5V','fontsize',50,'boxcolor','w');
    s5W = imread('shape5W.png'); s5W = insertText(s5W,[10,10],'5W','fontsize',50,'boxcolor','w');
    s5x = imread('shape5xt.png'); s5x = insertText(s5x,[10,10],'5x','fontsize',50,'boxcolor','w');
    
    wph = 255*ones(size(s1O));
    bigimg = [s1O s1P s1Q s1R s1s wph wph;
              s2O s2P s2Q s2R s2s wph wph;
              s2Op s2Pp s2Qp s2Rp s2sp wph wph;
              s3O s3P s3Q s3R s3s wph wph;
              s3T s3U s3V s3W s3x s3y s3z
              s4O s4P s4Q s4R s4s wph wph;
              s4Op s4Pp s4Qp s4Rp s4sp wph wph;
              s5O s5P s5Q s5R s5s wph wph;
              s5T s5U s5V s5W s5x wph wph];
    imshow(bigimg)
    imwrite(bigimg,'largeareashapes_winter.png');
%     fig2 = figure(); axis off;
%     cb = colorbar('southoutside'); cb.Label.String = 'mean curvature aH'; colormap('winter'); caxis([-40 40])
%     set(gcf,'color','w'); 
%     cbimg = frame2im(getframe(fig2)); [r c p] = size(cbimg); cbimg = cbimg(round(0.8*r):end,:,:);
%      [r c p] = size(cbimg);
%     [rbig cbig pbig] = size(bigimg);
%     ws = 255*ones([r ,round((cbig - c)/2),  p]);
%     cbbigimg = [bigimg; 255*ones([r ,round((cbig - c)/2),  p]) cbimg 255*ones([r ,round((cbig - c)/2)-1,  p])];
%     close(fig2);
%     imwrite(cbbigimg,'largeareashapes_winter.png')

end

function [] = cbfig()

    fig2 = figure(); axis off; caxis([-20 20])
    cb = colorbar('southoutside');  colormap('winter');
    cb.Label.String = '$Ha$';
    cb.Label.Interpreter = 'latex';
%     cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 18;
    
    set(cb,'TickLabels',20:-5:-20); %spheres have H < 0 in our convention
    set(cb,'YDir','reverse');
   
    
    set(gcf,'color','w'); 
    cbimg = frame2im(getframe(fig2)); [r c p] = size(cbimg); cbimg = cbimg(round(0.8*r):round(0.95*r),:,:);
    imwrite(cbimg,'colorbar_winter.png');
    disp('image saved')
%     close(fig2);


end


function out = snapshot(width,height,psi,psi_s,z,r,X,kappa,N)
        
         fig = figure;
        %[T,P] = meshgrid([1,1],linspace(0,2*pi,60));
        [T,P] = meshgrid(1:N,linspace(0,2*pi,60));
        C = (psi_s(T) + sin(psi(T))./r(T));
%        C = (psi_s(T) + sin(psi(T))./r(T)).^2*kappa/2;
        %r = [zeros(60,1) ones(60,1)]; z = zeros(60,2); C = zeros(60,2);
        %s=surf(z,r.*cos(P),r.*sin(P),log(C)/log(10)); hold on;
%         s=surf(z(T),r(T).*cos(P),r(T).*sin(P),log(C)/log(10)); hold on;
        s=surf(z(T),r(T).*cos(P),r(T).*sin(P),C); hold on;
        
%         cb = colorbar('southoutside'); cb.Label.String = 'log_{10} (\kappa/2)(2H)^2';
%         caxis([-2.0 2.0])
        colormap('winter')
        caxis([-20 20])
        s.EdgeColor = 'none'; 
        %s.FaceColor = 'interp';
        %X = 0
         plot3(ones([1,100])*(X+0),cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',4);
         plot3(-ones([1,100])*(X+0.01),cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',4);
         
         
         %text(-3*width,0,0,str,'fontsize',40)
         axis equal;
%         axis([-0.6 0.6 -1.35 1.35 -1.35 1.35]); %box on; 
        axis([-width width -height height -height height]); %box on; 
        camorbit(15,0);
        axis off;
        set(gcf,'color','w');set(gca,'FontSize',18);
        
        for q = 1:2:60
            plot3(z,r.*cos(2*pi*q/60),r.*sin(2*pi*q/60),'k');
            %plot3([0 0],[0,1*cos(2*pi*q/60)],[0,1*sin(2*pi*q/60)],'k')
        end
     
       frame = getframe(fig);
       out = frame2im(frame);

end

function [] = snapshotalt(width,height,psi,psi_s,z,r,X,kappa,N)
        
%          fig = figure;
        %[T,P] = meshgrid([1,1],linspace(0,2*pi,60));
        [T,P] = meshgrid(1:N,linspace(0,2*pi,60));
        C = (psi_s(T) + sin(psi(T))./r(T));
%        C = (psi_s(T) + sin(psi(T))./r(T)).^2*kappa/2;
        %r = [zeros(60,1) ones(60,1)]; z = zeros(60,2); C = zeros(60,2);
        %s=surf(z,r.*cos(P),r.*sin(P),log(C)/log(10)); hold on;
%         s=surf(z(T),r(T).*cos(P),r(T).*sin(P),log(C)/log(10)); hold on;
        s=surf(z(T),r(T).*cos(P),r(T).*sin(P),C); hold on;
        
%         cb = colorbar('southoutside'); cb.Label.String = 'log_{10} (\kappa/2)(2H)^2';
%         caxis([-2.0 2.0])
        colormap('winter')
        caxis([-20 20])
        s.EdgeColor = 'none'; 
        %s.FaceColor = 'interp';
        %X = 0
         plot3(ones([1,100])*(X+0),cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',1);
         plot3(-ones([1,100])*(X+0.01),cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',1);
         
         
         %text(-3*width,0,0,str,'fontsize',40)
         axis equal;
%         axis([-0.6 0.6 -1.35 1.35 -1.35 1.35]); %box on; 
        axis([-width width -height height -height height]); %box on; 
        camorbit(15,0);
        axis off;
        set(gcf,'color','w');set(gca,'FontSize',18);
        
        for q = 1:2:60
            plot3(z,r.*cos(2*pi*q/60),r.*sin(2*pi*q/60),'k');
            %plot3([0 0],[0,1*cos(2*pi*q/60)],[0,1*sin(2*pi*q/60)],'k')
        end
%      
%        frame = getframe(fig);
%        out = frame2im(frame);

end

%% cut whitespace 840 x 1120 px image 
function out = centercrop(img,wpercent,hpercent)
    [r,c,p] = size(img);
    out = img(round(wpercent*r):round((1-wpercent)*r),round(hpercent*c):round((1-hpercent)*c),:);
end
