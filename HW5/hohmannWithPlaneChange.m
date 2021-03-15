function hohmannWithPlaneChange(a0,af,inc0,incf,Omega0,Omegaf,mu,animate)
% ------------------------------------------------------------------- %
%                            hohmannWithPlaneChange.m                 %
% ------------------------------------------------------------------- %
% This code implements a Hohmann-type transfer with an inclination    %
% change to transfer a spacecraft between two circular non-coplanar   %
% orbits.  The longitudes of ascending node for each orbit may in     %
% general be different.  Thus, the transfer starts and terminates on  %
% the line of intersection of the two orbits.                         %
% ------------------------------------------------------------------- %
%   Inputs:                                                           %
%   a0:      initial semi-major axis                         (km)     %
%   af:      terminal semi-major axis                        (km)     %
%   inc0:    initial inclination                            (rad)     %
%   incf:    terminal inclination                           (rad)     %
%   Omega0:  initial longitude of ascending node            (rad)     %
%   Omegaf:  terminal longitude of ascending node           (rad)     %
%   mu:      planet gravitational parameter            (km^3/s^2)     %
%   animate: enable/disable animation of transfer orbit (boolean)     %
% ------------------------------------------------------------------- %
%   Outputs:  visualization of orbit transfer                         %
% ------------------------------------------------------------------- %
%   Exceptional Cases:                                                %
%      (1) Omega0 = Omegaf: line of intersection is set to lie along  %
%                           the line of nodes of the initial orbit.   %
%      (2) inc0=incf~=0:    line of intersection is set to lie along  %
%                           the line of nodes of the initial orbit.   %
%      (3) inc0=incf=0:     line of intersection is set to lie along  %
%                           the direction Ix = [1 0 0].               %
% ------------------------------------------------------------------- %
close all; clc;

if nargin<7,
    error('Not Enough Input Arguments to Function');
elseif isequal(nargin,7),
    fprintf('------------------------------------------------------------------\n');
    fprintf('Animation Flag Is Not An Input.    Animation Will Not Be Displayed\n');
    fprintf('------------------------------------------------------------------\n');
    animate = false;
    animate = false;
elseif ~(isequal(animate,0) | isequal(animate,1)),
    fprintf('------------------------------------------------------------------\n');
    fprintf('Animate Flag Not Equal to 0 or 1.  Animation Will Not Be Displayed\n');
    fprintf('------------------------------------------------------------------\n');
    animate = false;
end
% ----------------------------------------------------------------------------- %
%                                   Physical Constants                          %
% ----------------------------------------------------------------------------- %
Ez          = [0; 0; 1];
deg2rad     = pi/180;
rad2deg     = 1/deg2rad;

% ----------------------------------------------------------------------------- %
%                             Initial Circular Orbit                            %
% ----------------------------------------------------------------------------- %
e0            = 0;
omega0        = 0*deg2rad;
nu0           = 0*deg2rad;
tau0          = 2*pi*sqrt(a0^3/mu);
oe0           = [a0; e0; Omega0; inc0; omega0; nu0];
[rv00,vv00]   = oe2rv_BraceyIV_Edgar(oe0,mu);
hv0           = cross(rv00,vv00);
uhv0          = hv0/norm(hv0,2);
nv0           = cross(Ez,hv0);
tVector       = (0:tau0/2000:tau0).';
rvMatrix      = zeros(length(tVector),3);
vvMatrix      = zeros(length(tVector),3);
rv            = rv00;
vv            = vv00;
rvMatrix(1,:) = rv;
vvMatrix(1,:) = vv;
tPrev         = tVector(1);
for i = 2:length(tVector),
    tCurr                   = tVector(i);
    % [rv,vv,~,~,~,~] = propagateKepler_BraceyIV_Edgar(tPrev,tCurr,rv,vv,mu);
    [rv,vv,~,~,~,~] = propagateKepler_BraceyIV_Edgar(tVector(1),tCurr,rv00,vv00,mu);
    rvMatrix(i,:)           = rv.';
    vvMatrix(i,:)           = vv.';
    tPrev                   = tCurr;
end
rvInitialOrbit = rvMatrix;
vvInitialOrbit = vvMatrix;

% ----------------------------------------------------------------------------- %
%                             Terminal Circular Orbit                           %
% ----------------------------------------------------------------------------- %
ef               = 0;
omegaf        =  0*deg2rad;
nuf           =  0*deg2rad;
tauf          =  2*pi*sqrt(af^3/mu);
oef           = [af; ef; Omegaf; incf; omegaf; nuf];
[rvf0,vvf0]   = oe2rv_BraceyIV_Edgar(oef,mu);
hvf           = cross(rvf0,vvf0);
uhvf          = hvf/norm(hvf,2);
nvf           = cross(Ez,hvf);
tVector       = (0:tauf/5000:tauf).';
rMatrix       = zeros(length(tVector),3);
vMatrix       = zeros(length(tVector),3);
rv            = rvf0;
vv            = vvf0;
rvMatrix(1,:) = rvf0.';
vMatrix(1,:)  = vvf0.';
tPrev         = tVector(1);
for i = 2:length(tVector),
    tCurr                 = tVector(i);
    % [rv,vv,~,~,~,~] = propagateKepler_BraceyIV_Edgar(tPrev,tCurr,rv,vv,mu);
    [rv,vv,~,~,~,~] = propagateKepler_BraceyIV_Edgar(tVector(1),tCurr,rvf0,vvf0,mu);
    rvMatrix(i,:)          = rv.';
    vvMatrix(i,:)          = vv.';
    tPrev                  = tCurr;
end
rvTerminalOrbit = rvMatrix;
vvTerminalOrbit = vvMatrix;

% ----------------------------------------------------------------------------- %
%                          Characteristics of Transfer Orbit                    %
% ----------------------------------------------------------------------------- %
%   LOI        = Line of Intersection Between Initial and Terminal Orbits       %
%   r0Transfer = Initial Position on Transfer Orbit                             %
%   v0transfer = Initial Velocity on Transfer Orbit                             %
%   aTransfer  = Semi-Major Axis of Transfer Orbit                              %
% ----------------------------------------------------------------------------- %
if ~isequal(inc0,incf),
    LOI = cross(hv0,hvf)/norm(cross(hv0,hvf),2);
elseif isequal(inc0,0) || isequal(incf,0),
    LOI = [1; 0; 0];
else
    LOI = nv0/norm(nv0);        
end
rad0T      = a0;
speed0     = sqrt(mu/rad0T);
rv0T       = rad0T*LOI;
uv0        = cross(uhv0,rv0T)/norm(cross(uhv0,rv0T),2);  
vv0Tminus  = speed0*uv0;
aT         = (a0+af)/2;
if a0<af,
    eT = (af-a0)/(af+a0);
elseif a0>af,
    eT = (a0-af)/(a0+af);
else
    error('Semi-Major Axes Must Be Different Values');
end
    
DeltaV1mag = sqrt(mu)*sqrt(2/rad0T-1/aT)-speed0;
DeltaV1    = DeltaV1mag*uv0;
vv0T        = vv0Tminus + DeltaV1;

% ----------------------------------------------------------------------------- %
%            Elliptic Transfer Orbit Starting at Line of Intersection           %
%                       Between Initial and Terminal Orbits                     %
% ----------------------------------------------------------------------------- %
tauT       = 2*pi*sqrt(aT^3/mu);
tVector    = (0:(tauT/2)/500:tauT/2).';
rvT        = zeros(length(tVector),3);
vvT        = zeros(length(tVector),3);
rv         = rv0T;
vv         = vv0T;
rvCurr     = rv0T;
vvCurr     = vv0T;
rvT(1,:)   = rv0T;
vvT(1,:)   = vv0T;
for i = 2:length(tVector),
    [rvT(i,:),vvT(i,:),~,~,~,~] = propagateKepler_BraceyIV_Edgar(tVector(i-1),tVector(i),rvT(i-1,:).',vvT(i-1,:).',mu);
end

deltaInc    = incf - inc0;
rvfT        = rvT(end,:).';
vvfT        = vvT(end,:).';
speedfT     = norm(vvfT,2);
uvfT        = vvfT/speedfT.';
speedf      = sqrt(mu/af);
uvfT        = cross(uhvf,rvfT)/norm(rvfT,2);
vvf         = speedf*uvfT;
DeltaV2     = vvf - vvfT;

ctheta         = (vvf.'*vvfT)/(norm(vvf,2)*norm(vvfT,2));
theta          = acos(ctheta);
% ----------------------------------------------------------------------------- %
%                       Plot the Initial and Terminal Orbits                    %
% ----------------------------------------------------------------------------- %
rvInitialOrbit  = rvInitialOrbit/1000;
vvInitialOrbit  = vvInitialOrbit/1000;
rvTerminalOrbit = rvTerminalOrbit/1000;
vvTerminalOrbit = vvTerminalOrbit/1000;
pp = plot3(rvInitialOrbit(:,1),rvInitialOrbit(:,2),rvInitialOrbit(:,3),'.r');
set(pp,'LineWidth',2);
grid on
view((Omegaf-Omega0)*rad2deg,20);
hold on
pp = plot3(rvTerminalOrbit(:,1),rvTerminalOrbit(:,2),rvTerminalOrbit(:,3),'.b');
set(pp,'LineWidth',2);

% ----------------------------------------------------------------------------- %
%    Plot the Line of Intersection Between the Initial and Terminal Orbits      %
% ----------------------------------------------------------------------------- %
LOIPlot = [-rvfT rvfT].'/1000;
pp = plot3(LOIPlot(:,1),LOIPlot(:,2),LOIPlot(:,3),'-m');
set(pp,'LineWidth',2);

% ----------------------------------------------------------------------------- %
%                            Plot the Transfer Orbit                            %
% ----------------------------------------------------------------------------- %
rvT  = rvT/1000;
vvT  = vvT/1000;
vvf  = vvf/1000;
pp = plot3(rvT(:,1),rvT(:,2),rvT(:,3),'-k');
set(pp,'LineWidth',2);

% ----------------------------------------------------------------------------- %
%          Show the Impulses with Arrows Using the QUIVER3 Function             % 
% ----------------------------------------------------------------------------- %
vv0      = vv0T/1000;
rv0T     = rv0T/1000;
rvfT     = rvfT/1000;
vvfT     = vvfT/1000;
DeltaV1  = DeltaV1/1000;
DeltaV2  = DeltaV2/1000;
qq = quiver3(rv0T(1),rv0T(2),rv0T(3),DeltaV1(1),DeltaV1(2),DeltaV1(3),5000);
set(qq,'LineWidth',2);
qq = quiver3(rvfT(1),rvfT(2),rvfT(3),vvfT(1),vvfT(2),vvfT(3),5000);
set(qq,'LineWidth',2);
qq = quiver3(rvfT(1),rvfT(2),rvfT(3),vvf(1),vvf(2),vvf(3),5000);
set(qq,'LineWidth',2);
qq = quiver3(rvfT(1),rvfT(2),rvfT(3),DeltaV2(1),DeltaV2(2),DeltaV2(3),5000);
set(qq,'LineWidth',2);

ll = legend('Initial Orbit','Terminal Orbit','Line of Intersection','Transfer Orbit','$\Delta\mathbf{v}_1$','$\mathbf{v}(t_f)$ (Transfer Orbit)','$\mathbf{v}(t_f)$ (Terminal Orbit)','$\Delta\mathbf{v}_2$','Location','NorthEast');
set(ll,'Interpreter','LaTeX','FontSize',14);

% ----------------------------------------------------------------------------- %
%                         For Fun, Animate the Transfer Orbit                   %
% ----------------------------------------------------------------------------- %
set(ll,'AutoUpdate','off');
if animate,
    for i=1:size(rvT,1);
        pp = plot3(rvT(i,1),rvT(i,2),rvT(i,3),'-ko');
        pause(0.1);
    end;
end

xl = xlabel('$x$~(km$\times 1000$)');
yl = ylabel('$x$~(km$\times 1000$)');
zl = zlabel('$x$~(km$\times 1000$)');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(zl,'FontSize',16,'Interpreter','LaTeX');

hold on
qq = quiver3(0,0,0,hv0(1),hv0(2),hv0(3),0.0001);
set(qq,'LineWidth',2);
hold on
qq = quiver3(0,0,0,hvf(1),hvf(2),hvf(3),0.0001);
set(qq,'LineWidth',2);

theta = acos(dot(hv0,hvf)/(norm(hv0)*norm(hvf)));
fprintf('--------------------------------------------------------\n');
fprintf('----------Hohmann Transfer with a Plane Change----------\n');
fprintf('--------------------------------------------------------\n');
fprintf('--------------------------------------------------------\n');
fprintf('-----------------Data Supplied for This Run-------------\n');
fprintf('--------------------------------------------------------\n');
fprintf('Initial Semi-Major Axis \t\t= %8.4f (km)\n',a0);
fprintf('Terminal Semi-Major Axis \t\t= %8.4f (km)\n',af);
fprintf('Initial Inclination \t\t\t= %8.4f (deg)\n',inc0*180/pi);
fprintf('Terminal Inclination \t\t\t= %8.4f (deg)\n',incf*180/pi);
fprintf('Initial Longitude of Ascending Node \t= %8.4f (deg)\n',Omega0*180/pi);
fprintf('Terminal Longitude of Ascending Node \t= %8.4f (deg)\n',Omegaf*180/pi);
fprintf('--------------------------------------------------------\n');
fprintf('--------------------------------------------------------\n');
fprintf('---------------Features of Transfer Orbit---------------\n');
fprintf('--------------------------------------------------------\n');
fprintf('Semi-Major Axis \t\t\t= %8.4f (km)\n',aT);
fprintf('Eccentricity    \t\t\t= %8.4f (km)\n',eT);
fprintf('--------------------------------------------------------\n');
fprintf('--------------------------------------------------------\n');
fprintf('----------------Cranking Angle Required-----------------\n');
fprintf('--------------------------------------------------------\n');
fprintf('Cranking Angle \t\t\t\t= %8.4f (deg)\n',theta*180/pi);
fprintf('--------------------------------------------------------\n');
