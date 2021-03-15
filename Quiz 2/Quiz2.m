%Quiz 2
clc; clear;

%%Task 2 & 3
%Define given position and velocity vectors
rv0 = [-1217.39430415697 -3091.41210822807 -6173.40732877317]';
vv0 = [9.88635815507896 -0.446121737099303 -0.890884522967222]';
%Define Gauss points and known mu
mu = 398600;
N = 20;
%Calculate orbital elements with the given position and velocity
oe = rv2oe_BraceyIV_Edgar(rv0,vv0,mu)';
a = oe(1); %Define apoapsis
e = oe(2); %Define eccentricity
p = a*(1-e^2); %Define semi-latus rectum
tau = 2*pi*sqrt(a^3/mu); %Define the orbital period
%Define the integrand
f = @(nu) sqrt(p.^3/mu)./(1+e.*cos(nu)).^2;

%Calculate the value of t0
nu_ap = pi;
nu0 = oe(6);
[x, w] = GaussPointsWeights(nu_ap,nu0,N);
t0 = w'*f(x);
%Define the terminal time as five minutes after the initial time
t = t0+300;

%Define the number of intervals to run the root finder over
N_int = 2*round(tau)/round(t-t0); %Round to get integer values
%Preallocate the time and true anomaly storage
T = zeros(N_int+1,1); N = zeros(size(T));
%Define the first values that have already been calculated
T(1) = t0; N(1) = nu0;
%Preallocate the position and velocity storage 
R = zeros(length(T),3);  V = zeros(length(T),3); 
%Define the first values that have already been calculated
R(1,:) = rv0; V(1,:) = vv0;
%Run root finder for every five minute interval over two orbital periods
for ii = 2:length(T)
    %Calculate new nu value
    nu = rootFinder(f,t0,nu0,t,p,e,mu,N);
    %Redefine initial nu as the new nu
    nu0 = nu;
    %Store the new true anomaly and time values
    N(ii) = nu; T(ii) = T(ii-1) + 300;
    %Calculate the new position and velocity at the new nu
    oe(6) = nu;
    [rv, vv] = oe2rv_BraceyIV_Edgar(oe,mu);
    %Store the new position and velcoity values
    R(ii,:) = rv; V(ii,:) = vv;
end

%%Task 4
earthSphere %Plot the Earth
hold on
%Define the x, y, and z components of the position vectors
x = R(:,1); y = R(:,2); z = R(:,3);
plot3(x,y,z,'LineWidth',2) %Plot the orbit around the Earth
view(49.5,22.8) %View the orbit so that apoapsis and periapsis are visible

%%Task 5
lon_inertial = atan2(y,x); %Calculate the inertial longitude
lat = atan2(z,sqrt(x.^2+y.^2)); %Calculate the geocentric latitude
t = T-(t0+tau/2); %Calculate the time from apoapsis for each lon-lat
Earth_rot = 7.2921150e-05; %Rotation of the Earth [rad/s]
%Calculate Earth relative longitude
lonEarth = lon_inertial - Earth_rot*t; %Account for rotation of the Earth
lonEarth = mod(lonEarth,2*pi)-pi; %Put Earth relative longitude on [-pi,pi]
%Plot ground tracking
earth = imread('earth.jpg');
figure(2);
clf
image('CData',earth,'XData',[-180 180],'YData',[90 -90])
hold on
plot(lonEarth*180/pi,lat*180/pi,'*');

%%Task 6
%Part A
spyTime_n = []; spyTime_p = []; %Initialize storage for time over specified intervals
%Check each Earth relative longitude to see if it's in either specified range
for ll = 1:length(lonEarth)
    if lonEarth(ll) >= deg2rad(-100.5) && lonEarth(ll) <= deg2rad(-79.5) %[-100,-80]
        spyTime_n = [spyTime_n; T(ll)]; %Store each time where the lon is in range
    end
    if lonEarth(ll) >= deg2rad(79.5) && lonEarth(ll) <= deg2rad(100.5) %[80,100]
        spyTime_p = [spyTime_p; T(ll)]; %Store each time where the lon is in range
    end
end
%Calculate time spent in each range in hours
spyTime1 = (spyTime_p(length(spyTime_p))-spyTime_p(1))/3600;
spyTime2 = (spyTime_n(length(spyTime_n))-spyTime_n(1))/3600; 
fprintf('Task 6\n'); fprintf('(a)\n')
fprintf('The spacecraft spends %.3f hours on the Earth\n',spyTime1)
fprintf('relative longitude interval [80,100] degrees.\n')
fprintf('The spacecraft spends %.3f hours on the Earth\n',spyTime2)
fprintf('relative longitude interval [-100,-80] degrees.\n')
 
%Part B
fprintf('\n(b)\n')
fprintf('The segments of the orbit found in part (a)\n')
fprintf('are near the apoapsis since the time passed\n')
fprintf('on the interval is nearly all of the period.\n')

%Part C
%Find the max latitude and it's index
[maxlat, maxlatIndex] = max(lat);
%Find the one of the corresponding longitudes of the max latitude
lon_maxlat1 = lonEarth(maxlatIndex);
%Find the other longitude by adding half of the length of the
%latitude array to the longitude index. This works by taking 
%half of the array to a full orbital period.
lon_maxlat2 = lonEarth(maxlatIndex+(length(lat)-1)/2);
%Attach the lon-lats
max_loc1 = [lon_maxlat1 maxlat];
max_loc2 = [lon_maxlat2 maxlat];
%Convert the lon-lats to degrees for display
max_loc1_deg = rad2deg(max_loc1);
max_loc2_deg = rad2deg(max_loc2);
fprintf('\n(c)\n')
fprintf('The geocentric latitude reaches a maximum at the\n')
fprintf('coordinates (%.3f,%.3f) and (%.3f,%.3f) degrees.\n',max_loc1_deg,max_loc2_deg)

%Part D
fprintf('\n(d)\n')
fprintf('The two countries underneath the satellite at\n')
fprintf('apoapsis (in 1965) are Canada and the Soviet Union.\n')
fprintf('In Canada, the location is near the north side of\n')
fprintf('the Hudson Bay. In the Soviet Union, the location\n')
fprintf('is the small village of Turukhansk (%.3f,%.3f).\n',87.962222,65.793056)

%Part E
fprintf('\n(e)\n')
fprintf('The spacecraft would be spending the majority of it''s\n')
fprintf('orbit at the previously mentioned longitude ranges so that\n')
fprintf('it can observe North America for a long period of time\n')
fprintf('and then linger over the Soviet Union in order to transmit\n')
fprintf('the data that was collected. The high apoapsis lets the\n')
fprintf('spacecraft have a wide field-of-view at its height and the\n')
fprintf('shape allows the spacecraft to view the area of interest\n')
fprintf('over a wide range of angles. The direction of the eccentricity\n')
fprintf('vector and orientation of the orbit aligns the apoapsis\n')
fprintf('over the Northern hemisphere and the periapsis over the\n')
fprintf('Southern hemisphere. Since the majority of people live in\n')
fprintf('the Northern hemisphere, it is to the benefit of the Soviet\n')
fprintf('Union to observe that portion of the globe.\n')