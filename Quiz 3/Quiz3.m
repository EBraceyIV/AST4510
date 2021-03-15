%Quiz 3
clc; clear;
%%% TASK 1
% % Initializing Variables
%Define orbital period (converted to seconds)
tau_hr = 24; tau = tau_hr*3600;
%Define gravitational parameter
mu = 398600; %km^3/s^-2
%Define semi-major axis
a = ((tau^2*mu)/(4*pi^2))^(1/3); %km
%Define eccentricity
e = 0.25;
%Define orbital inclination (converted to radians)
inc = deg2rad(63.4);
%Define values of argument of periapsis (converted to radians)
omega = deg2rad([90 270]');
%Define values of longitude of the ascending node (converted to radians)
Omega = deg2rad([0 45 90 135 180]');
%Define initial nu
nu0 = 0; %radians
%Define timespan (converted to seconds)
tstep_m= 5; tstep = tstep_m*60; %Time step of five minutes (300 seconds)
t0 = 0; tf = tau; %Initial and final times
tspan = t0:tstep:tf; %Time span of whole orbit (seconds)
%Initialize storage for all positions
%The 3 in the column size accounts for the 3 components of position
rPCI_all = zeros(length(tspan),3*length(Omega),length(omega));
% % Plot orbits
for ii = 1:length(omega)
    figure %Open new figure for new omega value
    earthSphere %Plot the Earth
    hold on
    %Title figure appropriately
    if ii == 1; title('Orbital Paths, \omega = 90 degrees');
    else; title('Orbital Paths, \omega = 270 degrees');
    end
    %Calculate initial conditions for a specific Omega value
    for jj = 1:length(Omega)
        %Define initial orbital elements array
        oe = [a e Omega(jj) inc omega(ii) nu0];
        %Calculate initial position & velocity
        [rPCI0, vPCI0] = oe2rv_BraceyIV_Edgar(oe,mu);
        %Initialize position & velocity storage
        rPCI = [rPCI0'; zeros(length(tspan)-1,3)];
        %Use Kepler solver to find position & velocity every five minutes
        for kk = 2:length(tspan)
            [rPCIf,vPCIf,E0,nu0,E,nu] = propagateKepler_BraceyIV_Edgar(t0,tspan(kk),rPCI0,vPCI0,mu);
            %Store new position & velocity
            rPCI(kk,:) = rPCIf';
        end
        rPCI_all(:,((jj*3)-2):jj*3,ii) = rPCI;
        %Plot the orbit for the current Omega value
        plot3(rPCI(:,1),rPCI(:,2),rPCI(:,3))
    end
    %Label the orbits for each Omega value
    legend('Earth','\Omega = 0 degrees','\Omega = 45 degrees','\Omega = 90 degrees','\Omega = 135 degrees','\Omega = 180 degrees');    
end
%Print comments regarding orbits
fprintf('Task 1\n')
fprintf('------\n')
fprintf('The orbits shift around Earth due to increasing longitude of the\n')
fprintf('ascending node values such that the line of nodes moves farther\n')
fprintf('from the inertially fixed x-direction of the Earth but maintain\n')
fprintf('the same angle between the angular momentum vector and the\n')
fprintf('inertially fixed z-direction of the Earth.\n')
fprintf('Additionally, the change in the argument of the periapsis causes\n')
fprintf('a "reflection" to take place across eccentricity vector since the\n')
fprintf('angle between the eccentricity vector and the line of nodes is\n')
fprintf('increased by 180 degrees.\n')
%%% TASK 2 & 3
% % Groundtracking & Orbital Analysis
%Initialize longitude & latitude storage
lonEarth = zeros(length(tspan),length(Omega),length(omega));
lat = zeros(size(lonEarth));
%Define groundtracking markers
mark = ["o" "s" "d" "v" ">"];
%Print comments regarding groundtracking
fprintf('\nTask 2\n')
fprintf('------\n')
fprintf('When the argument of the periapsis is equal to 90 degrees\n')
fprintf('the spacecraft spends most of its time over the Southern\n')
fprintf('hemisphere. When the argument of the periapsis is equal\n')
fprintf('to 270 degrees the spacecraft spends most of its time\n')
fprintf('over the Northern hemisphere. In both of these instances\n')
fprintf('the spacecraft covers a specific region of longitude on\n')
fprintf('Earth. This may be to collect/transmit data between the\n')
fprintf('regions in different hemispheres with the same longitudes.\n')
fprintf('The difference between the two values of the argument of\n')
fprintf('the periapsis results in the groundtracks being reflections\n')
fprintf('of each other across the equator. They are also shifted\n')
fprintf('by 180 degrees past the Prime Meridian from each other.\n')
fprintf('\nTask 3\n')
fprintf('------\n')
fprintf('(a)\n')
for ll = 1:length(omega)
    %Plot Earth
    earth = imread('earth.jpg');
    figure; clf
    image('CData',earth,'XData',[-180 180],'YData',[90 -90])
    hold on
    %Attach appropriate title
    if ll == 1; title('Orbital Groundtracking, omega = 90 degrees');
    else; title('Orbital Groundtracking, omega = 270 degrees');
    end
    %%Calculations per Omega value
    for mm = 1:length(Omega)
        %%Groundtracking
        x = rPCI_all(:,(mm*3)-2,ll); y = rPCI_all(:,(mm*3)-1,ll); z = rPCI_all(:,mm*3,ll);
        lon_inertial = atan2(y,x); %Calculate the inertial longitude
        lat(:,mm,ll) = atan2(z,sqrt(x.^2+y.^2)); %Calculate the geocentric latitude
        t = tspan-(t0+tau/2); %Calculate the time from apoapsis for each lon-lat
        Earth_rot = 7.2921150e-05; %Rotation of the Earth [rad/s]
        %Calculate Earth relative longitude
        lonEarth(:,mm,ll) = lon_inertial - Earth_rot*t'; %Account for rotation of the Earth
        lonEarth(:,mm,ll) = mod(lonEarth(:,mm,ll),2*pi)-pi; %Put Earth relative longitude on [-pi,pi]
        plot(lonEarth(:,mm,ll)*180/pi,lat(:,mm,ll)*180/pi,mark(mm));
        hold on
        %%Orbital Analysis
        if ll == 1 %Looking at omega = 90, (b)
            for ww = 1:length(tspan) %Check lat-lon for each time interval
                %Check if spacecraft passes over Australia or South America
                %Rough bounds for Australia: -40 to -15 degrees latitude
                %                            115 to 155 degrees longitude
                %Rough bounds for S. America: 10 to -60 degrees latitude
                %                            -80 to -40 degrees longitude
                if rad2deg(lat(ww,mm,ll)) <= -15 && rad2deg(lat(ww,mm,ll)) >= -40 && rad2deg(lonEarth(ww,mm,ll)) >= 115 && rad2deg(lonEarth(ww,mm,ll)) <= 155
                    OmegaAust = Omega(mm); break
                elseif rad2deg(lat(ww,mm,ll)) <= 10 && rad2deg(lat(ww,mm,ll)) >= -60 && rad2deg(lonEarth(ww,mm,ll)) >= -80 && rad2deg(lonEarth(ww,mm,ll)) <= -40
                    OmegaSA = Omega(mm); break
                end
            end
        end
        if ll == 2 %Looking at omega = 270, (c)
            for ww = 1:length(tspan) %Check lat-lon for each time interval
                %Check if spacecraft passes over central Russia or central North America
                %Rough bounds for Russia: 50 to 70 degrees latitude
                %                         60 to 120 degrees longitude
                %Rough bounds for N. America: 35 to 55 degrees latitude
                %                           -110 to -70 degrees longitude
                if rad2deg(lat(ww,mm,ll)) <= 70 && rad2deg(lat(ww,mm,ll)) >= 50 && rad2deg(lonEarth(ww,mm,ll)) >= 60 && rad2deg(lonEarth(ww,mm,ll)) <= 120
                    OmegaRuss = Omega(mm); break
                elseif rad2deg(lat(ww,mm,ll)) <= 55 && rad2deg(lat(ww,mm,ll)) >= 35 && rad2deg(lonEarth(ww,mm,ll)) >= -110 && rad2deg(lonEarth(ww,mm,ll)) <= -70
                    OmegaNA = Omega(mm); break
                end
            end
        end
    end
    %Label the plotted orbits for each Omega value
    legend('\Omega = 0 degrees','\Omega = 45 degrees','\Omega = 90 degrees','\Omega = 135 degrees','\Omega = 180 degrees');
    %%Find which hemisphere the spacecraft stays over the longest
    %Initialize time spent in either hemisphere
    North_time = 0; South_time = 0;
    %Add up the time spent in either hemisphere per latitude value
    for uu = 1:length(tspan)
        if lat(uu,1,ll) > 0; North_time = North_time + tstep;
        else; South_time = South_time + tstep;
        end
    end
    %Label and output the print statements
    part = ["(i)" "(ii)"]; fprintf('%s\n',part(ll))
    if North_time > South_time
        fprintf('The spacecraft spends the majority of\n')
        fprintf('its time over the Northern hemisphere.\n')
    else
        fprintf('The spacecraft spends the majority of\n')
        fprintf('its time over the Southern hemisphere.\n')
    end
    %%Find how long it takes the spacecraft return to the crossing point
    %Sorted the longitudes in ascending order
    [lonEarthSorted, indicesSorted] = sort(lonEarth(:,1,ll));
    %Sort the values of time according to longitude
    timeSorted = tspan(indicesSorted);
    %Find the difference between the time values
    dTimeSorted = diff(timeSorted);
    %The absolute value of the maximum difference in time is how long it
    %takes for the spacecraft
    crossTime = max(abs(dTimeSorted));
    fprintf('The spacecraft takes %g seconds (or %g hours)\n',crossTime,crossTime/3600)
    fprintf('to return to the crossing point once it is passed.\n')
end
fprintf('\n(b)\n')
fprintf('For the value omega = 90 degrees the spacecraft passes\n')
fprintf('over Australia when Omega = %g degrees and passes over\n',rad2deg(OmegaAust))
fprintf('South America when Omega = %g degrees.\n',rad2deg(OmegaSA))
fprintf('\n(c)\n')
fprintf('For the value omega = 270 degrees the spacecraft passes\n')
fprintf('over Russia when Omega = %g degrees and passes over\n',rad2deg(OmegaRuss))
fprintf('North America when Omega = %g degrees.\n',rad2deg(OmegaNA))