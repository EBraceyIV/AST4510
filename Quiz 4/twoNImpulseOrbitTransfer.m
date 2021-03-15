function twoNImpulseOrbitTransfer(r0,rf,inc0,incf,Omega0,Omegaf,N)
% INITAL & FINAL VALUES %
mu = 398600;
e0 = 0;
ef = 0;
omega0 = 0;
omegaf = 0;
nu00 = 0;
nuf0 = 0;
% nu = linspace(0,2*pi,200);
% nut = linspace(0,2*pi,200); %heheh nut
tau0 = 2*pi*sqrt(r0^3/mu);
tauf = 2*pi*sqrt(rf^3/mu);

% INITIAL ORBIT %
oe0 = [r0 e0 Omega0 inc0 omega0 nu00]';
[rv0, vv0] = oe2rv_BraceyIV_Edgar(oe0,mu);
time0 = linspace(0,tau0,500);
rv0_hold = [rv0.'; zeros(length(time0)-1,3)];
vv0_hold = [vv0.'; zeros(length(time0)-1,3)];
for ii = 2:length(time0)
    [rv, vv, ~, ~, ~, ~] = propagateKepler_BraceyIV_Edgar(time0(1),time0(ii),rv0,vv0,mu);
    [rv0_hold(ii,:)] = rv.';
    [vv0_hold(ii,:)] = vv.';
end

% FINAL ORBIT %
oef = [rf ef Omegaf incf omegaf nuf0]';
[rvf, vvf] = oe2rv_BraceyIV_Edgar(oef,mu);
timef = linspace(0,tauf,500);
rvf_hold = [rvf.'; zeros(length(timef)-1,3)];
vvf_hold = [vvf.'; zeros(length(timef)-1,3)];
for ii = 2:length(timef)
    [rv, vv, ~, ~, ~, ~] = propagateKepler_BraceyIV_Edgar(timef(1),timef(ii),rvf,vvf,mu);
    [rvf_hold(ii,:)] = rv.';
    [vvf_hold(ii,:)] = vv.';
end

hv0 = cross(rv0,vv0);
hvf = cross(rvf,vvf);
hv_hold = [hv0 0; 0; 0 0; 0; 0; hvf]

for n = 1:N
    rn = 0; Omegan = 0; incn = 0;
    omegan = 0; nun0 = 0;
    valsn = [rn Omegan incn]';
    vals0 = [r0 Omega0 inc0]';
    valsf = [rf Omegaf incf]';
    valsn = vals0 + n/N*(valsf-vals0);
    rn = valsn(1); Omegan = valsn(2); incn = valsn(3);
    
    an = (r0 + rn)/2;
    en = (rn - r0)/(rn + r0);
    taun = 2*pi*sqrt(an^3/mu);
    
    oen = [rn 0 Omegan incn omegan nun0]';
    [rvn, vvn] = oe2rv_BraceyIV_Edgar(oen,mu);
    timen = linspace(0,taun*2,500);
    rvn_hold = [rvn.'; zeros(length(timen)-1,3)];
    vvn_hold = [vvn.'; zeros(length(timen)-1,3)];
    for ii = 2:length(timen)
        [rv, vv, ~, ~, ~, ~] = propagateKepler_BraceyIV_Edgar(timen(1),timen(ii),rvn,vvn,mu);
        [rvn_hold(ii,:)] = rv.';
        [vvn_hold(ii,:)] = vv.';
    end
    plot3(rvn_hold(:,1),rvn_hold(:,2),rvn_hold(:,3))%,'LineWidth',2)
    hold on
    
    hv_n_1 = cross(rv0,vv0);
    hv_n = cross(rvf,vvf);
    lv = cross(hv0,hvf)/norm(cross(hv0,hvf));
    
    oet = [an en Omegan incn omegan nun0]';
    [rvn, vvn] = oe2rv_BraceyIV_Edgar(oet,mu);
    timet = linspace(0,taun/2);
    rvn_hold = [rvn.'; zeros(length(timet)-1,3)];
    vvn_hold = [vvn.'; zeros(length(timet)-1,3)];
    for ii = 2:length(timet)
        [rv, vv, ~, ~, ~, ~] = propagateKepler_BraceyIV_Edgar(timet(1),timet(ii),rvn,vvn,mu);
        [rvn_hold(ii,:)] = rv.';
        [vvn_hold(ii,:)] = vv.';
    end
    plot3(rvn_hold(:,1),rvn_hold(:,2),rvn_hold(:,3),':','LineWidth',2)
    hold on
    
    rimp = rn/lv;
    
    
    
end



% PLOT ORBITS %
plot3(rv0_hold(:,1),rv0_hold(:,2),rv0_hold(:,3),'g','LineWidth',2)
hold on
plot3(rvf_hold(:,1),rvf_hold(:,2),rvf_hold(:,3),'r','LineWidth',2)