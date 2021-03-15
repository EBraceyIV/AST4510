%Homework 5 Problem 1 & 2
%%% Problem 1
clc; clear;
%Normalized Hohmann transfer impulses
H_Dv1 = @(R) sqrt((2.*R)./(1+R))-1;
H_Dv2 = @(R) sqrt(1./R).*(1-sqrt(2./(1+R)));
H_Dv = @(R) H_Dv1(R) + H_Dv2(R);
%Normalized bi-elliptic transfer impulses
BE_Dv1 = @(R,S) sqrt((2*R)./(1./S+R))-1;
BE_Dv2 = @(R,S) sqrt(1./(R.*S)).*(sqrt(2./(1+S))-sqrt(2./(1+R.*S)));
BE_Dv3 = @(R,S) sqrt(2.*S./(R+R.*S))-sqrt(1./R);
BE_Dv = @(R,S) BE_Dv1(R,S) + BE_Dv2(R,S) + BE_Dv3(R,S);
%Normalized bi-parabolic transfer impulses
BP_Dv1 = @(R) sqrt(2)-1;
BP_Dv2 = @(R) 0;
BP_Dv3 = @(R) sqrt(1./R).*(sqrt(2)-1);
BP_Dv = @(R) BP_Dv1(R) + BP_Dv2(R) + BP_Dv3(R);
%Define varying orbit paramters
R = linspace(1,20,30);
S = [2 5 10 11 12 15]';
%Color and mark different curves on plot
colorMarks = ["-ob" "-sg" "-dr" "-vy" "->m" "-<k"];
%Plot Hohman transfer impulses
plot(R,H_Dv(R),"k-*"); hold on
%Plot bi-ellipcti trasnfer impulses
for ii = 1:length(S)
    plot(R,BE_Dv(R,S(ii)),colorMarks(ii)); hold on
end
%Plot bi-parabolic transfer impulses
plot(R,BP_Dv(R),"-cp")
%Label plot
xlabel('$$R$$','interpreter','latex')
ylabel('$$\frac{\Delta v}{v_{c1}}$$','interpreter','latex')
title('Normalized impulses as functions of R')
legend('Hohmann Transfer','Bi-Elliptic Transfer (S=2)','Bi-Elliptic Transfer (S=5)','Bi-Elliptic Transfer (S=10)','Bi-Elliptic Transfer (S=11)','Bi-Elliptic Transfer (S=12)','Bi-Elliptic Transfer (S=15)','Bi-Parabolic Transfer','Location','best')
% %Zoomed in view
figure
%Plot Hohman transfer impulses
plot(R,H_Dv(R),"k-*"); hold on
%Plot bi-ellipcti trasnfer impulses
for ii = 1:length(S)
    plot(R,BE_Dv(R,S(ii)),colorMarks(ii)); hold on
end
%Plot bi-parabolic transfer impulses
plot(R,BP_Dv(R),"-cp")
%Label plot
xlabel('$$R$$','interpreter','latex')
ylabel('$$\frac{\Delta v}{v_{c1}}$$','interpreter','latex')
title('Normalized impulses as functions of R (zoomed)')
legend('Hohmann Transfer','Bi-Elliptic Transfer (S=2)','Bi-Elliptic Transfer (S=5)','Bi-Elliptic Transfer (S=10)','Bi-Elliptic Transfer (S=11)','Bi-Elliptic Transfer (S=12)','Bi-Elliptic Transfer (S=15)','Bi-Parabolic Transfer','Location','best')
%Zoom in on area of interest
xlim([10 16]); ylim([0.51 0.55]);

%%% Problem 2
%Find the point where bi-parabolic becomes more efficient than Hohmann
Hoh_BP_root = @(R) H_Dv(R) - BP_Dv(R);
Hoh_BP_cross = fsolve(Hoh_BP_root,R(1));
%Find the point where bi-elliptic becomes more efficient than Hohmann
%for varying values of S
Hoh_BE_cross = zeros(size(S));
for jj = 1:length(S)
    Hoh_BE_root = @(R) H_Dv(R) - BE_Dv(R,S(jj));
    Hoh_BE_cross(jj) = fsolve(Hoh_BE_root,R(1));
end
figure
plot(Hoh_BE_cross,S,'-*')
xlabel('$$S$$','interpreter','latex')
ylabel('$$R$$','interpreter','latex')
title('R for S at which BE is more efficient than Hoh')
legend('R vs. S')