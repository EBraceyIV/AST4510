%Homework 5 Problem 1 & 2
clc; clear;
%Define normalized Hohmann transfer impulses
H_Dv1_nodim = @(R) sqrt((2.*R)./(1+R))-1;
H_Dv2_nodim = @(R) sqrt(1./R).*sqrt(1-sqrt(2./(1+R)));
H_Dv_nodim = @(R) (H_Dv1_nodim(R) + H_Dv2_nodim(R));
%Define normalized bi-ellipctic transfer impulses
BE_Dv1_nodim = @(R,S) sqrt((2.*R)./(1./S+R))-1;
BE_Dv2_nodim = @(R,S) sqrt(1./(R.*S)).*(sqrt(2./(1+S))-sqrt(2./(1+R.*S)));
BE_Dv3_nodim = @(R,S) sqrt(2.*S./(R+R.*S))-sqrt(1./R);
BE_Dv_nodim = @(R,S) (BE_Dv1_nodim(R,S) + BE_Dv2_nodim(R,S) + BE_Dv3_nodim(R,S));
%Define normalized bi-parabolic transfer impulses
BP_Dv1_nodim = @(R) sqrt(2)-1;
BP_Dv2_nodim = @(R) 0;
BP_Dv3_nodim = @(R) sqrt(1./R).*(sqrt(2)-1);
BP_Dv_nodim = @(R) (BP_Dv1_nodim(R) + BP_Dv2_nodim(R) + BP_Dv3_nodim(R));
%Define varying orbit parameters
R = linspace(1,20,30);
S = [2 5 10 11 12 15]';
%Color and denote different curves on plot
colorMarks = ["-ob" "-sg" "-dr" "-vy" "->m" "-<k"];

%Plot Hohmann transfer impulses
plot(R,H_Dv_nodim(R),'k-*'); hold on
%Plot bi-elliptic transfer impulses
for ii = 1:length(S)
    plot(R,BE_Dv_nodim(R,S(ii)),colorMarks(ii)); hold on
end
%Plot bi-parabolic transfer impulses
plot(R,BP_Dv_nodim(R),"-cp");
%Label plot
xlabel('$$R$$','interpreter','latex')
ylabel('$$\frac{\Delta v}{v_{c1}}$$','interpreter','latex')
title('Normalized impulses as functions of R')
legend('Hohmann Transfer','Bi-Elliptic Transfer (S=2)','Bi-Elliptic Transfer (S=5)','Bi-Elliptic Transfer (S=10)','Bi-Elliptic Transfer (S=11)','Bi-Elliptic Transfer (S=12)','Bi-Elliptic Transfer (S=15)','Bi-Parabolic Transfer','Location','best')
% %Zoomed in view
figure
%Plot Hohmann transfer impulses
plot(R,H_Dv_nodim(R),'k-*'); hold on
%Plot bi-elliptic transfer impulses
for ii = 1:length(S)
    plot(R,BE_Dv_nodim(R,S(ii)),colorMarks(ii)); hold on
end
%Plot bi-parabolic transfer impulses
plot(R,BP_Dv_nodim(R),"-cp");
%Label plot
xlabel('$$R$$','interpreter','latex')
ylabel('$$\frac{\Delta v}{v_{c1}}$$','interpreter','latex')
title('Normalized impulses as functions of R')
legend('Hohmann Transfer','Bi-Elliptic Transfer (S=2)','Bi-Elliptic Transfer (S=5)','Bi-Elliptic Transfer (S=10)','Bi-Elliptic Transfer (S=11)','Bi-Elliptic Transfer (S=12)','Bi-Elliptic Transfer (S=15)','Bi-Parabolic Transfer','Location','best')
%Zoom in on area of interest
xlim([10 16]); ylim([0.51 0.55]);

% %%% Problem 2
% HohmannBi1 = intersect(H_Dv_nodim(R),BE_Dv_nodim(R,S(4)))
% chunk = intersect([1 4 2 8],[3 7 5 4]);
% 
% 

