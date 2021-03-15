%Quiz 4
clc; clear;

%Earth values
E_r = 6378.145; %km
mu = 386000; %km^3*s^-2

%Inital orbit
alt0 = 500; %km
r0 = E_r + alt0; %km
inc0 = deg2rad(28.5); %radians
Omega0 = deg2rad(60); %radians

%Final orbit
altf = 20200; %km
rf = E_r + altf; %km
incf = deg2rad(57); %radians
Omegaf = deg2rad(120); %radians

N = 3;

twoNImpulseOrbitTransfer(r0,rf,inc0,incf,Omega0,Omegaf,N)