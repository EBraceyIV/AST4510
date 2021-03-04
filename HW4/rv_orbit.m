function Y = rv_orbit(t,Z)
%Assign each value from Z to appropriate position/velocity component
x = Z(1); y = Z(2); z = Z(3); vx = Z(4); vy = Z(5); vz = Z(6);
%Calulate/define other values used in ODE
r = sqrt(x^2+y^2+z^2); mu = 398600; coef = -mu/r^3;
vxdot = coef*x; vydot = coef*y; vzdot = coef*z;
%Send values for ODE to be calculated with
Y = [vx vy vz vxdot vydot vzdot]';
end
