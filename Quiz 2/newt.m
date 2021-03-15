function x = newt(fx,dfx,x1,N)
%Hold on to the initial value & preallocate
x = [x1 zeros(1,N)];
for ii = 1:N %Run n iterations using Newton's method
    x = x1 - fx(x1)/dfx(x1);
    x1 = x2;
    x(ii+1) = x1; %Hold on to each iteration of x
end
