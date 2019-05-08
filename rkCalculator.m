function output = rkCalculator(theta, omega, T, g, l, k, b, m)
%This function calculates theta and omega using the runge kutta method

x11 = T *omega(k-1);
x21 = T *((-b/m)*omega(k-1) + (-g/l)*(sin(theta(k-1)))); 

x12 = T *(omega(k-1) + 0.5*x21);
x22 = T *((-b/m)*(omega(k-1) + 0.5*x21) + (-g/l)*(sin(theta(k-1) + 0.5*x11)));

x13 = T *(omega(k-1) + 0.5*x22);
x23 = T *((-b/m)*(omega(k-1) + 0.5*x22) + (-g/l)*(sin(theta(k-1) + 0.5*x12)));

x14 = T *(omega(k-1) + x23);
x24 = T *((-b/m)*(omega(k-1) + 0.5*x23) + (-g/l)*(sin(theta(k-1) + x13)));


output = zeros(1,2);
output(1) = theta(k-1) + (x11 + 2*x12 + 2*x13 + x14)/6;
output(2) = omega(k-1) + (x21 + 2*x22 + 2*x23 + x24)/6;