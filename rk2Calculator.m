function output = rk2Calculator(theta1, theta2, omega1, omega2, m1, m2, l1, l2, g, k, T)
    %This function calculates theta1/2 and omega1/2 using the runge kutta method

    %Runge Kutta Method
    x11 = theta1Calc(omega1, 0, k, T, 0);
    x21 = theta2Calc(omega2, 0, k, T, 0);
    x31 = omega1Calc(theta1, theta2, omega1, omega2, 0, 0, 0, 0, m1, m2, l1, l2, g, k, T, 0);
    x41 = omega2Calc(theta1, theta2, omega1, omega2, 0, 0, 0, 0, m1, m2, l1, l2, g, k, T, 0);
    
    x12 = theta1Calc(omega1, x31, k, T, 0.5);
    x22 = theta2Calc(omega2, x41, k, T, 0.5);
    x32 = omega1Calc(theta1, theta2, omega1, omega2, x11, x21, x31, x41, m1, m2, l1, l2, g, k, T, 0.5);
    x42 = omega2Calc(theta1, theta2, omega1, omega2, x11, x21, x31, x41, m1, m2, l1, l2, g, k, T, 0.5);
    
    x13 = theta1Calc(omega1, x32, k, T, 0.5);
    x23 = theta2Calc(omega2, x42, k, T, 0.5);
    x33 = omega1Calc(theta1, theta2, omega1, omega2, x12, x22, x32, x42, m1, m2, l1, l2, g, k, T, 0.5);
    x43 = omega2Calc(theta1, theta2, omega1, omega2, x12, x22, x32, x42, m1, m2, l1, l2, g, k, T, 0.5);
    
    x14 = theta1Calc(omega1, x33, k, T, 1);
    x24 = theta2Calc(omega2, x43, k, T, 1);
    x34 = omega1Calc(theta1, theta2, omega1, omega2, x13, x23, x33, x43, m1, m2, l1, l2, g, k, T, 1);
    x44 = omega2Calc(theta1, theta2, omega1, omega2, x13, x23, x33, x43, m1, m2, l1, l2, g, k, T, 1);
    
    output = zeros(1,4);
    
    output(1) = theta1(k-1) + (x11 + 2*x12 + 2*x13 + x14)/6;
    output(2) = theta2(k-1) + (x21 + 2*x22 + 2*x23 + x24)/6;
    output(3) = omega1(k-1) + (x31 + 2*x32 + 2*x33 + x34)/6;
    output(4) = omega2(k-1) + (x41 + 2*x42 + 2*x43 + x44)/6;
  
end

% x1 = theta1, x2 = theta2, x3 = omega1, x4 = omega2
function newTheta1 = theta1Calc(omega1, x3, k, T, multiplier)
    newTheta1 = T*(omega1(k-1) + (multiplier*x3));
end

function newTheta2 = theta2Calc(omega2, x4, k, T, multiplier)
    newTheta2 = T*(omega2(k-1) + (multiplier*x4));
end

function newOmega1 = omega1Calc(theta1, theta2, omega1, omega2, x1, x2, x3, x4, m1, m2, l1, l2, g, k, T, multiplier)
    %num
    part1 = -g*(2*m1 + m2)*sin(theta1(k-1) + multiplier*x1);
    part2 = -m2*g*sin((theta1(k-1) + multiplier*x1) - 2*(theta2(k-1) + multiplier*x2));
    part3 = -2*sin((theta1(k-1) + multiplier*x1) - (theta2(k-1) + multiplier*x2))*m2;
    part4 = ((omega2(k-1) + multiplier*x4)^2)*l2 + ((omega1(k-1) + multiplier*x3)^2)*l1*cos((theta1(k-1) + multiplier*x1) - (theta2(k-1) + multiplier*x2));

    %den
    part5 = l1*(2*m1 + m2 - m2*cos(2*(theta1(k-1) + multiplier*x1) - 2*(theta2(k-1) + multiplier*x2)));
    
    newOmega1 = T*((part1 + part2 + (part3*part4))/part5);
end
    
function newOmega2 = omega2Calc(theta1, theta2, omega1, omega2, x1, x2, x3, x4, m1, m2, l1, l2, g, k, T, multiplier)
    %num
    part1 = 2*sin((theta1(k-1) + multiplier*x1) - (theta2(k-1) + multiplier*x2));
    part2 = ((omega1(k-1) + multiplier*x3)^2)*l1*(m1+m2);
    part3 = g*(m1 + m2)*cos(theta1(k-1) + multiplier*x1);
    part4 = ((omega2(k-1) + multiplier*x4)^2)*l2*m2*cos((theta1(k-1) + multiplier*x1) - (theta2(k-1) + multiplier*x2));
   
    %den
    part5 = l2*(2*m1 + m2 - m2*(cos(2*(theta1(k-1) + multiplier*x1) - 2*(theta2(k-1) + multiplier*x2))));
    
    newOmega2 = T*((part1*(part2 + part3 + part4))/part5);
end