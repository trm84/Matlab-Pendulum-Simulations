% Tyler Matthews
% Pendulum Simulation / Animation
% ** NEED rk2Calculator in the same directory as this file
clc; close all; clear all;

disp('Double Pendulum Simulation - Tyler Matthews');

%% Changable Parameters
    % Simulation Parameters
        stopTime = 50;   % how long to run the simulation (seconds)
        steadyCam = 1;   % if(1) camera doesn't move, if(steadyCam!=1) camera follows pendulum
        traceLength = 10; % how long should the pendulum be traced (seconds)        

    % Model Parameters
        initialTheta1 = pi; % Starting Position of mass 1 (radians)
        initialTheta2 = pi;    % Starting Position of mass 2 (radians)
        initialOmega1 = 0.1;    % Starting Velocity of mass 1 (radians/second)
        initialOmega2 = 0;    % Starting Velocity of mass 2 (radians/second)
        m1 = 10;             % mass 1 (kg)
        m2 = 10;             % mass 2 (kg)
        l1 = 7;             % pendulum rod 1 length (meters)
        l2 = 3;             % pendulum rod 2 length (meters)
        g = 9.8;            % gravitational constant (m/ms^2)
 
%% Initializing
plotTitle = sprintf('Time = %i seconds, m1 = %i, m2 = %i, l1 = %i, l2 = %i, g = %0.2f', stopTime, m1, m2, l1, l2, g);

% Simulation Parameters
    startTime = 0;
    T = 0.05;                               % each step is 50ms
    steps = stopTime/T;                     % # of steps
    t = linspace(startTime,stopTime,steps); % time vector

%Initialize thetas and omegas
    theta1 = zeros(1, steps);     
    theta2 = zeros(1, steps);
    omega1 = zeros(1,steps);
    omega2 = zeros(1, steps);     

    theta1(1) = initialTheta1;    
    omega1(1) = initialOmega1;    
    theta2(1) = initialTheta2;
    omega2(1) = initialOmega2;

%Inital plotting
    xStartPoint = (l1+l2);    % Sets the x anchored point of the pendulum
    yStartPoint = (l1+l2);    % Sets the y anchored point of the pendulum

    x1_pos = l1*sin(theta1(1));
    y1_pos = l1*cos(theta1(1)); 
    x2_pos = x1_pos + l2*sin(theta2(1)); 
    y2_pos = y1_pos + l2*cos(theta2(1));

%Trace the pendulum
    % posArr1 = [x values of rod 1 ; y values of rod 1]
        posArr1 = [zeros(1, length(t)); zeros(1, length(t))];
        posArr1(1,1) = x1_pos;
        posArr1(1,2) = y1_pos;
    % posArr2 = [x values of rod 2 ; y values of rod 2]
        posArr2 = [zeros(1, length(t)); zeros(1, length(t))];
        posArr2(1,1) = x2_pos;
        posArr2(1,2) = y2_pos;
    index = 1;    

%Variables for drawing circles / masses
    th = 0:pi/50:2*pi;          
    circleSize1 = (1/20)*m1;       % Mass 1 size
    circleSize2 = (1/20)*m2;       % Mass 2 size
    stationaryCircleSize = 0.25;   % Anchored circle size
    
    %{
        Double Pendulum Equations Available Here: 
        http://web.mit.edu/jorloff/www/chaosTalk/double-pendulum/double-pendulum-en.html
    %}

%% Simulation
figure;
%double pedulum with damping
for k=2:steps
   %Runge Kutta Method  
       temp = rk2Calculator(theta1, theta2, omega1, omega2, m1, m2, l1, l2, g, k, T);
       theta1(k) = temp(1);
       theta2(k) = temp(2);
       omega1(k) = temp(3);
       omega2(k) = temp(4);
      
   % Current Position
       x1_pos = l1*sin(theta1(k));
       y1_pos = l1*cos(theta1(k));
       x2_pos = l2*sin(theta2(k)); 
       y2_pos = l2*cos(theta2(k));

   % Adding to tracing array
       posArr1(1, index) = xStartPoint - x1_pos;
       posArr1(2, index) = yStartPoint - y1_pos;
       posArr2(1, index) = posArr1(1, index) - x2_pos;
       posArr2(2, index) = posArr1(2, index) - y2_pos;
   
   if(index >= length(posArr1))
        lastIndex = index;
        index = 1;
   else
        lastIndex = index;
        index = index + 1;
   end
   
   % PLOTTING   
       clf;
       hold on
           % Plotting mass 1 and mass 2 traces
               if(k*T <= traceLength)
                   plot(posArr1(1,1:k-1), posArr1(2,1:k-1), 'b')
                   plot(posArr2(1,1:k-1), posArr2(2,1:k-1), 'r')
               else
                   plot(posArr1(1,k-(traceLength/T):k-1), posArr1(2,k-(traceLength/T):k-1), 'b')
                   plot(posArr2(1,k-(traceLength/T):k-1), posArr2(2,k-(traceLength/T):k-1), 'r')
               end
           % Plotting rod 1 and rod 2
               plot([xStartPoint, posArr1(1,lastIndex)],[yStartPoint, posArr1(2,lastIndex)], 'black');
               plot([posArr1(1,lastIndex), posArr2(1,lastIndex)],[posArr1(2,lastIndex), posArr2(2,lastIndex)], 'black'); 

          % Draw anchored point 
               circleX = stationaryCircleSize * cos(th) + xStartPoint;
               circleY = stationaryCircleSize * sin(th) + yStartPoint;
               plot(circleX, circleY, 'black');
               fill(circleX, circleY, 'black');
          % Draw the mass on the end of the pendulum rod 1
               circleX = circleSize1*cos(th) + posArr1(1, lastIndex);
               circleY = circleSize1*sin(th) + posArr1(2,lastIndex);
               plot(circleX, circleY, 'b');
               fill(circleX, circleY, 'g');
          % Draw the mass on the end of the pendulum rod 2
               circleX = circleSize2*cos(th) + posArr2(1,lastIndex);
               circleY = circleSize2*sin(th) + posArr2(2,lastIndex);
               plot(circleX, circleY, 'b');
               fill(circleX, circleY, 'g');
       hold off
       title(plotTitle)
       
    % Adjust camera View
       if(steadyCam  == 1)
           xlim([-1,2*(l1+l2) + 1]);
           ylim([-1,2*(l1+l2) + 1]);
       else
           xlim([posArr2(1, lastIndex) - (l1 + l2) - 1, posArr2(1, lastIndex) + (l1 + l2) + 1]);
           ylim([posArr2(2, lastIndex) - (l1 + l2) - 1, posArr2(2, lastIndex) + (l1 + l2) + 1]);
       end   
     
    % ANIMATE
       if(mod(k,1) == 0)
            pause(0.01) 
            %disp(k)
       end
end
