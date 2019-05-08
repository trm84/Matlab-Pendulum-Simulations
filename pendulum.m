% Tyler Matthews
% Single Pendulum Simulation / Animation
% 5/7/2019
% ** NEED rkCalculator in the same directory as this file
clc; close all; clear all;

disp('Single Pendulum Simulation - Tyler Matthews');

%% Changable Parameters
    % Simulation Parameters
        stopTime = 50;   % how long to run the simulation (seconds)
        steadyCam = 1;   % if(1) camera doesn't move, if(steadyCam!=1) camera follows pendulum
        traceLength = 3; % how long should the pendulum be traced (seconds)

    %Model Parameters
        initialTheta = 3*pi/4;    % starting position (radians)
        initialOmega = 0;       % starting velocity (radians / second)
        b = 0.3;                % damping factor
        m = 10;                 % mass (kg)
        g = 9.8;                % gravitational constant (m/s^2)
        l = 7;                  % pendulum length (meters)




%% Initializing
plotTitle = sprintf('Time = %i seconds, m = %i, l = %i, g = %0.2f', stopTime, m, l, g);

startTime = 0;                          % start simulation at 0 seconds
T = 0.05;                               % time step = 50ms
steps = stopTime/T;                     % total number of steps in simulation
t = linspace(startTime,stopTime,steps); % time vector

%Initialize
    theta = zeros(1, steps); %theta 
    omega = zeros(1, steps); %omega = dtheta/dt

    theta(1) = initialTheta;
    omega(1) = initialOmega;

    xStartPoint = l;        % Sets the x anchored point of the pendulum
    yStartPoint = l;    % Sets the y anchored point of the pendulum

%inital plot 
    x_pos = l*sin(theta(1));
    y_pos = l*cos(theta(1)); 

%Array to trace the pendulum
    posArr = [zeros(1, length(t)); zeros(1, length(t))];
    posArr(1,1) = x_pos;
    posArr(1,2) = y_pos;
    index = 1;    
   
 %Variables for drawing circles / masses
    th = 0:pi/50:2*pi;
    circleSize = (1/20)*m;
    stationaryCircleSize = 0.25;
    
    %{
    Single pendulum with damping:  theta'' = -b/m*theta' + -g/l*sin(theta)
    %}
    
%% Simulation
for k=2:steps
   %Runge Kutta Method  
       thetaOmega = rkCalculator(theta, omega, T, g, l, k, b, m); % Custom function
       theta(k) = thetaOmega(1);
       omega(k) = thetaOmega(2);
   
    
   % PLOTTING
   % Current Position
       x_pos = l*sin(theta(k));
       y_pos = l*cos(theta(k));
  
   % Adding to tracing array
       posArr(1, index) = xStartPoint - x_pos;
       posArr(2, index) = yStartPoint - y_pos;
   
   % Increment the array index and keep track of the last index b/c this
   % essentially works as a ring buffer
       if(index >= length(posArr))
            lastIndex = index;
            index = 1;
       else
            lastIndex = index;
            index = index + 1;
       end
   
   % Plotting
       clf;
       hold on
           % Plot the trace that follows the pendulum mass
               if(k*T <= traceLength)
                   plot(posArr(1,1:k-1), posArr(2,1:k-1))
               else
                   plot(posArr(1,k-(traceLength/T):k-1), posArr(2,k-(traceLength/T):k-1))
               end

           % Plot the pendulum rod
               plot([xStartPoint, posArr(1,lastIndex)],[yStartPoint, posArr(2,lastIndex)], 'black');

           % Draw anchored point
               circleX = stationaryCircleSize * cos(th) + xStartPoint;
               circleY = stationaryCircleSize * sin(th) + yStartPoint;
               plot(circleX, circleY, 'black');
               fill(circleX, circleY, 'black');
           % Draw the mass on the end of the pendulum
               circleX = circleSize * cos(th) + posArr(1, lastIndex);
               circleY = circleSize * sin(th) + posArr(2, lastIndex);
               plot(circleX, circleY, 'b');
               fill(circleX, circleY, 'g');
       hold off
       title(plotTitle)
       
   % Adjust camera View
       if(steadyCam  == 1)
           xlim([-1, 2*l + 1]);
           ylim([-1, 2*l + 1]);
       else
           xlim([posArr(1, lastIndex) - 2*l - 1, posArr(1, lastIndex) + 2*l + 1]);
           ylim([posArr(2, lastIndex) - 2*l - 1, posArr(2, lastIndex) + 2*l + 1]);
       end
    
   % ANIMATE
       if(mod(k,1) == 0)
            pause(0.01) 
            %disp(k)
       end
end
