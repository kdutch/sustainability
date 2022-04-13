
clear
close all
%space variables (ENTER)
Nx=100;           %number of columns
Ny=50;           %number of rows 
%discretization variables (ENTER)
dx = 1;          % x-grid size
dy = 1;          % y-grid size
%fluid density & viscosity (ENTER)
density = 1000;
viscosity = 1;
ratio = density/viscosity;       %ratio of density and dynamic viscosity 
%size and dimension of pressure and velocity (AUTO)
p = zeros(Ny,Nx);      %pressure
u = zeros(Ny+1,Nx+1);  %x-velocity  
v = zeros(Ny+1,Nx+1);  %y-velocity
residual = zeros(Ny*Nx,1); %residuals from continuity
dp = zeros(Ny*Nx,1);  %changes in pressures
%initial conditions (ENTER)
u = zeros(Ny+1,Nx+1)+1;
%kdutcher 0-27-2022 03:14

clear; clc; close all;
totalTimeLength = 1000;
dt = 0.75;
dx = 0.01;
simTimeSteps = int64(totalTimeLength/dt);

%It should be noted that heats move ONLY by conduction and convection in
%this simulation. In reality, heat rising(that is, air rising as it
%decreases in density) also contributes greatly to the distribution of heat
%in the atmosphere. This air is acting like a solid, and is not moving, and
%therefore is continuously acruing heat.
 writerObj = VideoWriter('normal_stove');
 writerObj.FrameRate = 30;
 writerObj.Quality = 100;
 open(writerObj);


I = imread('normal.png');
[rows, columns, numberOfColorChannels] = size(I);
redChannel = I(:, :, 1);
greenChannel = I(:, :, 2);
blueChannel = I(:, :, 3);
material = (redChannel == 255) & (greenChannel == 255) & (blueChannel == 255);
% https://www.researchgate.net/publication/304719050_An_Investigation_into_the_Thermal_Properties_of_Termite_Mound_Clay_Applicable_to_Grain_Silo_Construction#pf3
airSH = 1003.5;% Specific Heat of Air - J/kg/C
tmcSH = 2576.9;% Specific Heat of TMC (Termite Mound Clay) - J/kg/C 
tmcTC = 0.21;% Thermal conductivity of TMC - W/mK
airTC = 0.025;% Thermal conductivity of Air - W/mK
airDensity = 1.225;% Density of Air - kg/m3
tmcDensity = 1833;% Density of TMC - kg/m3
airResist = 1/(airSH * airDensity);% Air Material Resistance
tmcResist = 1/(tmcSH * tmcDensity);% TMC Material Resistance
airCon = airTC * (totalTimeLength/simTimeSteps); % Air Adjusted Conductivity
tmcCon = tmcTC * (totalTimeLength/simTimeSteps); % TMC Adjusted Conductivity
simX = columns;%x dimension of sim grid (adjusted to termite mound image)
simY = rows; %y dimension of sim grid (adjusted to termite mound image)
temp = zeros(simY, simX); %initialize temps based on image size
temp1 = zeros(simY, simX); %initialize temps
heatSource = 38; %cold Source temp
flow = zeros(simY, simX)
coldSource = 14; %heat source temp
map = zeros(simY, simX); %material map
simXhalf = int64(simX/2); %half the image
simYhalf = int64(simY/2); %half the image
map(~material) = 1; %"copies" the image to be a binary yes/no where the termite mound is
temp(:, :) = (heatSource + coldSource)/2;%average of the cold/hot source
temp(material) = 14;
temp(simY, :) = coldSource;
solar = 425*dt*dx^2*0.3; %energy re-radiated by the surface and absorbed by the atmosphere 
%solar = 5*dt*dx^2;
q_x_first = ones(simX); %this is to store where the the sun lands
%set up an initial condition of solar irradiation

%figure initializatin
f1 = figure('Position',[1000, 100, 500, 500]);
temp1 = temp;
imagesc(temp1, [14 38]);
colorbar;
drawnow;

%initialize variables
dTemp = 0;
q_right = 0;
q_down = 0;
q_up = zeros(simX, 1);
t_half = int64(simTimeSteps/2);
midSource = (coldSource+heatSource)/2;
% Start Calcuations
%with initial conditions, it acts like the sun just rose and is starting to
%add heat to the system
for t = 1:simTimeSteps
    %reset values for the "first" non-air material
    q_x_first = ones(simX, 1);
    temp(340:350, 150:160) = 35;
    temp(430:440, 125:135) = 35;
    temp(430:440, 175:185) = 35;
    %calculate heat for first row
    for x = 1:simX
        q_up(x) = conduct(airTC, temp(1, x), temp(2, x), dt, dx);
    end
    %calculate heat for all other rows
    for y = 2:simY-1
        if (map(y, 1) == 0)
            k = airTC;
        else
            k = tmcTC;
        end
        
        %calculate left-most column
        q_left = -conduct(k, midSource, temp(y, 1), dt, dx);
        %reset up/down (heat going up is opposite heat going down)
        q_up = -1*q_up;
        q_x_first(simX) = 0;
        %go across all columns
        for x = 1:simX-1
            %check materiaals
           if map(y, x) == 0
               k = airTC;
               c = airSH;
               p = airDensity;
               mat=0;
           else
               k = tmcTC;
               c = tmcSH;
               p = tmcDensity;
               mat=1;
           end
           %calculate heat from the up/left
           q0 = q_up(x) + q_left;
           if (x >= 43 && x <=65 && y ==96)
                %a negative heat signifies negative heat-flow, meaning heat is flowing in
                q0 = q0 - solar;  
           end
           %calculate the difference in temperature, and heat from down and
           %the right
           [dTemp, q_right, q_down] = heat(q0, dt, dx, temp(y, x), temp(y, x+1), temp(y+1, x), k, p, c, map(y, x), map(y, x+1), map(y+1, x));
           q_left = -q_right;
           q_up(x) = q_down;
           temp1(y, x) = temp(y, x) + dTemp;
        end
        %deal with the last row of x (needs to be done seperately)
        if (map(y, simX) == 0)
            k = airTC;
            p = airDensity;
            c = airSH;
        else
            k = tmcTC;
            p = tmcDensity;
            c = tmcSH;
        end
        %check if the materials are the same, and pick
        %conduction/convection
        if (map(y, simX) == map(y+1, simX))
            q_down = conduct(k, temp(y, simX), temp(y+1, simX), dt, dx);
        else
            q_down = convect(k, temp(y, simX), temp(y+1, simX), dt, dx);
        end
        %temp from the left/right (boundary sources)
        q_right = conduct(k, temp(y, simX), midSource, dt, dx);
        q_total = q_left + q_down + q_up(simX) + q_right;
        q_up(simX) = q_down;
        m = dx^3*p;
        %set new temp
        dTemp = q2T(q_total, m, c);
        temp1(y, simX) = temp(y, simX) + dTemp;
    end  
    %set new temps
    temp = fluid(temp1, map, t, simY, simX, midSource, flow); %fluid dynamics
    %temp = fluid(temp, map, t, simY, simX, midSource); %fluid dynamics
    imagesc(temp, [14 38]);
    colorbar
    frame = getframe;
    writeVideo(writerObj,frame);
    drawnow;

end


function [dT, q_right, q_down] = heat(q0, dt, dx, T1, T2, T3, k, p, c, m_self, m_right, m_down)
    %calculate right cell
    if (m_self == m_right)
        q_right = conduct(k, T1, T2, dt, dx);
    else
        q_right = convect(k, T1, T2, dt, dx);
    end
    %calculate bottom cell
    if (m_self == m_down)
        q_down = conduct(k, T1, T3, dt, dx);
    else
        q_down = convect(k, T1, T3, dt, dx);
    end
    q_total = q_down + q_right + q0;
    m = dx^3*p;
    dT = q2T(q_total, m, c);
end

%simple heat and temperatue equations
function dTemp = q2T(q, m, c)
    dTemp = -q/(m*c);
end

function q = conduct(k, T1, T2, dt, dx)
    q = k*dx*(T1-T2)*dt;
end
    
function q = convect(k, T1, T2, dt, dx)
    q = k*dx^2*(T1-T2)*dt;
end



%just a lot of if/else statements that swap temps to mimicing moving air
function new_temp = fluid(temp, map, t, simY, simX, midSource, flow)
    temp(1, :) = midSource;
    %temp(:, 1) = midSource;
    %temp(:, simX) = midSource;
    new_temp = temp;
    bias = (-1)^(t);
    for y = 2:simY-1
        
        if (map(y, 1) == 0) %check if air
            if (map(y-1, 1) == 0) %check if lower is also air
            	if (temp(y, 1) < temp(y-1, 1))
                	new_temp(y, 1) = temp(y-1, 1);
                	new_temp(y-1, 1) = temp(y, 1);
                	temp(y, 1) = new_temp(y, 1);
                	temp(y-1, 1) = new_temp(y-1, 1);
                end
            end
            if (map(y+1, 1) == 0) %check if lower is also air
            	if (temp(y, 1) < temp(y+1, 1))
                	new_temp(y, 1) = temp(y-1, 1);
                	new_temp(y+1, 1) = temp(y, 1);
                	temp(y, 1) = new_temp(y, 1);
                	temp(y+1, 1) = new_temp(y+1, 1);
                end
            end
        end
        for x = 3:simX-1
            
           if (map(y, x) == 0)
                
                if (map(y-1, x) == 0) %check if lower is also air
                    if (temp(y, x) > temp(y-1, x))
                        new_temp(y, x) = temp(y-1, x);
                        new_temp(y-1, x) = temp(y, x);
                        temp(y, x) = new_temp(y, x);
                        temp(y-1, x) = new_temp(y-1, x);
                    end
                    if (map(y, x+bias) == 0) %if it doesnt
                        if (temp(y, x) > temp(y, x+bias))
                            new_temp(y, x) = temp(y, x+bias);
                            new_temp(y, x+bias) = temp(y, x);
                            temp(y, x) = new_temp(y, x);
                            temp(y, x+bias) = new_temp(y, x+bias);
                        end
                    end
                    if (map(y, x+bias) == 0) %if it doesnt
                        if (temp(y, x) > temp(y, x+bias))
                            new_temp(y, x) = temp(y, x+bias);
                            new_temp(y, x+bias) = temp(y, x);
                            temp(y, x) = new_temp(y, x);
                            temp(y, x+bias) = new_temp(y, x+bias);

                        end
                    end
                end
                 %check if air
                if (map(y, x+bias) == 0) %if it doesnt
                    if (temp(y, x) > temp(y, x+bias))
                        new_temp(y, x) = temp(y, x+bias);
                        new_temp(y, x+bias) = temp(y, x);
                        temp(y, x) = new_temp(y, x);
                        temp(y, x+bias) = new_temp(y, x+bias);

                    end
                end
                if (map(y, x-bias) == 0 )
                    if (temp(y, x) > temp(y, x-bias))
                        new_temp(y, x) = temp(y, x-bias);
                        new_temp(y, x-bias) = temp(y, x);
                        temp(y, x) = new_temp(y, x);
                        temp(y, x-bias) = new_temp(y, x-bias);
                    end
               end 
                if (map(y+1, x) == 0)
                    if (temp(y+1, x) < temp(y, x))
                        midTemp = (temp(y+1, x) + temp(y, x))/2;
                        temp(y, x) = midTemp;
                        temp(y+1, x) = midTemp;
                        new_temp(y, x) = midTemp;
                        new_temp(y+1, x) = midTemp;
                    end
                end
            end
        %x = x+3;
        end
        if (map(y, simX) == 0) %check if air
            if (map(y+1, simX) == 0) %check if lower is also air
            	if (temp(y, simX) < temp(y+1, simX))
                	new_temp(y, simX) = temp(y-1, simX);
                	new_temp(y+1, simX) = temp(y, simX);
                	temp(y, simX) = new_temp(y, simX);
                	temp(y+1, simX) = new_temp(y+1, simX);
                end
            end
            if (map(y-1, simX) == 0) %check if lower is also air
            	if (temp(y, simX) > temp(y-1, simX))
                	new_temp(y, simX) = temp(y-1, simX);
                	new_temp(y-1, simX) = temp(y, simX);
                	temp(y, simX) = new_temp(y, simX);
                	temp(y-1, simX) = new_temp(y-1, simX);
                end
            end
        end
    %y = y+2;
    end
end