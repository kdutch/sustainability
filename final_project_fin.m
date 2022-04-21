
%kdutcher 0-27-2022 03:14

clear; clc; close all;
totalTimeLength = 10;
dt = 0.01;
dx = 0.01;
simTimeSteps = int64(totalTimeLength/dt);

%It should be noted that heats move ONLY by conduction and convection in
%this simulation. In reality, heat rising(that is, air rising as it
%decreases in density) also contributes greatly to the distribution of heat
%in the atmosphere. This air is acting like a solid, and is not moving, and
%therefore is continuously acruing heat.
 writerObj = VideoWriter('normal_8');
 writerObj.FrameRate = 30;
 writerObj.Quality = 100;
 open(writerObj);


I = imread('normal.png');
[rows, columns, numberOfColorChannels] = size(I);
redChannel = I(:, :, 1);
greenChannel = I(:, :, 2);
blueChannel = I(:, :, 3);
material = (redChannel ~= 255) & (greenChannel ~= 255) & (blueChannel ~= 255);
% https://www.researchgate.net/publication/304719050_An_Investigation_into_the_Thermal_Properties_of_Termite_Mound_Clay_Applicable_to_Grain_Silo_Construction#pf3
airSH = 718;% Specific Heat of Air - J/kg/C
metalSH = 903;% Specific Heat of aluminum) - J/kg/C 
metalTC = 251;% Thermal conductivity of aluminum - W/mK
airTC = 0.025;% Thermal conductivity of Air - W/mK
airDensity = 1.225;% Density of Air - kg/m3
metalDensity = 2713;% Density of aluminum - kg/m3
simX = columns;%x dimension of sim grid
simY = rows; %y dimension of sim grid
temp = zeros(simY, simX); %initialize temps based on image size
temp1 = zeros(simY, simX); %initialize temps
heatSource = 38; %cold Source temp
flow = zeros(simY, simX);
coldSource = 14; %heat source temp
map = zeros(simY, simX); %material map
simXhalf = int64(simX/2); %half the image
simYhalf = int64(simY/2); %half the image
map(material) = 1; %"copies" the image to be a binary yes/no
temp(:, :) = (heatSource + coldSource)/2;%average of the cold/hot source
temp(material) = 50;
combust = 16*0.1;
q_x_first = ones(simX); %this is to store where the the sun lands
%set up an initial condition of solar irradiation

%figure initializatin
f1 = figure('Position',[1000, 100, 500, 500]);
temp1 = temp;
imagesc(temp1, [14 100]);
colorbar;
drawnow;
flow = zeros(simY, simX);
%initialize variables
dTemp = 0;
q_right = 0;
q_down = 0;
q_up = zeros(simX, 1);
t_half = int64(simTimeSteps/2);
mid_source = (coldSource+heatSource)/2;
% Start Calcuations
%with initial conditions, it acts like the sun just rose and is starting to
%add heat to the system
for t = 1:simTimeSteps
    %reset values for the "first" non-air material
    q_x_first = ones(simX, 1);
    %calculate heat for first row
    for x = 1:simX
        q_up(x) = conduct(airTC, temp(1, x), temp(2, x), dt, dx);
    end
    %calculate heat for all other rows
    for y = 2:simY-1
        if (map(y, 1) == 0)
            k = airTC;
        else
            k = metalTC;
        end
        %calculate left-most column
        q_left = -conduct(k, mid_source, temp(y, 1), dt, dx);
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
               k = metalTC;
               c = metalSH;
               p = metalDensity;
               mat=1;
           end
           %calculate heat from the up/left
           q0 = q_up(x) + q_left;
           if (x >= 43 && x <= 62 && y ==96)
               if (flow(y, x) < 10)
                %a negative heat signi;fies negative heat-flow, meaning heat is flowing in
                    q0 = q0 - combust;
                    flow(y, x) = flow(y, x) + 10; %oxygen is 'consumed'
               end
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
            k = metalTC;
            p = metalDensity;
            c = metalSH;
        end
        %check if the materials are the same, and pick
        %conduction/convection
        if (map(y, simX) == map(y+1, simX))
            q_down = conduct(k, temp(y, simX), temp(y+1, simX), dt, dx);
        else
            q_down = convect(k, temp(y, simX), temp(y+1, simX), dt, dx);
        end
        %temp from the left/right (boundary sources)
        q_right = conduct(k, temp(y, simX), mid_source, dt, dx);
        q_total = q_left + q_down + q_up(simX) + q_right;
        q_up(simX) = q_down;
        m = dx^3*p;
        %set new temp
        dTemp = q2T(q_total, m, c);
        temp1(y, simX) = temp(y, simX) + dTemp;
    end  
    new_flow = fluid_flow(map, temp1, flow, simX, simY);
    temp = temp_flow(new_flow, temp1, simX, simY, mid_source);
    flow = new_flow;
    imagesc(temp, [0 100]);
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
function new_temp = temp_flow(flow, temp, simX, simY, mid_source)
    new_temp = zeros(simY, simX);
    new_temp(1, :) = mid_source;
    new_temp(:, 1) = mid_source;
    new_temp(:, simX) = mid_source;
    new_temp(simY, :) = mid_source;
    for y = 2:simY-1
        for x = 2:simX-1
            bias = 0;
            if (flow(y, x) >= 10)
                bias = 10;
            end
            switch flow(y, x)-bias
                case 0
                    new_temp(y, x) = temp(y, x);
                case 1
                    new_temp(y, x) = temp(y-1, x);
                case 2
                    new_temp(y, x) = temp(y, x-1);
                case 3
                    new_temp(y, x) = temp(y, x+1);
                case 4
                    new_temp(y, x) = temp(y+1, x);
                case 5 
                    new_temp(y, x) = temp(y, x);
            end
        end
    end
end

function new_flow = fluid_flow(mat, temp, flow, simX, simY)
   new_flow = zeros(simY, simX); %initalize
   new_flow(1, :) = 4;
   new_flow(:, 1) = 3;
   new_flow(:, simX) = 2;
   for y = 2:simY-1
       for x = 2:simX-1
           bias = 0;
           if (flow >= 10) %check if consumed
               bias = 10;
           end
           if mat(y, x) == 0 %check if air
                if (temp(y, x) > temp(y-1, x) && mat(y-1, x) == 0)
                    new_flow(y, x) = 1 + bias;
                    if (flow(y-1, x) < 10) %if not consumed
                            new_flow(y, x) = new_flow(y, x) - bias; %label not consumed
                    end
                elseif (temp(y, x) < temp(y+1, x) && mat(y+1, x) == 0) %if air and temp greater
                    new_flow(y, x) = 4 + bias; %push down
                    if (flow(y+1, x) < 10) %if air below is non consumed
                            new_flow(y, x) = new_flow(y, x) - bias; %set to non-consumed
                    end
                else
                    neighbor_total = 0; %initialize neighbor
                    if (flow(y,x-1) ~= 3) || (flow(y, x-1) ~= 13 && mat(y, x-1) == 0) %if left air isn't moving right
                        neighbor_total = neighbor_total - 1; %air wants to move left
                    elseif (flow(y, x+1) ~= 2) || (flow(y, x+1) ~= 12 && mat(y, x+1) == 0) %if air isn't moving left
                        neighbor_total = neighbor_total + 1;  %air wants to move right
                    end
                    if (neighbor_total < 0) %if nieghbor total is negative, air wants to move left
                        new_flow(y, x) = 2 + bias;
                        if (flow(y, x-1) < 10) %check consumed
                            new_flow(y, x) = new_flow(y, x) - bias;
                        end
                    elseif(neighbor_total > 0) %otherwise, air wants to move right
                        new_flow(y, x) = 3 + bias;
                        if (flow(y, x+1) < 10) %check consimed
                            new_flow(y, x) = new_flow(y, x) - bias;
                        end
                    else %if no preference (ie not left or right)
                        if abs(temp(y, x+1) - temp(y, x)) > abs(temp(y, x-1) - temp(y, x && mat(y, x-1) == 0)) %if the temperature diffence is greater
                            new_flow(y, x) = 3 + bias;
                            if (flow(y, x+1) < 10) %check consumed
                                new_flow(y, x) = new_flow(y, x) - bias;
                            end
                        elseif abs(temp(y, x-1) - temp(y, x)) > abs(temp(y, x+1) - temp(y, x)) && mat(y, x+1) == 0 %if temp dif is greater
                            new_flow(y, x) = 2 + bias;
                             if (flow(y, x-1) < 10) %check consumed
                                new_flow(y, x) = new_flow(y, x) - bias;
                            end
                        end
                        new_flow(y, x) = 0 + bias; %otherwise, air stays in place
                    end
                end
           end
       end
   end
   idx = mat==1; %get material
   new_flow(idx) = 5; %set material value ie. solid
end