 % Parameters
c_r = 500;      % Specific heat capacity of the radiator (J/(kg·°C))
m_r = 45;       % Mass of the radiator (kg)
h_r = 25;       % Heat transfer coefficient (W/(m^2·°C))
A_r = 1.5;        % Surface area of the radiator (m^2)
T_r0 = 14;      % Initial temperature of the radiator (°C)

T_room0 = 10;   % Initial temperature of the room (°C)
T_out = 10;     % Outside temperature (°C)
Length = 5;
Width = 5;
Height = 2.3;
V = Length * Width * Height;
c_a = 700;

m_a = V * 1.293;

n = 1.5;
target_T = 18;

ttotal= 24*3600;
dt = 10;
time = (0:10:ttotal);
Powerplot = zeros(size(time));


T_radiator = zeros(size(time));
T_radiator(1) = T_r0;

T_room = zeros(size(time));
T_room(1) = T_room0;

[heatingpower, ttarget] = findOptimumPower(target_T, T_out, T_out);
for i = 2:length(time)
    % Check if the current time is within the occupancy schedule
    if times(i) >= 9 * 3600 && times(i) <= 18 * 3600
        % unoccupied period
        power_profile = 0; % Implement your optimization algorithm
    else
        
        if T_room(i-1) < target_T
            power_profile = heatingpower;
        else
            power_profile = calculateHeatLoss(target_T);
        end
    end

    dT = (power_profile - h_r * A_r * (T_radiator(i-1) - T_room(i-1))) * dt / (m_r * c_r);
    
    HeatLoss = calculateHeatLoss(T_room(i-1)); % Implement your heat loss calculation
    
    dTroom = (h_r * A_r * (T_radiator(i-1) - T_room(i-1)) - HeatLoss) * dt / (m_a * c_a);

    T_radiator(i) = T_radiator(i-1) + dT;
    T_room(i) = T_room(i-1) + dTroom;
    Powerplot(i) = power_profile;
end



figure;
% Radiator Temperature subplot
subplot(3, 1, 1);
plot(time, T_radiator, 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Temperature of Radiator (°C)');
title('Radiator Temperature Over Time');
grid on;

% Room Temperature subplot
subplot(3, 1, 2);
plot(time, T_room, 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Room Temperature (°C)');
title('Room Temperature Over Time');
grid on;

subplot(3,1,3);
plot(time, Powerplot, 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Power (W)');
title('Power Over Time');
grid on;
ylim([0, 1600])

totalEnergyConsumed = trapz(time, Powerplot) / (3.6*10^6); % Convert from watts to kilowatt-hours
fprintf('Total energy consumed: %.2f kWh\n', totalEnergyConsumed);

function HeatLoss =  calculateHeatLoss(T_rm)
    n = 1.5;
    Length = 5;
    Width = 5;
    Height = 2.3;
    h_walls = 0.6;
    h_floor = 0.26;
    h_roof = 0.16;
    h_win = 5;
    target_T = 18;
    V = Length * Width * Height;
    T_out = 10;
    A_win = 1;
    A_floor = 25;
    A_roof = 25;
    A_w = Length * Height * 2 + Width * Height * 2 - A_win;
    HeatLoss = (h_walls * A_w + h_floor * A_floor + h_roof * A_roof + h_win * A_win + 0.33 * n * V) * (T_rm - T_out);
end 

