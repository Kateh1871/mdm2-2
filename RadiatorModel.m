% Parameters
c_r = 500;      % Specific heat capacity of the radiator (J/(kg·°C))
m_r = 45;       % Mass of the radiator (kg)
h_r = 25;       % Heat transfer coefficient (W/(m^2·°C))
A_r = 1.5;      % Surface area of the radiator (m^2)
T_r0 = 14;      % Initial temperature of the radiator (°C)

T_room0 = 14;   % Initial temperature of the room (°C)
T_out = 8;      % Outside temperature (°C)
T_target = 18;  % Target room temperature

% Room dimensions
Length = 5;    
Width = 5;
Height = 2.3;
V = Length * Width * Height;
disp(['Volume of room: ', num2str(V)])

c_a = 700;  % Specific heat capacity of air (?)
% Heat transfer coefficient of room (?)
h_walls = 0.6;
h_floor = 0.26;
h_roof = 0.16;
h_win = 5;

% I have no idea about any of this, I think it's a window
A_win = 1;
A_floor = 25;
A_roof = 25;
A_w = Length * Height * 2 + Width * Height * 2 - A_win;
m_a = V * 1.293;

n = 1.5;
% Time Settings
T_total = 7200; % Total simulation time (seconds)
dt = 10;        % Time step for simulation (seconds)
time = 0:dt:T_total;

% Input Power
P_in = 600 * ones(size(time));

% Temperature Simulation
T_radiator = zeros(size(time));
T_radiator(1) = T_r0;

T_room = zeros(size(time));
T_room(1) = T_room0;

% Manual Power and Energy Search
minPower = 0;
minTotalEnergy = inf;

% Simulate until the room temperature reaches target temperature degrees Celsius
for power = 0:10:1600  % Adjust the step size as needed
    % Simulate with current power
    for i = 2:length(time)
        dT = (power - h_r * A_r * (T_radiator(i-1) - T_room(i-1))) * dt / (m_r * c_r);

        HeatLoss = (h_walls * A_w + h_floor * A_floor + h_roof * A_roof + h_win * A_win + 0.33 * n * V) * (T_room(i-1) - T_out);
        dTroom = (h_r * A_r * (T_radiator(i-1) - T_room(i-1)) - HeatLoss) * dt / (m_a * c_a);

        T_radiator(i) = T_radiator(i-1) + dT;
        T_room(i) = T_room(i-1) + dTroom;

        % Check if the room temperature reaches the target
        if T_room(i) >= T_target
            % Calculate the total energy consumption until this point
            totalEnergy = trapz(time(1:i), P_in(1:i));

            % Update minimum power and energy if the current power yields a smaller total energy
            if totalEnergy < minTotalEnergy
                minTotalEnergy = totalEnergy;
                minPower = power;
                timetotarget = time(i);
            end

            % Break the inner loop as we reached the target temperature
            break;
        end
    end
end

% Adjust power to maintain target temperature degrees Celsius
adjustedPower = zeros(size(time));
for i = find(time > time(i-1), 1):length(time)
    % Calculate the difference between the room temperature and the target
    temperatureDifference = T_target - T_room(i-1);

    % Adjust power based on the difference
    powerAdjustment = h_r * A_r * temperatureDifference * m_a * c_a / (dt * A_w);

    % Ensure the adjusted power is non-negative
    adjustedPower(i) = max(minPower + powerAdjustment, 0);

    % Simulate with the adjusted power
    dT = (adjustedPower(i) - h_r * A_r * (T_radiator(i-1) - T_room(i-1))) * dt / (m_r * c_r);

    HeatLoss = (h_walls * A_w + h_floor * A_floor + h_roof * A_roof + h_win * A_win + 0.33 * n * V) * (T_room(i-1) - T_out);
    dTroom = (h_r * A_r * (T_radiator(i-1) - T_room(i-1)) - HeatLoss) * dt / (m_a * c_a);

    T_radiator(i) = T_radiator(i-1) + dT;
    T_room(i) = T_room(i-1) + dTroom;
end


% Plotting Results
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

% Power subplot
subplot(3, 1, 3);
plot(time, adjustedPower, 'LineWidth', 2);
hold on;
plot([0, timetotarget], [minPower, minPower], '--', 'LineWidth', 2, 'Color', 'r');
xlabel('Time (seconds)');
ylabel('Power (W)');
title('Power Over Time');
legend('Adjusted Power', 'Optimal Power');
grid on;
hold off;

% Display the optimal power and energy
disp(['Optimal Power: ' num2str(minPower)]);    % Why is this optimal power?
disp(['Min Total Energy: ' num2str(minTotalEnergy)]);
totalEnergyUsed = (trapz(time, adjustedPower) + minTotalEnergy)*2.777*10^-7;    % Where does this number come from?

disp(['Total Energy Used: ' num2str(totalEnergyUsed) ' kWh']);
