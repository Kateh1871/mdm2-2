function [optimumPower, timetotarget] = findOptimumPower(target_T, T_rad, T_rm)
    c_r = 500;      % Specific heat capacity of the radiator (J/(kg·°C))
    m_r = 45;       % Mass of the radiator (kg)
    h_r = 25;       % Heat transfer coefficient (W/(m^2·°C))
    A_r = 1.5;        % Surface area of the radiator (m^2)
    T_r0 = 14;      % Initial temperature of the radiator (°C)
    
    T_room0 = 14;   % Initial temperature of the room (°C)
    T_out = 8;     % Outside temperature (°C)
    Length = 5;
    Width = 5;
    Height = 2.3;
    V = Length * Width * Height;
    c_a = 700;
    h_walls = 0.6;
    h_floor = 0.26;
    h_roof = 0.16;
    h_win = 5;
    
    A_win = 1;
    A_floor = 25;
    A_roof = 25;
    A_w = Length * Height * 2 + Width * Height * 2 - A_win;
    m_a = V * 1.293;
    
    n = 1.5;
    minTotalEnergy = Inf;
    timetotarget = 0;
    
    
    time = 3600; % Total simulation time (seconds)
    dt = 10;        % Time step for simulation (seconds)
    time = 0:dt:time;
    adjustedPower = zeros(size(time));
    % Temperature Simulation
    T_radiator = zeros(size(time));
    T_radiator(1) = T_rad;
    
    T_room = zeros(size(time));
    T_room(1) = T_rm;
    
    % Manual Power and Energy Search
    minPower = 0;

    % Outer Loop
    for power = 0:10:1600
        P_in = power * ones(size(time));
        for i = 2:length(time)
            dT = (power - h_r * A_r * (T_radiator(i-1) - T_room(i-1))) * dt / (m_r * c_r);
    
        HeatLoss = (h_walls * A_w + h_floor * A_floor + h_roof * A_roof + h_win * A_win + 0.33 * n * V) * (T_room(i-1) - T_out);
        dTroom = (h_r * A_r * (T_radiator(i-1) - T_room(i-1)) - HeatLoss) * dt / (m_a * c_a);
    
        T_radiator(i) = T_radiator(i-1) + dT;
        T_room(i) = T_room(i-1) + dTroom;
    
        % Check if the room temperature reaches the target
        if T_room(i) >= target_T
            % Calculate the total energy consumption until this point            
            totalEnergy = trapz(time(1:i), P_in(1:i));
            
            % Update minimum power and energy if the current power yields a smaller total energy
            if totalEnergy < minTotalEnergy
                minTotalEnergy = totalEnergy;
                
                timetotarget = time(i);
                adjustedPower(1:i) = power;
            end
            
            break;        
        end
    end
    
    % Return the optimum power
    optimumPower = adjustedPower(1);
    end
end
