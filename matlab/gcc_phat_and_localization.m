% Localization process with GCC-PHAT

if gunshot_detected
    % Example usage for loading and processing microphone signals
    fs = 44100;  % Sampling frequency
    num_mics = 5;  % Number of microphones
    mic_signals = cell(1, num_mics);
    
    % Load the microphone signals
    for i = 1:num_mics
        [mic_signals{i}, fs] = audioread(sprintf('digital_mic%d.wav', i));
    end
    
    % Calculate TDOAs using GCC-PHAT
    tdoas = zeros(1, num_mics - 1);
    for i = 2:num_mics
        tdoas(i-1) = gcc_phat(mic_signals{i}, mic_signals{1}, fs);
    end
    
    % Microphone positions in 3D space
    mic_positions = [0, 10.2, 0;        % Microphone 1
                     0, -4.8, -8.66;    % Microphone 2
                     10, -4.8, 8.66;    % Microphone 3
                     -10, -4.8, 8.66;   % Microphone 4
                     -40, 0, 0];        % Microphone 5
    
    speed_of_sound = 343;  % Speed of sound in m/s
    estimated_position = tdoa_localization(mic_positions, tdoas, speed_of_sound);
    
    % Example usage of cuckoo search for optimization
    objective_function = @(pos) norm(tdoa_localization(mic_positions, tdoas, speed_of_sound) - pos);
    bounds = [[estimated_position(1)-5,estimated_position(2)-5,estimated_position(3)-5]; [estimated_position(1)+5,estimated_position(2)+5,estimated_position(3)+5]];  % Example bounds for x, y, z
    initial_guess = estimated_position;
    num_nests = 50;
    num_iterations = 100;
    best_position = cuckoo_search(objective_function, initial_guess, bounds, num_nests, num_iterations);
    
    %% Angle Calculation and Display
    [azimuth, elevation] = calculate_angles(best_position);
    disp(['Azimuth: ', num2str(azimuth), ' degrees']);
    disp(['Elevation: ', num2str(elevation), ' degrees']);
    
    %% Distance Calculation
    droneX = 0; droneY = 0; droneZ = 0;  % Assuming drone at origin
    distance = sqrt((best_position(1) - droneX)^2 + ...
                    (best_position(2) - droneY)^2 + ...
                    (best_position(3) - droneZ)^2);
    
    fprintf('Distance from Drone: %.2f meters\n', distance);
    
    %% Plot Azimuth and Elevation
    figure;
    polarplot([0 deg2rad(azimuth)], [0 1], 'r', 'LineWidth', 2);  % Azimuth plot
    hold on;
    polarplot([0 deg2rad(elevation)], [0 1], 'b', 'LineWidth', 2);  % Elevation plot (polar is approximation)
    title('Azimuth (red) and Elevation (blue) Angles');
    legend('Azimuth', 'Elevation');
    hold off;
    
    %% Display Estimated Position
    disp('Initially Estimated Position (x, y, z):');
    disp(estimated_position);
    disp('Optimized Estimated Position (x, y, z):');
    disp(best_position);
end

%% Function Definitions - Move these OUTSIDE the if-block

% GCC-PHAT function
function tdoa = gcc_phat(x1, x2, fs)
    N = length(x1) + length(x2) - 1;
    
    % Cross power spectrum
    X1 = fft(x1, N);
    X2 = fft(x2, N);
    R = X1 .* conj(X2);
    
    % GCC-PHAT
    cc = ifft(R ./ abs(R));
    
    % Finding the peak
    [~, maxIndex] = max(abs(cc));
    
    % Adjust for MATLAB 1-based indexing
    if maxIndex > N/2
        maxIndex = maxIndex - N;
    end
    
    % Calculate TDOA
    tdoa = maxIndex / fs;
end

% TDOA Localization Function
function estimated_position = tdoa_localization(mic_positions, tdoas, speed_of_sound)
    % Cost function for least squares optimization
    function residuals = cost_function(position)
        distances = sqrt(sum((mic_positions - position).^2, 2));
        predicted_tdoas = (distances - distances(1)) / speed_of_sound;
        residuals = tdoas - predicted_tdoas;
    end

    % Initial guess based on average mic position
    initial_guess = mean(mic_positions);
    
    % Options for optimization
    options = optimoptions('lsqnonlin', 'Display', 'off');
    
    % Run optimization to find estimated position
    estimated_position = lsqnonlin(@cost_function, initial_guess, [], [], options);
end

% Cuckoo Search Algorithm
function best_position = cuckoo_search(objective_function, initial_guess, bounds, num_nests, num_iterations)
    % Parameters
    pa = 0.25;  % Discovery rate of alien eggs
    alpha = 0.01;  % Step size
    beta = 1.5;  % Lévy flight exponent
    
    % Initialize nests
    nests = initial_guess + rand(num_nests, length(initial_guess)) .* (bounds(2,:) - bounds(1,:));
    fitness = arrayfun(@(i) objective_function(nests(i,:)), 1:num_nests);
    
    for iter = 1:num_iterations
        for i = 1:num_nests
            % Lévy flight for generating new solutions
            step = alpha * (nests(i,:) - initial_guess) .* randn(1, length(initial_guess)) .* abs(randn(1).^(-1 / beta));
            new_nest = nests(i,:) + step;
            new_nest = min(max(new_nest, bounds(1,:)), bounds(2,:));  % Ensure bounds
            
            % Evaluate fitness of the new solution
            new_fitness = objective_function(new_nest);
            if new_fitness < fitness(i)
                nests(i,:) = new_nest;
                fitness(i) = new_fitness;
            end
        end
        
        % Discovery of alien eggs
        for i = 1:num_nests
            if rand() < pa
                nests(i,:) = bounds(1,:) + rand(1, length(initial_guess)) .* (bounds(2,:) - bounds(1,:));
            end
        end
        
        fitness = arrayfun(@(i) objective_function(nests(i,:)), 1:num_nests);  % Update fitness
    end
    
    % Return the best solution
    [~, idx] = min(fitness);
    best_position = nests(idx, :);
end

% Angle Calculation
function [azimuth, elevation] = calculate_angles(position)
    azimuth = atan2d(position(2), position(1));
    elevation = atan2d(position(3), sqrt(position(1)^2 + position(2)^2));
end
