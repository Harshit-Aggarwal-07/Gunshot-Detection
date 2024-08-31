% Parameters
fs = 44100; % Sampling frequency in Hz
duration = 1; % Duration of the signals in seconds
c = 343; % Speed of sound in m/s
t = 0:1/fs:duration-1/fs; % Time vector

% Microphone positions (5 reference mics and 1 for localization)
center = [0, 0, 1]; % Center position at (0,0,1)
radius = 0.4; % 40 cm from the center
angles = linspace(0, 2*pi, 6); % Equal angles for 5 mics

micPos = [radius*cos(angles(1:5))', radius*sin(angles(1:5))', ones(5,1)]; % Reference mics
micPos = [micPos; center]; % Add localization mic at the center

% Generate synthetic gunshot signal
gunshotSignal = zeros(1, length(t));
gunshotSignal(round(0.1 * fs)) = 1; % Impulse at 0.1 seconds

% Gunshot source initial position
sourcePosStart = [1.5, 1.5, 1]; % Initial source position
sourcePosEnd = [0.5, 0.5, 1]; % Final source position

% Animation parameters
numFrames = 100; % Number of frames in the animation
movementStep = (sourcePosEnd - sourcePosStart) / numFrames; % Movement per frame

figure('Name', 'Sound Source Tracking Animation', 'NumberTitle', 'off', 'Position', [100, 100, 1600, 1200]); % Large figure size

% Animation loop
for tFrame = 1:numFrames
    clf; % Clear the current figure
    
    % Calculate current source position
    sourcePos = sourcePosStart + (tFrame - 1) * movementStep;
    
    % Calculate distances from the current source position to each microphone
    distances = sqrt(sum((micPos - sourcePos).^2, 2));
    timeDelays = distances / c;
    
    % Generate signals for each microphone with appropriate delays
    micSignals = zeros(6, length(t));
    for i = 1:6
        delaySamples = round(timeDelays(i) * fs);
        micSignals(i, delaySamples+1:end) = gunshotSignal(1:end-delaySamples);
    end
    
    % Cross-correlation and time differences (TDOA)
    timeDiffs = zeros(5, 1);
    for i = 1:5
        [corr, lags] = xcorr(micSignals(6,:), micSignals(i,:));
        [~, idx] = max(corr);
        timeDiffs(i) = lags(idx) / fs; % Time difference in seconds
    end
    
    % Estimate DoA (simplified for coplanar microphones)
    doaEstimates = asin(c * timeDiffs / radius); % Direction of arrival in radians
    
    % Estimate source position (simple average method)
    xEstimate = mean(radius * cos(doaEstimates));
    yEstimate = mean(radius * sin(doaEstimates));
    zEstimate = 1; % Assuming the source is at the same height as the microphones
    
    estimatedPos = [xEstimate, yEstimate, zEstimate];
    
    % Plot Environment and Estimated Position
    hold on;
    
    % Plot terrain (simple example: flat ground with elevation)
    [x, y] = meshgrid(-2:0.1:2, -2:0.1:2);
    z = zeros(size(x)); % Flat ground
    surf(x, y, z, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Plot microphones
    scatter3(micPos(:,1), micPos(:,2), micPos(:,3), 150, 'g', 'filled', 'DisplayName', 'Microphones');
    
    % Plot actual and estimated source positions
    scatter3(sourcePos(1), sourcePos(2), sourcePos(3), 200, 'r', 'filled', 'DisplayName', 'Actual Gunshot');
    scatter3(estimatedPos(1), estimatedPos(2), estimatedPos(3), 200, 'b', 'filled', 'DisplayName', 'Estimated Source');
    
    % Plot connections from microphones to estimated source
    for i = 1:5
        plot3([micPos(i,1), estimatedPos(1)], [micPos(i,2), estimatedPos(2)], [micPos(i,3), estimatedPos(3)], ...
            'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', ['Mic' num2str(i)]);
    end
    
    % Annotations
    text(sourcePos(1), sourcePos(2), sourcePos(3) + 0.1, 'Actual Gunshot', 'FontSize', 10, 'Color', 'r');
    text(estimatedPos(1), estimatedPos(2), estimatedPos(3) + 0.1, 'Estimated Source', 'FontSize', 10, 'Color', 'b');
    
    % Customize Appearance
    title('Sound Source Tracking with Sensor Network - Animation', 'FontSize', 14);
    xlabel('X (m)', 'FontSize', 12);
    ylabel('Y (m)', 'FontSize', 12);
    zlabel('Z (m)', 'FontSize', 12);
    legend('show', 'FontSize', 10, 'Location', 'northeastoutside');
    axis equal;
    grid on;
    view(3); % Set view angle
    
    % Update the plot
    drawnow;
end

