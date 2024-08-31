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

% Calculate distances from the gunshot to each microphone
distances = sqrt(sum((micPos - center).^2, 2));
timeDelays = distances / c;

% Generate signals for each microphone with appropriate delays
micSignals = zeros(6, length(t));
for i = 1:6
    delaySamples = round(timeDelays(i) * fs);
    micSignals(i, delaySamples+1:end) = gunshotSignal(1:end-delaySamples);
end

% Cross-correlation between the localization mic and each reference mic
corrResults = zeros(5, length(micSignals(1,:)) * 2 - 1);
lags = zeros(5, length(micSignals(1,:)) * 2 - 1);
timeDiffs = zeros(5, 1);

for i = 1:5
    [corrResults(i,:), lags(i,:)] = xcorr(micSignals(6,:), micSignals(i,:));
    [~, idx] = max(corrResults(i,:));
    timeDiffs(i) = lags(i,idx) / fs; % Time difference in seconds
end

% Display time differences
disp('Time differences between the localization mic and reference mics:');
disp(timeDiffs);

% 3D Plot of the Environment with Sensor Network and Environmental Features
figure('Name', '3D Environment with Sensor Network and Terrain', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
hold on;

% Plot terrain (simple example: flat ground with elevation)
[x, y] = meshgrid(-2:0.1:2, -2:0.1:2);
z = zeros(size(x)); % Flat ground
z(x.^2 + y.^2 > 1.5^2) = NaN; % Circular terrain feature (e.g., a pond)

% Plot terrain
surf(x, y, z, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Plot microphone coverage areas (as semi-transparent cones)
coneHeight = 0.2; % Height of the cone
coneRadius = 0.4; % Radius of the cone
theta = linspace(0, 2*pi, 50);
x_cone = coneRadius * cos(theta);
y_cone = coneRadius * sin(theta);

for i = 1:5
    % Plot coverage area for each reference mic
    fill3(x_cone + micPos(i,1), y_cone + micPos(i,2), coneHeight * ones(size(x_cone)) + micPos(i,3), 'g', 'FaceAlpha', 0.2);
end

% Plot connections and positions
plot3([center(1), micPos(:,1)'], [center(2), micPos(:,2)'], [center(3), micPos(:,3)'], '-o', ...
    'LineWidth', 2, 'MarkerSize', 10, 'Color', [0 0.4470 0.7410]); % Connections
scatter3(center(1), center(2), center(3), 200, 'r', 'filled', 'DisplayName', 'Gunshot'); % Gunshot
scatter3(micPos(:,1), micPos(:,2), micPos(:,3), 150, 'g', 'filled', 'DisplayName', 'Microphones'); % Microphones

% Annotations
text(center(1), center(2), center(3) + 0.1, 'Gunshot', 'FontSize', 10, 'Color', 'r');

% Customize Appearance
title('3D Environment: Sensor Network and Terrain Features', 'FontSize', 14);
xlabel('X (m)', 'FontSize', 12);
ylabel('Y (m)', 'FontSize', 12);
zlabel('Z (m)', 'FontSize', 12);
legend('show', 'FontSize', 10, 'Location', 'northeastoutside');
axis equal;
grid on;
view(3);

% Add lighting and shading for a more realistic effect
camlight('headlight');
lighting gouraud;

% Animate Gunshot Movement
figure('Name', 'Gunshot Animation', 'NumberTitle', 'off', 'Position', [100, 100, 1600, 1200]); % Increased figure size
numFrames = 50; % Number of frames in the animation (increased speed)
startDistance = 2.0; % Start distance of the gunshot (outside coverage area)
endDistance = 0.4; % End distance of the gunshot (inside coverage area)
movementStep = (startDistance - endDistance) / numFrames; % Movement per frame

% Initial gunshot position (start outside the coverage area)
gunshotPos = [0, 0, 1 + startDistance]; 

% Animation loop
for t = 0.5:numFrames
    clf; % Clear the current figure
    
    hold on;
    % Plot terrain
    surf(x, y, z, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Plot coverage areas
    for i = 1:5
        fill3(x_cone + micPos(i,1), y_cone + micPos(i,2), coneHeight * ones(size(x_cone)) + micPos(i,3), 'g', 'FaceAlpha', 0.2);
    end

    % Update the gunshot position
    gunshotPos(3) = gunshotPos(3) - movementStep; % Move gunshot towards coverage area
    
    % Plot the gunshot wavefront
    waveRadius = (numFrames - t) * 0.02; % Expand the wavefront radius (adjusted for speed)
    [x_sphere, y_sphere, z_sphere] = sphere; % Create sphere data
    surf(waveRadius*x_sphere + gunshotPos(1), waveRadius*y_sphere + gunshotPos(2), waveRadius*z_sphere + gunshotPos(3), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.1 0.1 0.9]); % Plot wavefront
    
    % Plot microphones and connections
    plot3([gunshotPos(1), micPos(:,1)'], [gunshotPos(2), micPos(:,2)'], [gunshotPos(3), micPos(:,3)'], '-o', 'LineWidth', 2, 'MarkerSize', 10, 'Color', [0 0.4470 0.7410]);
    scatter3(gunshotPos(1), gunshotPos(2), gunshotPos(3), 200, 'r', 'filled', 'DisplayName', 'Gunshot'); % Gunshot
    scatter3(micPos(:,1), micPos(:,2), micPos(:,3), 150, 'g', 'filled', 'DisplayName', 'Microphones'); % Microphones

    % Add annotations and labels
    text(gunshotPos(1), gunshotPos(2), gunshotPos(3) + 0.1, 'Gunshot', 'FontSize', 10, 'Color', 'r');
    
    % Customize appearance
    title('3D Environment: Gunshot Propagation Animation', 'FontSize', 14);
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

% Plot distance insights using color-coded spheres
figure('Name', 'Distance Insights', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]); % Adjusted figure size
hold on;
% Create a grid for distance visualization
[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);
Z = zeros(size(X));
distances = sqrt(X.^2 + Y.^2);
scatter3(X(:), Y(:), Z(:), 50, distances(:), 'filled', 'MarkerEdgeColor', 'k');
colorbar;
caxis([0 max(distances(:))]); % Set color axis limits

% Plot microphones
scatter3(micPos(:,1), micPos(:,2), micPos(:,3), 200, 'g', 'filled', 'DisplayName', 'Microphones'); % Microphones
scatter3(center(1), center(2), center(3), 300, 'r', 'filled', 'DisplayName', 'Gunshot'); % Gunshot

% Annotations
text(center(1), center(2), center(3) + 0.1, 'Gunshot', 'FontSize', 12, 'Color', 'r');

% Customize Appearance
title('Distance Insights from Gunshot to Sensors', 'FontSize', 16);
xlabel('X (m)', 'FontSize', 14);
ylabel('Y (m)', 'FontSize', 14);
zlabel('Z (m)', 'FontSize', 14);
legend('show', 'FontSize', 12, 'Location', 'northeastoutside');
axis equal;
grid on;
view(3); % Set view angle




