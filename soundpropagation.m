% Parameters
fs = 44100; % Sampling frequency in Hz
duration = 1; % Duration of the signal in seconds
c = 343; % Speed of sound in m/s

% Microphone positions (5 reference mics and 1 for localization)
center = [0, 0, 1]; % Center position at (0,0,1)
radius = 0.4; % 40 cm from the center
angles = linspace(0, 2*pi, 6); % Equal angles for 5 mics

micPos = [radius*cos(angles(1:5))', radius*sin(angles(1:5))', ones(5,1)]; % Reference mics
micPos = [micPos; center]; % Add localization mic at the center

% Gunshot source position
sourcePos = [0.5, 0.5, 1]; % Position of the gunshot

% Animation parameters
numFrames = 200; % Number of frames in the animation
maxDistance = max(sqrt(sum((micPos - sourcePos).^2, 2))) + 0.5; % Max distance for wavefront

figure('Name', 'Sound Propagation Animation', 'NumberTitle', 'off', 'Position', [100, 100, 1600, 1200]);

% Animation loop
for tFrame = 1:numFrames
    clf; % Clear the current figure
    
    % Calculate current wavefront radius
    wavefrontRadius = (tFrame / numFrames) * maxDistance;
    
    % Plot the sound source
    scatter3(sourcePos(1), sourcePos(2), sourcePos(3), 200, 'r', 'filled', 'DisplayName', 'Gunshot Source');
    hold on;
    
    % Plot the expanding wavefront
    [x, y, z] = sphere(50);
    surf(sourcePos(1) + wavefrontRadius * x, sourcePos(2) + wavefrontRadius * y, sourcePos(3) + wavefrontRadius * z, ...
        'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    % Plot the microphones
    scatter3(micPos(:,1), micPos(:,2), micPos(:,3), 150, 'g', 'filled', 'DisplayName', 'Microphones');
    
    % Customize the plot
    axis equal;
    grid on;
    xlabel('X (m)', 'FontSize', 12);
    ylabel('Y (m)', 'FontSize', 12);
    zlabel('Z (m)', 'FontSize', 12);
    title('Sound Propagation from Gunshot Source', 'FontSize', 14);
    legend('show', 'FontSize', 10, 'Location', 'northeastoutside');
    
    % Set view angle
    view(3);
    
    % Update the plot
    drawnow;
end


