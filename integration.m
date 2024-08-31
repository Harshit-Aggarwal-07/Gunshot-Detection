% Parameters
fs = 44100; % Sampling frequency in Hz
duration = 1; % Duration of the signal in seconds
c = 343; % Speed of sound in m/s

% Frequency range for the gunshot sound
f_min = 250;  % Minimum frequency in Hz
f_max = 3000; % Maximum frequency in Hz

% Generate a random frequency within the specified range for each component
numComponents = 10; % Number of frequency components
frequencies = f_min + (f_max - f_min) * rand(1, numComponents);

% Generate the gunshot waveform
t = 0:1/fs:duration-1/fs; % Time vector
waveform = zeros(size(t));
for i = 1:numComponents
    amplitude = rand; % Random amplitude for each component
    phase = 2*pi*rand; % Random phase for each component
    waveform = waveform + amplitude * sin(2*pi*frequencies(i)*t + phase);
end

% Apply a decay to simulate the fading sound of a gunshot
decayFactor = 5; % Decay factor to simulate fading
waveform = waveform .* exp(-decayFactor * t);

% Normalize the waveform to avoid clipping
waveform = waveform / max(abs(waveform));

% Microphone positions (5 reference mics and 1 for localization)
center = [0, 0, 1]; % Center position at (0,0,1)
radius = 0.4; % 40 cm from the center
angles = linspace(0, 2*pi, 6); % Equal angles for 5 mics

micPos = [radius*cos(angles(1:5))', radius*sin(angles(1:5))', ones(5,1)]; % Reference mics
micPos = [micPos; center]; % Add localization mic at the center

% Gunshot source position
sourcePos = [0.5, 0.5, 1]; % Position of the gunshot

% Calculate distances from the gunshot to each microphone
distances = sqrt(sum((micPos - sourcePos).^2, 2));
timeDelays = distances / c;

% Generate signals for each microphone with appropriate delays
micSignals = zeros(6, length(t));
for i = 1:6
    delaySamples = round(timeDelays(i) * fs);
    micSignals(i, delaySamples+1:end) = waveform(1:end-delaySamples);
end

% Simulation parameters
numFrames = 200; % Number of frames in the animation
maxDistance = max(distances) + 0.5; % Max distance for wavefront

figure('Name', 'Realistic Sound Propagation & Detection Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 1600, 1200]);

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
        'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    % Plot the microphones
    scatter3(micPos(:,1), micPos(:,2), micPos(:,3), 150, 'g', 'filled', 'DisplayName', 'Microphones');
    
    % Highlight microphones when the wavefront reaches them
    for i = 1:6
        if wavefrontRadius >= distances(i)
            scatter3(micPos(i,1), micPos(i,2), micPos(i,3), 300, 'm', 'filled', 'DisplayName', ['Mic ' num2str(i)]);
        end
    end
    
    % Customize the plot
    axis equal;
    grid on;
    xlabel('X (m)', 'FontSize', 12);
    ylabel('Y (m)', 'FontSize', 12);
    zlabel('Z (m)', 'FontSize', 12);
    title('Sound Propagation and Detection from Gunshot Source', 'FontSize', 14);
    legend('show', 'FontSize', 10, 'Location', 'northeastoutside');
    
    % Set view angle
    view(3);
    
    % Update the plot
    drawnow;
end

% Plot the propagated waveforms for each microphone
figure('Name', 'Received Signals at Microphones', 'NumberTitle', 'off', 'Position', [100, 100, 1600, 800]);
for i = 1:6
    subplot(3, 2, i);
    plot(t, micSignals(i,:), 'b');
    title(['Mic ' num2str(i) ' Signal'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Amplitude', 'FontSize', 12);
    grid on;
end
