% BEAMFORMING AND MUSIC ALGORITHM FOR GUNSHOT DETECTION

% Load signals
num_mics = 5;
mic_signals = cell(1, num_mics);
for i = 1:num_mics
    [mic_signals{i}, fs] = audioread(sprintf('Anc_output_mic%d.wav', i));
    mic_signals{i} = mic_signals{i}(:); % Ensure column vector
end

% Ensure all mic signals are the same length
signal_length = min(cellfun(@length, mic_signals)); % Use the shortest length
for i = 1:num_mics
    mic_signals{i} = mic_signals{i}(1:signal_length);
end

% Beamforming (optional, for enhancing signals)
mic_positions = [0, 10.2, 0; 0, -4.8, -8.66; 10, -4.8, 8.66; -10, -4.8, 8.66; -40, 0, 0];
theta = -30:1:30; % Reduced azimuth angles
phi = -30:1:30; % Reduced elevation angles
N_theta = length(theta);
N_phi = length(phi);

% Initialize the beamformed signal
beamformed_signal = zeros(signal_length, N_theta, N_phi);

% Time vector
t = (0:signal_length-1) / fs;

% Beamforming computation
for theta_idx = 1:N_theta
    for phi_idx = 1:N_phi
        % Convert angle indices to actual angles
        theta_rad = deg2rad(theta(theta_idx));
        phi_rad = deg2rad(phi(phi_idx));
        
        % Steering vector calculation
        direction = [cos(phi_rad)*cos(theta_rad); cos(phi_rad)*sin(theta_rad); sin(phi_rad)];
        steering_vector = exp(-1j * 2 * pi * (mic_positions * direction));
        
        % Accumulate the beamformed signal
        for i = 1:num_mics
            x = mic_signals{i}; % x is already a column vector
            beamformed_signal(:, theta_idx, phi_idx) = beamformed_signal(:, theta_idx, phi_idx) + real(x .* conj(steering_vector(i)));
        end
    end
end

% After the loop, sum across all angles to get the final beamformed signal
beamformed_signal = sum(sum(beamformed_signal, 2), 3);

% Normalize beamformed signal
beamformed_signal = beamformed_signal / max(abs(beamformed_signal)); % Normalize

% Save beamformed signal (optional, for debugging)
audiowrite('beamformed_signal.wav', beamformed_signal, fs);

% Gunshot Detection: Simple peak detection
[~, gunshot_indices] = findpeaks(beamformed_signal, 'MinPeakHeight', 0.5); % Adjust MinPeakHeight as needed

% MUSIC Algorithm for DOA estimation
R = zeros(num_mics);
for i = 1:num_mics
    for j = i:num_mics
        % Compute correlation for each pair of microphones
        R(i, j) = sum(mic_signals{i} .* conj(mic_signals{j}));
        R(j, i) = conj(R(i, j)); % Ensure Hermitian symmetry
    end
end
R = R / signal_length;

% Eigenvalue decomposition
[eig_vectors, eig_values] = eig(R);
[~, sorted_indices] = sort(diag(eig_values), 'descend');
eig_vectors = eig_vectors(:, sorted_indices);

% Assume number of sources
num_sources = 2; % Adjust based on your scenario

% MUSIC Spectrum calculation
music_spectrum = zeros(N_theta, N_phi);

for i = 1:N_theta
    for j = 1:N_phi
        % Convert angles to radians
        theta_rad = deg2rad(theta(i));
        phi_rad = deg2rad(phi(j));
        
        % Steering vector for MUSIC
        direction = [cos(phi_rad)*cos(theta_rad); cos(phi_rad)*sin(theta_rad); sin(phi_rad)];
        steering_vector = exp(-1j * 2 * pi * (mic_positions * direction));
        
        % MUSIC Spectrum calculation
        noise_subspace = eig_vectors(:, num_sources+1:end);
        music_spectrum(i, j) = 1 / (steering_vector' * (noise_subspace * noise_subspace') * steering_vector);
    end
end

% Convert MUSIC Spectrum to magnitude for peak detection
music_spectrum_mag = abs(music_spectrum);

% Find Peaks in MUSIC Spectrum Magnitude
[music_peaks, peak_indices] = findpeaks(music_spectrum_mag(:));
[peak_theta_idx, peak_phi_idx] = ind2sub(size(music_spectrum_mag), peak_indices);

% Convert indices to angles
detected_angles = [theta(peak_theta_idx), phi(peak_phi_idx)];

% Display detected angles
disp('Detected angles (Azimuth, Elevation):');
disp(detected_angles);

% Plot MUSIC Spectrum
figure;
imagesc(theta, phi, 10*log10(music_spectrum_mag));
title('MUSIC Spectrum');
xlabel('Azimuth Angle (degrees)');
ylabel('Elevation Angle (degrees)');
colorbar;
axis xy; % To flip the y-axis and have low values at the bottom
