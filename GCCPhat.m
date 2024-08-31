% GCC-PHAT Localization Script with Specified Microphone Coordinates
% Description: This script implements GCC-PHAT for localizing gunshots using
% a microphone array on a drone with 5 specified microphone coordinates and predefined audio files.

clc;
clear;
close all;

% Parameters
fs = 44100; % Sampling rate in Hz
window_size_ms = 500; % Window size in milliseconds
overlap = 0.5; % 50% overlap
speed_of_sound = 343; % Speed of sound in m/s

% Microphone coordinates for 5 microphones
mic_coords = [0, 10.2, 0;        % Microphone 1
              0, -4.8, -8.66;    % Microphone 2
              10, -4.8, 8.66;    % Microphone 3
              -10, -4.8, 8.66;   % Microphone 4
              -40, 0, 0];        % Microphone 5

% Generate microphone pairs (all combinations)
mic_pairs = nchoosek(1:5, 2);

% Load audio data for each microphone from predefined audio files
num_mics = size(mic_coords, 1);
mic_data = cell(num_mics, 1);

% Load audio data from files named microphone1.wav to microphone5.wav
disp('Loading audio data for each microphone...');
for i = 1:num_mics
    % Construct the filename for each microphone
    filename = sprintf('microphone%d.wav', i);
    
    % Load the audio file
    if isfile(filename)
        [audio, fs_audio] = audioread(filename);
        
        % Check if the sample rate matches the desired fs
        if fs_audio ~= fs
            error('Sampling rate of the audio file %s does not match the specified fs of %d Hz.', filename, fs);
        end
        
        mic_data{i} = audio(:, 1); % Assuming mono audio or using the first channel
    else
        error('Audio file %s not found.', filename);
    end
end

% Ensure all audio signals are of the same length
min_length = min(cellfun(@length, mic_data));
for i = 1:num_mics
    mic_data{i} = mic_data{i}(1:min_length);
end

% Motion Compensation: Adjust audio_data based on imu_data (Placeholder)
adjusted_audio_data = cell2mat(mic_data'); % Placeholder for motion compensation

% Initialize TDOA estimates
tdoa_estimates = zeros(size(mic_pairs, 1), 1);

% Loop through each microphone pair
for pair_idx = 1:size(mic_pairs, 1)
    mic_i = adjusted_audio_data(mic_pairs(pair_idx, 1), :);
    mic_j = adjusted_audio_data(mic_pairs(pair_idx, 2), :);
    
    % Apply STFT to both microphones
    win_len = round(window_size_ms * fs / 1000); % Window length in samples
    hop_size = round(win_len * (1 - overlap)); % Hop size for 50% overlap
    [S_i, f, t] = stft(mic_i, fs, 'Window', hamming(win_len, 'periodic'), 'OverlapLength', hop_size, 'FFTLength', win_len);
    [S_j, ~, ~] = stft(mic_j, fs, 'Window', hamming(win_len, 'periodic'), 'OverlapLength', hop_size, 'FFTLength', win_len);
    
    % Compute Cross-Power Spectrum
    G_ij = S_i .* conj(S_j);
    
    % Apply PHAT weighting
    G_phat = G_ij ./ (abs(G_ij) + eps); % Avoid division by zero
% Inverse STFT to get GCC-PHAT function
    [R_phat, tau] = istft(G_phat, fs, 'Window', hamming(win_len, 'periodic'), 'OverlapLength', hop_size, 'FFTLength', win_len);
    
    % Find peak to estimate TDOA
    [~, peak_idx] = max(abs(R_phat));
    tau_peak = tau(peak_idx); % TDOA in seconds
    
    % Estimate relative velocity (v_r) between drone and source
    % Placeholder for velocity estimation: v_r could be derived from IMU data analysis
    v_r = 0; % Placeholder (actual velocity estimation required)
    
    % Doppler Correction
    tau_corrected = tau_peak * (1 - v_r / speed_of_sound);
    
    % Store corrected TDOA estimate
    tdoa_estimates(pair_idx) = tau_corrected;
end

% Output TDOA estimates
disp('TDOA Estimates (corrected for Doppler effect):');
disp(tdoa_estimates);

% Triangulation or localization calculation using TDOA estimates and mic_coords
% Placeholder: Implement your localization algorithm using tdoa_estimates and mic_coords