% Main file to handle the complete process
num_mics = 6;

for i = 1:num_mics  % Loop through different values
    inputIndex = i; % Set the input value for the other script
    run('bandpass.m'); % Run the other script
end


for i = 1:num_mics-1  % Loop through different values
    inputIndex = i; % Set the input value for the other script
    run('anc_rls.m'); % Run the other script
end


% 1. Window creation for each mic signal
window_maker;

% 2. Acoustic event detection
acoustic_event_detection;

% 3. Gunshot detection
gunshot_detected = false;  % Initialize flag

test;

% 4. If gunshot detected, proceed with localization and GCC-PHAT
if gunshot_detected
    gcc_phat_and_localization;
    % 5. Apply beamforming and MUSIC algorithm for enhanced localization
    beamforming_and_music;
else
    disp('No gunshot detected in any microphone signal.');
end
