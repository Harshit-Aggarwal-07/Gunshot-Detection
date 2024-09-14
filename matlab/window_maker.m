% Updated window_maker.m to process multiple files (mic1 to mic5)
for mic_no = 1:5
    % Load the audio file for the given microphone
    file_name = sprintf('Anc_output_mic%d.wav', mic_no);
    [audio_signal, sr] = audioread(file_name);

    % Define parameters
    window_duration = 0.5;  % Each window duration in seconds
    overlap_duration = 0.25;  % 50% overlap duration in seconds
    window_samples = floor(sr * window_duration);  % Samples per window
    overlap_samples = floor(sr * overlap_duration);  % Samples per overlap
    step_size = window_samples - overlap_samples;  % Step size between windows

    % Create overlapping windows
    num_windows = floor((length(audio_signal) - window_samples) / step_size) + 1;
    windows = zeros(window_samples, num_windows);

    for i = 1:num_windows
        start_idx = (i-1) * step_size + 1;
        windows(:, i) = audio_signal(start_idx:start_idx + window_samples - 1);
    end

    % Save the windows to a .mat file
    save_file_name = sprintf('overlapping_windows_mic%d.mat', mic_no);
    save(save_file_name, 'windows', 'sr');

    fprintf('Number of overlapping windows for mic%d: %d\n', mic_no, num_windows);
end
