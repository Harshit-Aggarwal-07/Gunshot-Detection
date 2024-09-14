% Updated acoustic_event_detection.m for multiple microphones
for mic_no = 1:5
    % Load the overlapping windows and sample rate from the .mat file
    load(sprintf('overlapping_windows_mic%d.mat', mic_no), 'windows', 'sr');

    % Parameters for event detection
    N0 = 1e-4;  % Example noise power spectral density
    T = 1 / sr;  % Sampling period
    W = sr / 2;  % Positive bandwidth (Nyquist frequency)

    % Define the chi2inv function in MATLAB to calculate the inverse chi-square
    calculate_adaptive_threshold = @(P_fa, k) chi2inv(1 - P_fa, k);

    % Function to calculate the decision statistic Z
    calculate_decision_statistic = @(r, N0, T, W) (2 / N0) * sum(r.^2) * T;

    % Apply the acoustic event detection on each overlapping window
    events_detected = false(1, size(windows, 2));  % Pre-allocate for speed

    for i = 1:size(windows, 2)
        window = windows(:, i);

        % Calculate adaptive threshold for this window
        adaptive_threshold = calculate_adaptive_threshold(0.01, 2 * W * T);

        % Calculate decision statistic Z
        Z = calculate_decision_statistic(window, N0, T, W);

        % Decision rule
        events_detected(i) = Z > adaptive_threshold;
    end

    % Save the detected events to a .mat file
    save_file_name = sprintf('events_detected_mic%d.mat', mic_no);
    save(save_file_name, 'events_detected');

    fprintf('Events detected in mic%d: ', mic_no);
    disp(events_detected);
end
