% Load the audio file
[audioData, fs] = audioread('gunshotsample.wav'); % Replace 'gunshot.wav' with your actual file name

% Time Difference of Arrival (TDOA) relative to M1 in seconds
tdoaRelativeToM1 = [0, 0.00034, 0.0259, 0.0086, 0.0281, 0.1076]; 

% Number of microphones
numMics = length(tdoaRelativeToM1);

% Pre-allocate space for delayed audio signals
maxDelaySamples = round(max(tdoaRelativeToM1) * fs); % Maximum delay in samples
delayedAudio = zeros(length(audioData) + maxDelaySamples, numMics); 

% Apply delays to each microphone
for i = 1:numMics
    delaySamples = round(tdoaRelativeToM1(i) * fs); % Convert delay to samples
    endIndex = delaySamples + length(audioData);
    
    % Ensure endIndex does not exceed the size of delayedAudio
    if endIndex > size(delayedAudio, 1)
        endIndex = size(delayedAudio, 1);
    end
    
    % Assign the audio data with delay
    delayedAudio(delaySamples + 1:endIndex, i) = audioData(1:(endIndex - delaySamples));
end

% Normalize the delayed audio signals to prevent clipping
for i = 1:numMics
    delayedAudio(:, i) = delayedAudio(:, i) / max(abs(delayedAudio(:, i)));
end

% Save the delayed audio for each microphone to separate files
for i = 1:numMics
    filename = sprintf('gunshot_delayed_M%d.wav', i);
    audiowrite(filename, delayedAudio(:, i), fs);
end

% Play the delayed sounds for each microphone (optional)
for i = 1:numMics
    fprintf('Playing delayed sound for Microphone M%d\n', i);
    sound(delayedAudio(:, i), fs);
    pause((length(delayedAudio) / fs) + 1); % Wait until the sound finishes
end


