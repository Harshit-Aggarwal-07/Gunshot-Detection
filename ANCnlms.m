% Load Primary and Reference Audio Signals
[primarySignal, fsPrimary] = audioread('Compile/filtered_primary_audio.wav');
[referenceNoise, fsReference] = audioread('Compile/filtered_reference_audio.wav');

% Ensure both signals have the same sampling rate
if fsPrimary ~= fsReference
    error('Sampling rates of the primary and reference signals must be the same.');
end

% Ensure both signals have the same length
minLength = min(length(primarySignal), length(referenceNoise));
primarySignal = primarySignal(1:minLength);
referenceNoise = referenceNoise(1:minLength);

% Parameters for NLMS Algorithm
L = 64;                % Filter length
mu = 0.05;             % Step size for NLMS
epsilon = 1e-6;        % Small constant to prevent division by zero

% Initialize NLMS Variables
W = zeros(L, 1);       % Filter coefficients (initially zeros)
xBuffer = zeros(L, 1); % Buffer for the reference noise
outputSignal = zeros(minLength, 1); % Output signal (error signal)

% NLMS Adaptive Filtering
for n = L:minLength
    % Update buffer with the current reference noise sample
    xBuffer = [referenceNoise(n); xBuffer(1:end-1)];
    
    % Filter output (predicted noise)
    y = W' * xBuffer;
    
    % Error signal (desired signal - predicted noise)
    e = primarySignal(n) - y;
    outputSignal(n) = e;
    
    % Update filter coefficients using NLMS update rule
    normFactor = (xBuffer' * xBuffer) + epsilon; % Normalization factor
    W = W + (mu / normFactor) * e * xBuffer;
end

% Save the filtered output to a WAV file
audiowrite('Compile/filtered_signal_nlms.wav', outputSignal, fsPrimary);

% Plot the Original and Filtered Signals
figure;
subplot(3, 1, 1);
plot(primarySignal);
title('Original Primary Signal with Noise');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(referenceNoise);
title('Reference Noise Signal');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(outputSignal);
title('Filtered Signal (Output of NLMS ANC)');
xlabel('Sample Index');
ylabel('Amplitude');
