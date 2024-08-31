% Fast RLS Adaptive Noise Cancellation in MATLAB

% Parameters
L = 64;          % Filter length
LAMBDA = 0.99;   % Forgetting factor
DELTA = 0.01;    % Regularization parameter

% Load audio signals
[primarySignal, fs1] = audioread('NewCompile/gunshot_delayed_M1.wav');
[referenceNoise, fs2] = audioread('NewCompile/gunshot_delayed_M6.wav');

% Ensure the sample rates are the same
if fs1 ~= fs2
    error('Sample rates of the two audio files do not match.');
end

% Make sure the length of both input clips is the same by padding the shorter one
N1 = length(primarySignal);
N2 = length(referenceNoise);

if N1 > N2
    referenceNoise = [referenceNoise; zeros(N1 - N2, 1)];
elseif N2 > N1
    primarySignal = [primarySignal; zeros(N2 - N1, 1)];
end

N = max(N1, N2);

% Allocate memory for the output signal
outputSignal = zeros(N, 1);

% Initialize filter coefficients and inverse correlation matrix
W = zeros(L, 1);        % Filter coefficients
P = DELTA * eye(L);     % Inverse correlation matrix
xBuffer = zeros(L, 1);  % Circular buffer for reference noise

% Process the signal sample by sample
lambdaInverse = 1 / LAMBDA;

for n = L:N
    % Update circular buffer
    xBuffer(2:end) = xBuffer(1:end-1);
    xBuffer(1) = referenceNoise(n);

    % Compute the output of the adaptive filter
    y = W' * xBuffer;

    % Calculate the error signal
    outputSignal(n) = primarySignal(n) - y;

    % Calculate the Kalman gain vector
    K = P * xBuffer / (lambdaInverse + xBuffer' * P * xBuffer);

    % Update the inverse correlation matrix P
    P = lambdaInverse * (P - K * xBuffer' * P);

    % Update the filter coefficients
    W = W + K * outputSignal(n);
end

% Save the filtered output
audiowrite('NewCompile/filtered_output_fast_rls_2.wav', outputSignal,fs1);