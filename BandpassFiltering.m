%BANDPASS FILTERING

[audio_signal, fs] = audioread('Compile/combined.wav');

% Define the bandpass filter parameters
low_cutoff = 250;  % Lower cutoff frequency (2.5 kHz)
high_cutoff = 3100; % Upper cutoff frequency (3.5 kHz)
filter_order = 600;  % Order of the FIR filter

% Design the bandpass filter
bp_filter = designfilt('bandpassfir', 'FilterOrder', filter_order, ...
                       'CutoffFrequency1', low_cutoff, 'CutoffFrequency2', high_cutoff, ...
                       'SampleRate', fs);

% Apply the bandpass filter to the audio signal
filtered_signal = filter(bp_filter, audio_signal);

%amplitude plot
figure;
subplot(2, 1, 1);
plot(audio_signal(:, 1));
title('Original Signal');
xlabel('Sample Number');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(filtered_signal(:, 1));
title('Bandpass Filtered Signal (Mic 1)');
xlabel('Sample Number');
ylabel('Amplitude');

% Plot the original and filtered signals in the frequency domain

% Calculate the FFT of the original and filtered signals
N = length(audio_signal);
f = (0:N-1)*(fs/N);

% FFT of the original signal
input_fft = abs(fft(audio_signal));

% FFT of the filtered signal
output_fft = abs(fft(filtered_signal));

% Plot the Frequency Spectrum of the Original Signal
figure;
subplot(2, 1, 1);
plot(f(1:round(N/2)), input_fft(1:round(N/2)));
title('Frequency Spectrum of Original Audio Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 5000]); % Limit the x-axis to the 100-5 kHz range

% Plot the Frequency Spectrum of the Filtered Signal
subplot(2, 1, 2);
plot(f(1:round(N/2)), output_fft(1:round(N/2)));
title('Frequency Spectrum of Filtered Audio Signal (2.5k - 3.5k Hz)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 5000]); % Limit the x-axis to the 100-5 kHz range

% Save the filtered signal as a new audio file 

audiowrite('Compile/filtered_primary_audio.wav', filtered_signal,fs);