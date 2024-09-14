% Updated test.m for multiple microphones
load('best_svm_gunshot_model.mat', 'bestModel');

for mic_no = 1:5
    % Load events detected for the mic
    load(sprintf('events_detected_mic%d.mat', mic_no), 'events_detected');
    
    % If no event detected, skip to next mic
    if ~any(events_detected)
        fprintf('No events detected in mic%d, skipping gunshot detection.\n', mic_no);
        continue;
    end

    % Load the audio file for gunshot detection if any event detected
    [sampleAudio, sr] = audioread(sprintf('digital_mic%d.wav', mic_no));
    windowSamples = floor(0.5 * sr);  % 0.5 seconds window
    stepSamples = floor(0.25 * sr);  % 0.25 seconds overlap
    lpcOrder = 8;  % LPC order

    % Feature extraction and SVM prediction
    [gunshotTemplate, sr] = audioread('positive gunshots/1 (1).wav');
    features = [];
    for j = 1:stepSamples:(length(sampleAudio) - windowSamples + 1)
        window = sampleAudio(j:j + windowSamples - 1);
        feat = extract_features(window, gunshotTemplate, lpcOrder);
        features = [features; feat];
    end

    % Predict gunshot events
    predictions = predict(bestModel, features);

    if any(predictions == 1)
        fprintf('Gunshot detected in mic%d audio.\n', mic_no);
        gunshot_detected = true;
    else
        fprintf('No gunshot detected in mic%d.\n', mic_no);
    end
end

function feat = extract_features(segment, gunshotTemplate, lpcOrder)
    % Ensure the segment is mono
    if size(segment, 2) > 1
        segment = mean(segment, 2);  % Convert to mono if necessary
    end
    
    % Ensure the gunshot template is mono
    if size(gunshotTemplate, 2) > 1
        gunshotTemplate = mean(gunshotTemplate, 2);  % Convert to mono if necessary
    end
    
    % Adjust the length of the segment to match the gunshotTemplate
    if length(segment) < length(gunshotTemplate)
        segment = [segment; zeros(length(gunshotTemplate) - length(segment), 1)];
    elseif length(segment) > length(gunshotTemplate)
        segment = segment(1:length(gunshotTemplate));
    end
    
    % Cross-correlation with the gunshot template
    correlation = xcorr(segment, gunshotTemplate, 'coeff');
    crossCorrelationMax = max(correlation);  % Max value of cross-correlation
    
    % Extract LPC coefficients
    lpcCoeffs = lpc(segment, lpcOrder);  % LPC analysis
    lpcCoeffs = lpcCoeffs(:)';  % Ensure row vector
    
    % Combine features into a single vector
    feat = [crossCorrelationMax, lpcCoeffs];  % Combine into one feature set
end


