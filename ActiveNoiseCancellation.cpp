#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define L 64  // Filter length
#define LAMBDA 0.99  // Forgetting factor
#define DELTA 0.01  // Regularization parameter

// Function to read WAV file header and data
void loadWavFile(const char *filename, int16_t **signal, int *length, int *sampleRate) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    // Read WAV header
    fseek(file, 24, SEEK_SET);
    fread(sampleRate, 4, 1, file); // Sample rate
    fseek(file, 40, SEEK_SET);
    int dataSize;
    fread(&dataSize, 4, 1, file); // Data size
    *length = dataSize / sizeof(int16_t);

    fseek(file, 44, SEEK_SET); // Move to start of data
    *signal = (int16_t *)malloc(*length * sizeof(int16_t));
    fread(*signal, sizeof(int16_t), *length, file);

    fclose(file);
}

// Function to write WAV file header
void writeWavHeader(FILE *file, int sampleRate, int numChannels, int bitsPerSample, int numSamples) {
    int byteRate = sampleRate * numChannels * bitsPerSample / 8;
    int blockAlign = numChannels * bitsPerSample / 8;
    int dataSize = numSamples * numChannels * bitsPerSample / 8;

    // Write RIFF header
    fwrite("RIFF", 1, 4, file); // Chunk ID
    uint32_t chunkSize = 36 + dataSize;
    fwrite(&chunkSize, 4, 1, file); // Chunk Size
    fwrite("WAVE", 1, 4, file); // Format

    // Write fmt subchunk
    fwrite("fmt ", 1, 4, file); // Subchunk 1 ID
    uint32_t subchunk1Size = 16;
    fwrite(&subchunk1Size, 4, 1, file); // Subchunk 1 Size
    uint16_t audioFormat = 1;
    fwrite(&audioFormat, 2, 1, file); // Audio Format (PCM)
    fwrite(&numChannels, 2, 1, file); // Number of Channels
    fwrite(&sampleRate, 4, 1, file); // Sample Rate
    fwrite(&byteRate, 4, 1, file); // Byte Rate
    fwrite(&blockAlign, 2, 1, file); // Block Align
    fwrite(&bitsPerSample, 2, 1, file); // Bits per Sample

    // Write data subchunk
    fwrite("data", 1, 4, file); // Subchunk 2 ID
    fwrite(&dataSize, 4, 1, file); // Subchunk 2 Size
}

// Function to save the processed audio data
void saveWavFile(const char *filename, int16_t *signal, int length, int sampleRate) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    int numChannels = 1; // Mono
    int bitsPerSample = 16; // 16-bit samples

    // Write WAV header
    writeWavHeader(file, sampleRate, numChannels, bitsPerSample, length);

    // Write audio data
    fwrite(signal, sizeof(int16_t), length, file);

    fclose(file);
}

// Main Fast RLS ANC function
void fastRLS(int16_t *primarySignal, int16_t *referenceNoise, int16_t *outputSignal, int N) {
    float W[L] = {0};  // Filter coefficients
    float P[L][L] = {0};  // Inverse correlation matrix
    float xBuffer[L] = {0};  // Circular buffer for reference noise
    float K[L] = {0};  // Kalman gain vector
    float lambdaInverse = 1.0f / LAMBDA;

    // Initialize the inverse correlation matrix P
    for (int i = 0; i < L; i++) {
        P[i][i] = DELTA;
    }

    // Process the signal sample by sample
    for (int n = L; n < N; n++) {
        // Update circular buffer
        for (int i = L-1; i > 0; i--) {
            xBuffer[i] = xBuffer[i-1];
        }
        xBuffer[0] = referenceNoise[n];

        // Compute the output of the adaptive filter
        float y = 0.0f;
        for (int i = 0; i < L; i++) {
            y += W[i] * xBuffer[i];
        }

        // Calculate the error signal
        outputSignal[n] = (int16_t)(primarySignal[n] - y);

        // Calculate the Kalman gain vector
        float xPx = 0.0f;
        for (int i = 0; i < L; i++) {
            K[i] = 0.0f;  // Initialize Kalman gain
            for (int j = 0; j < L; j++) {
                K[i] += P[i][j] * xBuffer[j];
            }
            xPx += K[i] * xBuffer[i];
        }
        float denominator = lambdaInverse + xPx;
        for (int i = 0; i < L; i++) {
            K[i] /= denominator;
        }

        // Update the inverse correlation matrix P
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                P[i][j] = lambdaInverse * (P[i][j] - K[i] * xBuffer[j] * P[j][i]);
            }
        }

        // Update the filter coefficients
        for (int i = 0; i < L; i++) {
            W[i] += K[i] * (primarySignal[n] - y);
        }
    }
}

int main() {
    int16_t *primarySignal, *referenceNoise, *outputSignal;
    int N_primary, N_reference;
    int sampleRate;

    // Load the audio signals
    loadWavFile("C:/Users/Hp/Documents/Audacity/primary_input.wav", &primarySignal, &N_primary, &sampleRate);
    loadWavFile("C:/Users/Hp/Documents/Audacity/reference_noise.wav", &referenceNoise, &N_reference, &sampleRate);

    // Determine the maximum length and pad the shorter signal
    int maxLength = (N_primary > N_reference) ? N_primary : N_reference;
    if (N_primary < maxLength) {
        primarySignal = (int16_t *)realloc(primarySignal, maxLength * sizeof(int16_t));
        for (int i = N_primary; i < maxLength; i++) {
            primarySignal[i] = 0;  // Padding with zeros
        }
        N_primary = maxLength;
    }
    if (N_reference < maxLength) {
        referenceNoise = (int16_t *)realloc(referenceNoise, maxLength * sizeof(int16_t));
        for (int i = N_reference; i < maxLength; i++) {
            referenceNoise[i] = 0;  // Padding with zeros
        }
        N_reference = maxLength;
    }

    // Allocate memory for the error signal
    outputSignal = (int16_t *)malloc(maxLength * sizeof(int16_t));

    // Run the Fast RLS ANC algorithm
    fastRLS(primarySignal, referenceNoise, outputSignal, maxLength);

    // Save the filtered output
    saveWavFile("C:/Users/Hp/Documents/Audacity/filtered_output_fast_rls.wav", outputSignal, maxLength, sampleRate);

    // Free allocated memory
    free(primarySignal);
    free(referenceNoise);
    free(outputSignal);

    return 0;
}