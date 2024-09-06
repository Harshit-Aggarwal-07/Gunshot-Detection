#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <fftw3.h>
#include <sndfile.h>
#include <Eigen/Dense>

// GCC-PHAT for TDOA Determination
std::vector<double> gcc_phat(const std::vector<double>& x1, const std::vector<double>& x2, double fs) {
    size_t N = x1.size() + x2.size() - 1;
    std::vector<std::complex<double>> X1(N), X2(N), R(N), cc(N);

    fftw_complex* in1 = reinterpret_cast<fftw_complex*>(X1.data());
    fftw_complex* in2 = reinterpret_cast<fftw_complex*>(X2.data());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(R.data());

    fftw_plan plan1 = fftw_plan_dft_r2c_1d(N, x1.data(), in1, FFTW_ESTIMATE);
    fftw_plan plan2 = fftw_plan_dft_r2c_1d(N, x2.data(), in2, FFTW_ESTIMATE);
    fftw_plan ifft_plan = fftw_plan_dft_c2r_1d(N, out, cc.data(), FFTW_ESTIMATE);

    fftw_execute(plan1);
    fftw_execute(plan2);

    for (size_t i = 0; i < N; ++i) {
        R[i] = X1[i] * std::conj(X2[i]);
    }

    fftw_execute(ifft_plan);

    std::vector<double> correlation(cc.size());
    for (size_t i = 0; i < cc.size(); ++i) {
        correlation[i] = std::abs(cc[i]);
    }

    size_t maxIndex = std::distance(correlation.begin(), std::max_element(correlation.begin(), correlation.end()));
    double tdoa = static_cast<double>(maxIndex) / fs;

    return tdoa;
}

// Localization using TDOA
Eigen::Vector3d tdoa_localization(const std::vector<Eigen::Vector3d>& mic_positions, const std::vector<double>& tdoas, double speed_of_sound) {
    Eigen::MatrixXd A(mic_positions.size(), 3);
    Eigen::VectorXd b(mic_positions.size());

    for (size_t i = 0; i < mic_positions.size(); ++i) {
        Eigen::Vector3d pos = mic_positions[i];
        A.row(i) = (pos - mic_positions[0]).transpose();
        b(i) = tdoas[i] * speed_of_sound;
    }

    Eigen::Vector3d position = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    return position;
}

// Cuckoo Search Algorithm
std::vector<double> cuckoo_search(std::function<double(const std::vector<double>&)> objective_function,
                                  const std::vector<double>& initial_guess,
                                  const std::vector<std::pair<double, double>>& bounds,
                                  int num_nests, int num_iterations) {
    std::vector<std::vector<double>> nests(num_nests, initial_guess);
    std::vector<double> fitness(num_nests);

    // Compute initial fitness
    for (int i = 0; i < num_nests; ++i) {
        fitness[i] = objective_function(nests[i]);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> dis_int(0, num_nests - 1);

    for (int iter = 0; iter < num_iterations; ++iter) {
        std::vector<std::vector<double>> new_nests(num_nests, initial_guess);

        // Generate new nests
        for (int i = 0; i < num_nests; ++i) {
            int j = dis_int(gen);
            for (size_t k = 0; k < initial_guess.size(); ++k) {
                new_nests[i][k] = nests[i][k] + 0.01 * (nests[j][k] - nests[i][k]);
                // Apply bounds
                new_nests[i][k] = std::max(bounds[k].first, std::min(bounds[k].second, new_nests[i][k]));
            }
        }

        // Evaluate new nests
        std::vector<double> new_fitness(num_nests);
        for (int i = 0; i < num_nests; ++i) {
            new_fitness[i] = objective_function(new_nests[i]);
        }

        // Replace some nests
        for (int i = 0; i < num_nests; ++i) {
            if (dis(gen) < 0.25) {
                nests[i] = new_nests[dis_int(gen)];
                fitness[i] = new_fitness[i];
            }
        }
    }

    // Find best solution
    auto min_it = std::min_element(fitness.begin(), fitness.end());
    int min_idx = std::distance(fitness.begin(), min_it);
    return nests[min_idx];
}

// Calculate Azimuth and Elevation Angles
std::pair<double, double> calculate_angles(const Eigen::Vector3d& position) {
    double azimuth = std::atan2(position[1], position[0]) * 180.0 / M_PI;
    double elevation = std::atan2(position[2], std::sqrt(position[0] * position[0] + position[1] * position[1])) * 180.0 / M_PI;
    return {azimuth, elevation};
}

// Load audio from file
std::vector<double> load_audio(const std::string& filename, int& sample_rate) {
    SF_INFO sf_info;
    SNDFILE* file = sf_open(filename.c_str(), SFM_READ, &sf_info);
    if (!file) {
        throw std::runtime_error("Failed to open audio file.");
    }

    sample_rate = sf_info.samplerate;
    std::vector<double> audio(sf_info.frames);
    sf_readf_double(file, audio.data(), sf_info.frames);
    sf_close(file);

    return audio;
}

int main() {
    int fs;
    int num_mics = 5;
    std::vector<std::vector<double>> mic_signals(num_mics);
    std::vector<double> tdoas(num_mics - 1);
    std::vector<Eigen::Vector3d> mic_positions = {
        Eigen::Vector3d(0, 0, 0),
        Eigen::Vector3d(1, 0, 0),
        Eigen::Vector3d(0, 1, 0),
        Eigen::Vector3d(0, 0, 1),
        Eigen::Vector3d(1, 1, 0)
    };
    double speed_of_sound = 343;

    // Load the microphone signals
    for (int i = 0; i < num_mics; ++i) {
        mic_signals[i] = load_audio("mic" + std::to_string(i + 1) + ".wav", fs);
    }

    // Calculate TDOA
    for (int i = 1; i < num_mics; ++i) {
        tdoas[i - 1] = gcc_phat(mic_signals[i], mic_signals[0], fs);
    }

    // Perform localization
    Eigen::Vector3d estimated_position = tdoa_localization(mic_positions, tdoas, speed_of_sound);

    // Define the objective function for Cuckoo Search
    auto objective_function = [&](const std::vector<double>& pos) -> double {
        Eigen::Vector3d test_pos(pos[0], pos[1], pos[2]);
        Eigen::Vector3d estimated = tdoa_localization(mic_positions, tdoas, speed_of_sound);
        return (test_pos - estimated).norm();
    };

    // Cuckoo Search
    std::vector<std::pair<double, double>> bounds = { {0, 10}, {0, 10}, {0, 10} };
    std::vector<double> initial_guess = { 1, 1, 1 };
    int num_nests = 50;
    int num_iterations = 100;
    std::vector<double> best_position = cuckoo_search(objective_function, initial_guess, bounds, num_nests, num_iterations);

    // Convert best_position to Eigen::Vector3d
    Eigen::Vector3d best_position_vec(best_position[0], best_position[1], best_position[2]);

    // Calculate angles
    auto [azimuth, elevation] = calculate_angles(best_position_vec);
    std::cout << "Azimuth: " << azimuth << " degrees" << std::endl;
    std::cout << "Elevation: " << elevation << " degrees" << std::endl;

    return 0;
}