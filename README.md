# Signal Filtering and Analysis Project

## Description

This MATLAB project demonstrates signal processing techniques for analyzing and filtering noisy signals using Butterworth and Chebyshev filters. The project includes:
- Time-domain and frequency-domain signal analysis
- Design and implementation of analog filters
- Comparison of different filter types and orders
- Visualization of original, noisy, and filtered signals

## Key Features

- **Signal Analysis**:
  - Time-domain visualization of pure and noisy signals
  - Frequency-domain analysis using FFT
- **Filter Design**:
  - Butterworth filters (low-pass and high-pass)
  - Chebyshev Type I filters (low-pass and high-pass)
  - Band-pass filter creation through cascade combinations
- **Performance Evaluation**:
  - Comparison of filtered vs original signals
  - Frequency response analysis of filters
- **Data Management**:
  - Option to save processed signals for further analysis

## Getting Started

### Prerequisites

- MATLAB (with Signal Processing Toolbox)
- Input data file: `time_pure_noisy.txt`

### Installation

1. Place all files in the same directory
2. Ensure `time_pure_noisy.txt` is in the working directory
3. Run `hadi_khaled_mohamed_sec2_bn28.m`

## Usage

The script automatically processes the input data and generates multiple figures:

1. **Time Domain Plots**:
   - Original pure signal vs noisy signal
   - Filtered signals vs original for each filter type

2. **Frequency Domain Plots**:
   - Power spectrum of original and noisy signals
   - Comparison of filtered vs noisy signals in frequency domain

3. **Filter Characteristics**:
   - Frequency response of each designed filter

### Customization Options

1. Modify filter parameters:
   - `n_butter_l_1`, `fc_butter_l_1`: Order and cutoff for Butterworth low-pass
   - `n_cheb_h_1`, `fc_cheb_h_1`: Order and cutoff for Chebyshev high-pass
   - `rp`: Passband ripple for Chebyshev filters

2. Toggle data saving:
   ```matlab
   save_data = true; % or false to disable saving
   ```

## File Descriptions

- `hadi_khaled_mohamed_sec2_bn28.m`: Main analysis script
- `time_pure_noisy.txt`: Input data file (time, pure signal, noisy signal)
- `all signals.txt`: Output file containing processed signals (if saved)

## Technical Details

### Filter Design Approach

1. **Butterworth Filters**:
   - 8th order low-pass (20×2π rad/s)
   - 8th order high-pass (10×2π rad/s)
   - Combined to create band-pass characteristics

2. **Chebyshev Type I Filters**:
   - 5th/6th order designs
   - 3dB passband ripple
   - Cutoff frequencies optimized for noise reduction

### Signal Processing

- Sampling frequency calculated from time vector
- Proper FFT normalization for accurate spectral representation
- Single-sided spectrum computation for visualization

## Results Interpretation

- **Time Domain**: Evaluate how well filters reconstruct the original signal
- **Frequency Domain**: Assess noise reduction in different frequency bands
- **Filter Responses**: Compare roll-off characteristics and passband behavior

## License

This project is available for academic and research purposes. For commercial use, please contact the author.

## Acknowledgments

This implementation demonstrates fundamental signal processing concepts using MATLAB's Signal Processing Toolbox functions.
