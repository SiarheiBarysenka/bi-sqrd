# Bi²: Real-Time Single-Channel Speech Enhancement in Bispectral Domain Using Two Slices of Bispectral Phase

## Listening Examples

Audio samples comparing Bi² with baseline methods are available at: https://siarheibarysenka.github.io/bi-sqrd/

## Usage

```matlab
%% Set parameters
p.fs                    = 16000;                    % Sampling frequency (Hz)
p.win_len               = 512;                      % Analysis window length
p.Nfft_STFT             = 512;                      % FFT size for STFT
p.overlap_factor_STFT   = 7 / 8;                    % Overlap factor for windowing
p.SNR_smoothing_alpha   = 0.98;                     % Alpha for SNR estimation
p.win                   = @(L) tukeywin(L, 0.15);   % Window type
p.Estimator             = 'LSA';                    % Gain function for SNR estimation
p.kappa_window_length   = p.win_len / p.fs / 4;     % Kappa estimation window (s)
p.p_H0_Theta2           = 0.60;                     % Null hypothesis prob. (diagonal slice)
p.p_H0_Theta3           = 0.75;                     % Null hypothesis prob. (next-to-diag. slice)

%% Bi² speech enhancement
result = bi_sqrd(noisy_speech, p);
bi_sqrd_enhanced_speech = result.xhat;

%% Loizou-Kim baseline (blind)
result_lk = loizou_kim_blind(noisy_speech, p);
lk_enhanced_speech = result_lk.xhat;

%% Loizou-Kim baseline (oracle)
result_lk_oracle = loizou_kim_oracle(noisy_speech, clean_speech, noise, p);
lk_oracle_enhanced_speech = result_lk_oracle.xhat;
```

## References

[1] S.Y. Barysenka and P. Mowlaee, “Bi²: Phase-Aware Binary Mask in Bispectral Domain for Single-Channel Speech Enhancement,” Submitted to IEEE Transactions on Audio, Speech and Language Processing.

[2] P. C. Loizou and G. Kim, “Reasons why current speech-enhancement algorithms do not improve speech intelligibility and suggested solutions,” IEEE Transactions on Audio, Speech and Language Processing, vol. 19, no. 1, pp. 47–56, Jan 2011. [DOI](https://doi.org/10.1109/TASL.2010.2045180)
