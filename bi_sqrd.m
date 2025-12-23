function result = bi_sqrd(y, p)
%BI_SQRD Performs speech enhancement using phase-aware binary mask in bispectral domain.
%
%   Input variables: 
%
%          y                        :  input signal in time domain
%          p.fs                     :  sampling frequency (in hertz)
%          p.win_len                :  analysis window length
%          p.Nfft_STFT              :  FFT size for STFT
%                                      (if greater than p.win_len,
%                                       zero padding will be applied) 
%          p.overlap_factor_STFT    :  overlap factor (< 1) for windowing
%                                      used for STFT calculation
%          p.SNR_smoothing_alpha    :  alpha for SNR estimation
%                                      (typically, 0.96 ... 0.98)
%          p.win                    :  window type (e.g. @blackman):
%                                      if [], default will be
%                                      sqrt hamming window
%          p.kappa_window_length    :  defines the sliding window length
%                                      (in seconds) for collecting the
%                                      circular statistics for kappa
%                                      estimation.
%          p.Estimator              :  defines gain function to use
%                                      when estimating the SNR.
%                                      See STFT_GetGain function for
%                                      possible cases.
%          p.p_H0_Theta2            :  null hypothesis probability for 
%                                      diagonal slice of bi-phase
%                                      (typically 0.6)
%          p.p_H0_Theta3            :  null hypothesis probability for 
%                                      next-to-diagonal slice of bi-phase 
%                                      (typically 0.75)
%
%   Output variables:
%                  
%          result.xhat              :  enhanced signal in time domain
%
%          result.Y                 :  STFT of input signal y
%          result.Xhat              :  STFT of enhanced signal xhat
%          result.STFT_time         :  time axis of STFT
%          result.STFT_frequency    :  frequency axis of STFT
%          
%          result.B2                :  diagonal slice of bi-spectrum (complex-valued) of input signal
%          result.B3                :  next-to-diagonal slice of bi-spectrum (complex-valued) of input signal
%          
%          result.SNR_prio          :  a priori SNR estimate
%          result.SNR_post          :  a posteriori SNR estimate
%          
%          result.mask_B2           :  binary mask of the diagonal slice of bi-phase
%          result.mask_B3           :  binary mask of the next-to-diagonal slice of bi-phase
%          result.mask              :  resulting binary mask
%          
%          result.RTF               :  real-time factor
%          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Bi^2: Phase-Aware Binary Mask in Bispectral Domain
%          for Single-Channel Speech Enhancement
%
%    Copyright (C) 2026  Siarhei Y. Barysenka
%    
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%    
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%    
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Start processing time measurement
    tic

    %% Compute STFT
    [Y, STFT_time, STFT_frequency, Nfft, tinc, anaWin, ~] = calculate_STFT(y, p.fs, p.win_len, p.Nfft_STFT, p.overlap_factor_STFT, p.win);
    
    %% Estimate SNR
    [prio, post] = estimate_SNR(Y, tinc, p.SNR_smoothing_alpha, p.Estimator);
    
    %% Retain STFT and SNR columns at f = 0 and remove them from bi-spectrum and statistics computation
    [S0, S] = split_f0_column(Y);
    [prio0, prio] = split_f0_column(prio);
    [post0, post] = split_f0_column(post);
    
    %% Perform bi-spectral processing
    [B2, ~, angle_Theta2, k_1_B2, k_2_B2, k_3_B2] = bispectrum_diagonal_slice(S);
    [B3, ~, angle_Theta3, k_1_B3, k_2_B3, k_3_B3] = bispectrum_next_to_diagonal_slice(S);

    %% Make non-even number of window length for kappa calculation
    Nwin_kappa = fix(p.kappa_window_length / tinc / 2) * 2 + 1;

    %% Binary mask of the diagonal slice of bi-phase
    [mask_B2, ~, ~, ~] = bi_phase_slice_binary_mask(angle_Theta2, k_1_B2, k_2_B2, k_3_B2, Nwin_kappa, prio, post, p.p_H0_Theta2);

    %% Binary mask of the next-to-diagonal slice of bi-phase
    [mask_B3, ~, ~, ~] = bi_phase_slice_binary_mask(angle_Theta3, k_1_B3, k_2_B3, k_3_B3, Nwin_kappa, prio, post, p.p_H0_Theta3);

    %% Bi^2 mask
    mask_magnitude = loizou_kim_mask(prio);
    mask_bi_sqrd = bi_sqrd_mask(mask_magnitude, mask_B2', mask_B3');

    %% Transform the enhanced signal back to time domain
    [xhat, Xhat] = reconstruct_from_mask(mask_bi_sqrd, S0, S, Nfft, anaWin, tinc, p.fs, y);

    %% Stop processing time measurement
    t_proc = toc;

    %% Compute RTF
    t_sig = length(y) / p.fs;
    RTF = t_proc / t_sig;

    %% Assign output variables
    result.xhat = xhat;
    result.Y = Y;
    result.Xhat = Xhat;
    result.STFT_time = STFT_time;
    result.STFT_frequency = STFT_frequency;
    
    result.B2 = B2;
    result.B3 = B3;
    
    result.SNR_prio = [prio0 prio]; % With restored column at f = 0
    result.SNR_post = [post0 post]; % With restored column at f = 0

    result.mask_B2 = mask_B2;
    result.mask_B3 = mask_B3;
    result.mask = mask_bi_sqrd;

    result.RTF = RTF;
end