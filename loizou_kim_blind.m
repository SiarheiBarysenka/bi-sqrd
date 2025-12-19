function result = loizou_kim_blind(y, p)
%LOIZOU_KIM_BLIND Performs speech enhancement using Loizou-Kim mask (SNR_prio > 1/3)
% in blind setting
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
%          p.Estimator              :  defines gain function to use
%                                      when estimating the SNR.
%                                      See STFT_GetGain function for
%                                      possible cases.
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
%          result.SNR_prio          :  a priori SNR estimate
%          
%          result.mask              :  binary mask by Loizou-Kim method (no phase awareness)
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
    [prio, ~] = estimate_SNR(Y, tinc, p.SNR_smoothing_alpha, p.Estimator);

    %% Retain STFT and SNR columns at f = 0 and remove them from transformations
    [S0, S] = split_f0_column(Y);
    [prio0, prio] = split_f0_column(prio);

    %% Loizou-Kim mask
    mask_magnitude = loizou_kim_mask(prio);

    %% Transform the enhanced signal back to time domain
    [xhat, Xhat] = reconstruct_from_mask(mask_magnitude, S0, S, Nfft, anaWin, tinc, p.fs, y);

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
    
    result.SNR_prio = [prio0 prio];

    result.mask = mask_magnitude;

    result.RTF = RTF;
end