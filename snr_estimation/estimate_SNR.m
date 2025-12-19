function [prio, post] = estimate_SNR(Y, tinc, snr_smoothing_alpha, Estimator)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATE_SNR estimates prior and posterior SNR with default parameters
%
% This is a convenience wrapper function for snr_decision_directed that
% provides sensible default parameters for typical speech enhancement
% applications. It uses MMSE noise estimation and configurable gain functions.
%
% Usage:
%   [prio, post] = estimate_SNR(Y, tinc, snr_smoothing_alpha, Estimator);
%
% Inputs:
%   Y                    Complex STFT of noisy speech (frames x frequency bins)
%   tinc                 Frame increment in seconds (used for noise estimation)
%   snr_smoothing_alpha  Smoothing factor for prior SNR (0 < alpha < 1)
%                        Typical value: 0.98. Higher values = smoother estimates
%   Estimator            String specifying gain function estimator (see STFT_GetGain)
%                        Options: 'WienerFilter', 'MMSESTSA', 'LSA', 'JMAPLV',
%                                 'JMAPGW', 'MAPGW', 'MMSEGW'
%
% Outputs:
%   prio                 Prior SNR estimate (frames x frequency bins)
%   post                 Posterior SNR estimate (frames x frequency bins)
%
% Default Parameters:
%   p.ne  = 1            Use MMSE noise estimation (Gerkmann & Hendriks method)
%   p.ri  = 1            Round frame increment to nearest power of 2 samples
%   p.gx  = 10000000     Maximum posterior SNR (70 dB)
%   p.gn  = 1            Min posterior SNR when estimating prior SNR (0 dB)
%   p.gz  = 0.0001       Min posterior SNR (-40 dB floor)
%   p.xn  = 0            Minimum prior SNR (0, equivalent to -Inf dB)
%   p.xb  = 1            Bias compensation factor
%
% The default parameters are chosen to:
%   - Use robust MMSE noise tracking
%   - Prevent numerical issues with very high/low SNR values
%   - Provide unbiased prior SNR estimates
%
% For more control over parameters, call snr_decision_directed directly.
%
% See also: snr_decision_directed, STFT_GetGain, estnoiseg
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
    p.ne = 1;
    p.ri = 1;        
    p.gx = 10000000;
    p.gn = 1;        
    p.gz = 0.0001;
    p.xn = 0;
    p.xb = 1;
    p.Estimator=Estimator;
    
    [prio, post] = snr_decision_directed(Y, tinc, snr_smoothing_alpha, p);
end