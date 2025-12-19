function [xi_mat,gam]=snr_decision_directed(Y, tinc, snr_smoothing_alpha, pp)
%SNR_DECISION_DIRECTED estimates prior and posterior SNR using decision-directed approach
%
% This function implements the decision-directed SNR estimation algorithm for
% speech enhancement. It is adapted from the VOICEBOX ssubmmse function but
% focuses solely on SNR estimation with configurable gain functions.
%
% Usage: 
%   [xi_mat, gam] = snr_decision_directed(Y, tinc, snr_smoothing_alpha, pp);
%
% Inputs:
%   Y                    Complex STFT of noisy speech (frames x frequency bins)
%   tinc                 Frame increment in seconds (used for noise estimation)
%   snr_smoothing_alpha  Smoothing factor for prior SNR estimation (0 < alpha < 1)
%                        Typical value: 0.98. Higher values = more smoothing
%   pp                   Parameter structure with fields:
%        pp.gx           Maximum posterior SNR as power ratio [default: 1000 = +30dB]
%        pp.gz           Minimum posterior SNR as power ratio [default: 0.001 = -30dB]
%        pp.gn           Min posterior SNR ratio when estimating prior SNR [default: 1 = 0dB]
%        pp.xn           Minimum prior SNR [default: 0]
%        pp.xb           Bias compensation factor for prior SNR [default: 1]
%        pp.ne           Noise estimation method: 0=minimum statistics [3], 1=MMSE [7]
%        pp.Estimator    Gain function estimator type (passed to STFT_GetGain)
%
% Outputs:
%   xi_mat              Prior SNR estimate (frames x frequency bins)
%                       xi = E[|S(f)|^2] / E[|N(f)|^2] where S=clean, N=noise
%   gam                 Posterior SNR estimate (frames x frequency bins)
%                       gamma = |Y(f)|^2 / E[|N(f)|^2] where Y=noisy
%
% Algorithm:
%   1. Compute power spectrum from STFT: |Y|^2
%   2. Estimate noise PSD using minimum statistics [3] or MMSE [7]
%   3. Calculate posterior SNR: gamma = |Y|^2 / noise_PSD
%   4. Estimate prior SNR using decision-directed approach [1]:
%      xi(t,f) = max(alpha*xb*xu(t-1,f) + (1-alpha)*max(gamma(t,f)-1, gn-1), xn)
%      where xu = gamma * gain^2 (unsmoothed prior SNR from previous frame)
%   5. Compute gain using specified estimator (via STFT_GetGain)
%
% The prior SNR estimation smooths between:
%   - Previous frame's estimate (weighted by alpha)
%   - Current frame's instantaneous estimate (weighted by 1-alpha)
% and applies floor constraints (xn, gn) and bias compensation (xb).
%
% Parameters for noise estimation algorithms:
%
% Minimum statistics noise estimate [3]: pp.ne=0 
%        pp.taca      % (11): smoothing time constant for alpha_c [0.0449 seconds]
%        pp.tamax     % (3): max smoothing time constant [0.392 seconds]
%        pp.taminh    % (3): min smoothing time constant (upper limit) [0.0133 seconds]
%        pp.tpfall    % (12): time constant for P to fall [0.064 seconds]
%        pp.tbmax     % (20): max smoothing time constant [0.0717 seconds]
%        pp.qeqmin    % (23): minimum value of Qeq [2]
%        pp.qeqmax    % max value of Qeq per frame [14]
%        pp.av        % (23)+13 lines: fudge factor for bc calculation  [2.12]
%        pp.td        % time to take minimum over [1.536 seconds]
%        pp.nu        % number of subwindows to use [3]
%        pp.qith      % Q-inverse thresholds to select maximum noise slope [0.03 0.05 0.06 Inf]
%        pp.nsmdb     % corresponding noise slope thresholds in dB/second [47 31.4 15.7 4.1]
%
% MMSE noise estimate [7]: pp.ne=1 
%        pp.tax      % smoothing time constant for noise power estimate [0.0717 seconds](8)
%        pp.tap      % smoothing time constant for smoothed speech prob [0.152 seconds](23)
%        pp.psthr    % threshold for smoothed speech probability [0.99] (24)
%        pp.pnsaf    % noise probability safety value [0.01] (24)
%        pp.pspri    % prior speech probability [0.5] (18)
%        pp.asnr     % active SNR in dB [15] (18)
%        pp.psini    % initial speech probability [0.5] (23)
%        pp.tavini   % assumed speech absent time at start [0.064 seconds]
%
% See also: STFT_GetGain, estnoisem, estnoiseg
%
% Original source: VOICEBOX ssubmmse function by Mike Brookes
% Refs:
%    [1] Ephraim, Y. & Malah, D.
%        Speech enhancement using a minimum-mean square error short-time spectral amplitude estimator
%        IEEE Trans Acoustics Speech and Signal Processing, 32(6):1109-1121, Dec 1984
%    [2] Ephraim, Y. & Malah, D.
%        Speech enhancement using a minimum mean-square error log-spectral amplitude estimator
%        IEEE Trans Acoustics Speech and Signal Processing, 33(2):443-445, Apr 1985
%    [3] Rainer Martin.
%        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
%        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.
%    [4] O. Cappe.
%        Elimination of the musical noise phenomenon with the ephraim and malah noise suppressor.
%        IEEE Trans Speech Audio Processing, 2 (2): 345???349, Apr. 1994. doi: 10.1109/89.279283.
%    [5] J. Erkelens, J. Jensen, and R. Heusdens.
%        A data-driven approach to optimizing spectral speech enhancement methods for various error criteria.
%        Speech Communication, 49: 530???541, 2007. doi: 10.1016/j.specom.2006.06.012.
%    [6] R. Martin.
%        Statistical methods for the enhancement of noisy speech.
%        In J. Benesty, S. Makino, and J. Chen, editors,
%        Speech Enhancement, chapter 3, pages 43???64. Springer-Verlag, 2005.
%    [7] Gerkmann, T. & Hendriks, R. C.
%        Unbiased MMSE-Based Noise Power Estimation With Low Complexity and Low Tracking Delay
%        IEEE Trans Audio, Speech, Language Processing, 2012, 20, 1383-1393

% Bugs/suggestions:
%   (1) sort out behaviour when si() is a matrix rather than a vector
%
%      Copyright (C) Mike Brookes 2004-2011
%      Version: $Id: ssubmmse.m 2460 2012-10-29 22:20:45Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yPSD=Y.*conj(Y);    % power spectrum of input speech
    
    nPSD = estimate_noise_PSD(yPSD,tinc, pp.ne);
    [xi_mat, gam] = estimate_SNR(yPSD, nPSD, snr_smoothing_alpha, pp);
end

function nPSD = estimate_noise_PSD(yPSD,tinc,ne)
    if ne>0
        [nPSD,~]=estnoiseg(yPSD,tinc);	% estimate the noise using MMSE
    else
        [nPSD,~]=estnoisem(yPSD,tinc);	% estimate the noise using minimum statistics
    end
end

function [prio, post] = estimate_SNR(yPSD, nPSD, snr_smoothing_alpha, pp)
    xu=1;                           % dummy unsmoothed SNR from previous frame
    xi_mat=zeros(size(yPSD));
    % yp PSD of speech, dp PSD hat of noise, gx max a priori SNR, gz min a
    % priori SNR
    gam=max(min(yPSD./nPSD,pp.gx),pp.gz);     % gamma = posterior SNR
    for i=1 : size(xi_mat, 1)
        gami=gam(i,:);
        % a...smoothing factor
        % xb..bias compensation?
        % xu..unsmoothed prior SNR
        % gn1.floor for a priori SNR
        % xn..floor for xi
        a=snr_smoothing_alpha; % SNR smoothing coefficient
        gn1=max(pp.gn-1,0); % floor for posterior SNR when estimating prior SNR
        xi=max(a*pp.xb*xu+(1-a)*max(gami-1,gn1),pp.xn);  % prior SNR
        xi_mat(i,:)=xi;
        gi=STFT_GetGain(xi, gami, pp.Estimator);
        xu=gami.*gi.^2;         % unsmoothed prior SNR
    end
    
    prio = xi_mat;
    post = gam;
end