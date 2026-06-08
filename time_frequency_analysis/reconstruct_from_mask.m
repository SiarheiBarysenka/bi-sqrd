function [xhat, Xhat] = reconstruct_from_mask(mask, S0, S, Nfft, anaWin, tinc, fs, y)
    % Reconstructs time-domain signal from magnitude mask by preserving phase
    %
    %   Input variables:
    %
    %          mask                    :  binary or soft magnitude mask (without f=0 column)
    %          S0                      :  STFT column at f = 0
    %          S                       :  STFT without f = 0 column
    %          Nfft                    :  FFT size
    %          anaWin                  :  analysis window
    %          tinc                    :  time increment between frames (in samples)
    %          fs                      :  sampling frequency (in hertz)
    %          y                       :  original noisy signal for gain normalization
    %
    %   Output variables:
    %
    %          xhat                    :  enhanced signal in time domain
    %          Xhat                    :  STFT of enhanced signal
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

    % Apply mask to magnitude while preserving f=0 column
    mag_xhat = [abs(S0) mask .* abs(S)];

    % Reconstruct full STFT from S0 and S
    Y = [S0 S];
    
    % Keep phase from noisy signal
    Xhat = mag_xhat .* (Y ./ abs(Y));
    
    % Inverse STFT
    xhat = overlapadd(irfft(Xhat, Nfft, 2), anaWin, tinc * fs);
end

