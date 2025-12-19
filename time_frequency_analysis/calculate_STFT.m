function [X, time, frequency, Nfft, tinc, w, x] = calculate_STFT(speech, fs, win_len, Nfft, overlap_factor, win)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE_STFT computes the Short-Time Fourier Transform of a speech signal
%
% This function performs time-frequency analysis by splitting the input speech
% signal into overlapping frames, applying a window function, and computing
% the real-valued FFT. It uses efficient real FFT computation and returns
% both the STFT and associated time/frequency axes.
%
% Usage:
%   [X, time, frequency, Nfft, tinc, w, x] = ...
%       calculate_STFT(speech, fs, win_len, Nfft, overlap_factor, win);
%
% Inputs:
%   speech          Input speech signal (samples x 1) or (1 x samples)
%   fs              Sampling frequency in Hz
%   win_len         Window length in samples
%   Nfft            FFT size (can be >= or < win_len for zero-padding/truncation)
%                   If Nfft > win_len: zero-padding is applied
%                   If Nfft < win_len: window is truncated
%   overlap_factor  Overlap between consecutive frames (0 < overlap < 1)
%                   Common values: 0.5 (50%), 0.75 (75%), 0.875 (87.5%)
%   win             Window function handle or empty for default (see split_into_frames)
%                   Examples: @hamming, @hann, @blackman, [] (sqrt-Hamming default)
%
% Outputs:
%   X               Complex STFT matrix (frames x frequency_bins)
%                   Contains only positive frequencies (real FFT)
%   time            Time vector for frame centers in seconds (frames x 1)
%   frequency       Frequency vector in Hz (1 x frequency_bins)
%                   Ranges from 0 to fs/2 (Nyquist frequency)
%   Nfft            FFT size (pass-through from input)
%   tinc            Frame increment (hop size) in seconds
%   w               Window function applied to frames (Nfft x 1)
%                   Adjusted to match Nfft size (zero-padded or truncated)
%   x               Framed signal before FFT (frames x win_len)
%                   Windowed and overlapped frames
%
% Notes:
%   - Uses real FFT (rfft) which is more efficient for real-valued signals
%   - Output X contains only positive frequencies (0 to fs/2)
%   - Frame increment: tinc = win_len * (1 - overlap_factor) / fs
%   - Number of frequency bins: Nfft / 2 + 1 (for real FFT)
%
% Example:
%   % Compute STFT with 32ms Hamming window, 75% overlap, 512-point FFT
%   fs = 16000;
%   win_len = 512;  % 32ms at 16kHz
%   [X, t, f] = calculate_STFT(speech, fs, 512, 512, 0.75, @hamming);
%
% See also: split_into_frames, rfft, reconstruct_from_mask
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
    winlen_sec = win_len / fs;
    
    [x, tt_x, ~, tinc, w] = split_into_frames(speech, fs, win, winlen_sec, overlap_factor, 0);
    
    X = rfft(x, Nfft, 2); % Speech spectrogram
    time = tt_x / fs; % frame times
    frequency = (1:size(X, 2)) ./ size(X, 2) .* (fs / 2);
    
    if Nfft > win_len
        % Add trailing zeros before returning the window
        w_temp = zeros(Nfft, 1);
        w_temp(1:win_len, 1) = w;
        w = w_temp;
    elseif Nfft < win_len
        % Truncate before returning the window
        w(Nfft + 1 : end, 1) = [];
    end
end