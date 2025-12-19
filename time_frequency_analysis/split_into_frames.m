function [y, tt, nf, tinc, w] = split_into_frames(s, fs, win, winlen_sec, overlap_factor, ri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPLIT_INTO_FRAMES segments a speech signal into overlapping windowed frames
%
% This function performs frame-based segmentation of a continuous speech signal
% with configurable overlap and window functions. It is designed for STFT
% analysis and ensures proper normalization for overlap-add reconstruction.
%
% Usage:
%   [y, tt, nf, tinc, w] = split_into_frames(s, fs, win, winlen_sec, ...
%                                             overlap_factor, ri);
%
% Inputs:
%   s                   Input speech signal (samples x 1) or (1 x samples)
%   fs                  Sampling frequency in Hz
%   win                 Window function handle or string
%                       - If empty or not provided: sqrt-Hamming window (default)
%                       - If function handle: called as w = win(nf)
%                       - Common options: @hamming, @hann, @blackman, etc.
%   winlen_sec          Window length in seconds
%                       Common values: 0.016, 0.032 (16ms, 32ms)
%   overlap_factor      Frame overlap as fraction (0 < overlap_factor < 1)
%                       - 0.5   = 50% overlap (hop = 50% of window)
%                       - 0.75  = 75% overlap (hop = 25% of window)
%                       - 0.875 = 87.5% overlap (hop = 12.5% of window)
%   ri                  Round frame increment flag [optional, default: 0]
%                       - 0: Use exact frame increment from overlap_factor
%                       - 1: Round to nearest power of 2 samples (for efficiency)
%
% Outputs:
%   y                   Framed signal (frames x window_length)
%                       Each row is a windowed frame
%   tt                  Frame center times in samples (frames x 1)
%   nf                  Window length in samples (FFT length)
%   tinc                Actual frame increment (hop size) in seconds
%   w                   Window function applied (1 x nf)
%                       Normalized for unity gain with overlap-add
%
% Algorithm:
%   1. Calculate window length: winlen = winlen_sec * fs
%   2. Calculate overlap in samples: overlap = round(winlen * overlap_factor)
%   3. Calculate hop size: hoptime = winlen - overlap
%   4. Optionally round hop to power of 2 (if ri=1)
%   5. Generate or normalize window function
%   6. Normalize window for overlap-add: w = w / sqrt(sum(w(1:ni:nf).^2))
%   7. Split signal into frames using enframe
%
% Window Normalization:
%   The window is normalized so that overlap-add reconstruction preserves
%   the signal amplitude: w = w / sqrt(sum(w(1:ni:nf).^2))
%   This ensures perfect reconstruction when using synthesis with same overlap.
%
% Frame Increment:
%   - If ri=0: ni = round(ti*fs) where ti = hoptime / fs
%   - If ri=1: ni = 2^(nextpow2(ti*fs*sqrt(0.5)))
%   The power-of-2 rounding can improve FFT efficiency at slight cost to timing.
%
% Example:
%   % Split signal into frames with 32ms window, 75% overlap
%   [frames, t, nf, tinc, w] = split_into_frames(speech, 16000, ...
%                                                 @hamming, 0.032, 0.75, 0);
%
% See also: enframe, calculate_STFT, reconstruct_from_mask
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

    winlen = winlen_sec * fs;
    overlap = round(winlen * overlap_factor);
    hoptime = winlen - overlap;
    ti = hoptime / fs;         % desired frame increment
    
    if exist('ri','var') && ri == 1
        ni=pow2(nextpow2(ti*fs*sqrt(0.5)));
    else
        ni=round(ti*fs);    % frame increment in samples
    end
    tinc=ni/fs;         % true frame increment time
    no = round(winlen / hoptime);     % overlap factor = (fft length)/(frame increment)
    nf = ni*no;                           % fft length

    if ~exist('win','var') || isempty(win)
        w = sqrt(hamming(nf+1))'; w(end)=[];  % for default use sqrt hamming window
    else
        w = win(nf);
    end
    w=w/sqrt(sum(w(1:ni:nf).^2));       % normalize to give overall gain of 1

    [y,tt]=enframe(s,w,ni,'r');
end