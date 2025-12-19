function result = bi_phase_mask_STFT_domain(mask_B2, mask_B3)
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

    % Set the size of the resulting mask:
    %   - time dimension (1-st dimension) corresponds to the size of the 1-st dimension of `mask_B2`
    %   - frequency dimension (2-nd dimension) corresponds to the doubled size of the 2-nd dimension of `mask_B2`
    mask_size = [size(mask_B2, 1), size(mask_B2, 2) * 2];
    
    % Initialize with all `false`
    result = false(mask_size);
    
    % Populate all even-indexed frequency bin vectors (2, 4, 6, ...) from the `mask_B2` (diagonal slice).
    result(:,2:2:end) = mask_B2;
    
    % Populate all odd-indexed frequency bin vectors (3, 5, 7, ...) from the `mask_B3` (next-to-diagonal slice).
    result(:,3:2:end) = mask_B3;
end