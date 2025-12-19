function [kappa_result] = calculate_kappa(input, Nwinfilt)
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

    filter_b = ones(Nwinfilt, 1);
    filter_a = 1;
    
    sin_filt = filter(filter_b, filter_a, sin(input));
    cos_filt = filter(filter_b, filter_a, cos(input));
    
    r = abs(cos_filt + 1i .* sin_filt) / Nwinfilt; % Matrix of the mean resultant lengths
    kappa = kappa_from_mean_resultant_lengths(r, Nwinfilt); % Kappa
    
    %% Compensate the delay from the `filter` operation by re-assigning the indices
    kappa_result = zeros(size(input));
    shift = Nwinfilt - ceil(Nwinfilt / 2);
    kappa_result(1 : end - shift, :) = kappa(shift + 1 : end, :);
end