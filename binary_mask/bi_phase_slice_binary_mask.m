function [mask, kappa_bi_phase, cosBiPhDev, null_hypothesis_threshold] = bi_phase_slice_binary_mask(bi_phase, k1_vector, k2_vector, k3_vector, Nwinfilt, prio, post, p_H0)
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

    kappa_bi_phase = calculate_kappa(bi_phase, Nwinfilt);
    cosBiPhDev = cosine_bi_phase_dev(k1_vector, k2_vector, k3_vector, prio', post');
    null_hypothesis_threshold = bi_phase_null_hypothesis_threshold(kappa_bi_phase', p_H0);
    mask = decision_to_reject_hull_hypothesis(null_hypothesis_threshold, cosBiPhDev);
end