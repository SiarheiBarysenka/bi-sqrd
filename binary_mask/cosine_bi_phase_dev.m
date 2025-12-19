function result = cosine_bi_phase_dev(k1_vector, k2_vector, k3_vector, prio, post)
    % This function provides a result equivalent to the following:
    %
    %         cos1 = cosine_phase_dev(prio(k1, :), post(k1, :));
    %         cos2 = cosine_phase_dev(prio(k2, :), post(k2, :));
    %         cos3 = cosine_phase_dev(prio(k3, :), post(k3, :));
    %
    %         sin1 = sqrt(1 - cos1 .^ 2);
    %         sin2 = sqrt(1 - cos2 .^ 2);
    %         sin3 = sqrt(1 - cos3 .^ 2);
    %
    %         result = cos1 .* cos2 .* cos3...
    %             + cos1 .* sin2 .* sin3...
    %             + sin1 .* cos2 .* sin3...
    %             - sin1 .* sin2 .* cos3;
    %
    % where
    %       function result = cosine_phase_dev(prio, post)
    %           result = (prio + post - 1) ./ (2 * sqrt(prio .* post));
    %       end

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

    prio_k1 = prio(k1_vector, :);
    prio_k2 = prio(k2_vector, :);
    prio_k3 = prio(k3_vector, :);
    
    post_k1 = post(k1_vector, :);
    post_k2 = post(k2_vector, :);
    post_k3 = post(k3_vector, :);
    
    result = 1 ./ (8 * alpha_snr(prio_k1, prio_k2, prio_k3, post_k1, post_k2, post_k3)) .* ( ...
          beta_snr(prio_k1, prio_k2, prio_k3, post_k1, post_k2, post_k3) ...
        + beta_snr(prio_k2, prio_k1, prio_k3, post_k2, post_k1, post_k3) ...
        - beta_snr(prio_k3, prio_k1, prio_k2, post_k3, post_k1, post_k2) ...
        + gamma_snr(prio_k1, prio_k2, prio_k3, post_k1, post_k2, post_k3) ...
    );
end

function result = xi_snr(prio, post)
    expression_under_square_root = ...
          2 * post .* (1 + prio) ...
        - (prio - 1) .^ 2 ...
        - post .^ 2 ...
    ;
    result = sqrt(expression_under_square_root);
end

function result = rho_snr(prio, post)
    result = prio + post - 1;
end

function result = gamma_snr(prio_k1, prio_k2, prio_k3, post_k1, post_k2, post_k3)
    result = rho_snr(prio_k1, post_k1) ...
          .* rho_snr(prio_k2, post_k2) ...
          .* rho_snr(prio_k3, post_k3);
end

function result = beta_snr(prio_k1, prio_k2, prio_k3, post_k1, post_k2, post_k3)
    result = rho_snr(prio_k1, post_k1) ...
          .* xi_snr(prio_k2, post_k2) ...
          .* xi_snr(prio_k3, post_k3);
end

function result = alpha_snr(prio_k1, prio_k2, prio_k3, post_k1, post_k2, post_k3)
    result = sqrt(prio_k1 .* post_k1) ...
          .* sqrt(prio_k2 .* post_k2) ...
          .* sqrt(prio_k3 .* post_k3);
end