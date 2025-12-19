function kappa = kappa_from_mean_resultant_lengths(R, N)
    %% This is the same math as in `circ_kappa`, but rewritten in the matrix form.
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

    kappa_less_053 = 2*R + R.^3 + 5*R.^5/6;
    kappa_less_085 = -.4 + 1.39*R + 0.43./(1-R);
    kappa_other = 1./(R.^3 - 4*R.^2 + 3*R);

    I_less_053 = false(size(R));
    I_less_053(R < 0.53) = true;
    
    I_less_085 = false(size(R));
    I_less_085(and(R>=0.53, R<0.85)) = true;
    
    I_other = ~or(I_less_053, I_less_085);
    
    kappa = I_less_053 .* kappa_less_053 ...
          + I_less_085 .* kappa_less_085 ...
          + I_other .* kappa_other;

    if N<15 && N>1
        kappa_less_2 = max(kappa-2*(N*kappa).^-1,0);
        kappa_more_2 = (N-1)^3*kappa./(N^3+N);
        
        I_less_2 = false(size(kappa));
        I_less_2(kappa < 2) = true;
        
        kappa = I_less_2 .* kappa_less_2 + ~I_less_2 .* kappa_more_2;
    end
end