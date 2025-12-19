function G_l= STFT_GetGain(xi_l,zeta_l,Method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STFT_GetGain computes spectral gain for various speech enhancement estimators
%
% This function calculates the time-frequency gain function used to suppress
% noise in speech signals. It implements multiple classical and modern gain
% estimation methods based on prior and posterior SNR estimates.
%
% Adapted by Siarhei Y. Barysenka from the original source code by Johannes Stahl:
% https://gitlab.com/johannesstahl/PACO/-/blob/master/Functions/STFT_GetGain.m
%
% Usage:
%   G_l = STFT_GetGain(xi_l, zeta_l, Method);
%
% Inputs:
%   xi_l    Prior SNR (a priori SNR) at frame l
%           Can be scalar or matrix (frames x frequency bins)
%   zeta_l  Posterior SNR (a posteriori SNR) at frame l  
%           Must have same dimensions as xi_l
%   Method  String specifying the gain estimator to use:
%           'WienerFilter' - Classical Wiener filter gain
%           'MMSESTSA'     - Minimum mean-square error short-time spectral amplitude [2]
%           'LSA'          - Log spectral amplitude estimator [3]
%           'JMAPLV'       - Joint MAP (Lotter and Vary) [4]
%           'JMAPGW'       - Joint MAP (Wolfe and Godsill), Eq. 29 [1]
%           'MAPGW'        - MAP (Wolfe and Godsill), Eq. 36 [1]
%           'MMSEGW'       - MMSE (Wolfe and Godsill), Eq. 39 [1]
%
% Output:
%   G_l     Spectral gain for frame l (same dimensions as inputs)
%           The enhanced spectrum is obtained by: S_hat = G_l .* Y
%           NaN values are automatically set to 0
%
% References:
% [1] P.J. Wolfe, S.J. Godsill "Efficient alternatives to the Ephraim and
%     Malah suppression rule for audio signal enhancement" (2003), pp. 1043–1051
% [2] Y. Ephraim and D. Malah, "Speech enhancement using a minimum-
%     mean square error short-time spectral amplitude estimator," IEEE Trans.
%     Audio, Speech, and Language Process., vol. 32, no. 6, pp. 1109–1121,
%     1984.
% [3] Y. Ephraim and D. Malah, "Speech enhancement using a minimum mean-square error log-
%     spectral amplitude estimator," IEEE Trans. Audio, Speech, and Language
%     Process., vol. 33, no. 2, pp. 443–445, 1985.
% [4] T. Lotter, P. Vary "Speech enhancement by MAP spectral amplitude
%     estimation using a super-gaussian speech model" EURASIP J. Adv.Signal Process., 2005 (7) (2005), pp. 1110–1126
%
% See also: snr_decision_directed, estimate_SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STFT_GetGain.m computes Gain for various estimators.
%     Copyright (C) 2017 Johannes Stahl
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch Method
    case 'WienerFilter'
        G_l = xi_l./(1+xi_l);
	case 'MMSESTSA' % [1]
        % default compute gain as proposed in [2]
        vk = xi_l.*zeta_l./(1+xi_l); % eq. 13 [1]
        G_l_EM = (0.277 + vk)./zeta_l; % accurate to 0.02 dB for v>1
        if any(vk<1)
           G_l_EM(vk<1) = gamma(1.5).*sqrt(vk(vk<1)).*((1+vk(vk<1))...
                          .*besseli(0,vk((vk<1))/2) + vk(vk<1).*besseli(1,vk(vk<1)/2))./(zeta_l(vk<1).*exp(vk(vk<1)/2)); % eq. 7 [1] 
        end
        G_l = G_l_EM; % [2]
    case 'LSA'
        G_l = xi_l./(1+xi_l).*exp(0.5*expint(xi_l./(1+xi_l).*zeta_l));
	case 'JMAPLV' % [3]
        nu = 0.2;
        mu = sqrt((nu+1)*(nu+2));
        u = 1/2 - mu./(4*sqrt(zeta_l.*xi_l));
        G_l = u + sqrt(u.^2 + nu./(2*zeta_l));
	case 'JMAPGW' % eq. 29 [1]
    	G_l = (xi_l + sqrt(xi_l.^2 +2*(1+xi_l).*xi_l./zeta_l))./(2*(1+xi_l));
    case 'MAPGW' % eq. 36 [1]
        G_l = (xi_l + sqrt(xi_l.^2 +(1+xi_l).*xi_l./zeta_l))./(2*(1+xi_l));
    case 'MMSEGW' % eq. 39 [1]
        vk = xi_l.*zeta_l./(1+xi_l); % eq. 13 [1]
        G_l = sqrt(xi_l./(1+xi_l).*(1+vk)./zeta_l);
end

G_l(isnan(G_l)) = 0;

end