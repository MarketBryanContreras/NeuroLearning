function DcmtAmp=BC_ampSpdDcmt(filteredSig, pos)
%% The function BC_ampSpdDcmt:
%            [Provide a general description of the function here]
%   Inputs:
%            -input1[type]: description 
% 
% 
%   Outputs:
%             -output1[type]: description 
% 
% 
%  First version BC 01-Aug-2024 
%% 
%% decimate to match sampling frequency
cfg_deci = [];
cfg_deci.decimateFactor = 20; 

DcmtAmp = decimate_tsd(cfg_deci,filteredSig); 
% interpolate the spd to match the downsampled LFP amp
DcmtAmp.data(2,:) = interp1(pos.tvec, pos.data(5,:), theta_amp_d.tvec, 'next', 'extrap');
DcmtAmp.label= {'Amp', 'Spd'};