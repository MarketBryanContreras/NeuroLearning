function iv_in_min=BC_iv2min(iv_in_sec)
%% The function BC_iv2min:
%            Converts the time stamps contained in an Iv in the vandermer
%            lab format to minutes
%   Inputs:
%            -iv_in_sec[type]: description 
% 
% 
%   Outputs:
%             -in_min[type]: description 
% 
% 
%  First version BC 06-Nov-2023 
%% Quick convertion by division

iv_in_min=iv_in_sec;
iv_in_min.tstart=iv_in_min.tstart/60;
iv_in_min.tend=iv_in_min.tend/60;
