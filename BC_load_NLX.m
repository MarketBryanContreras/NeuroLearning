function [csc, evts, pos] = BC_load_NLX(cfg_in)
%% BC_load_NLX:
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-11-03   initial version 
%
%
%
%% initialize


evts = LoadEvents([]);

pos = MS_DLC2TSD(cd, [], [7.2 7.2]);


csc = MS_LoadCSC(cfg_in);


end