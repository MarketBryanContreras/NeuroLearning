%% Sandbox power spectrum
function amp= BC_power_phase(sig)
%% MS_ModIdx_win: Simple fx to obtain the power of a already filtered signal
%                 and calculate the mean of it. Indended for theta - gamma modulation as per Tort et al. 2010/2009/2008
%
%    Inputs: 
%    - sig: [struct] filtered LFP data in the csc format (vdmlab codebase) to be used for the phase component.  
%
% 
%    Outputs: 
%    - amp: 
%
%
%
% BC 2023-08-24   initial version 
%
amp=sig;
amp.data=smooth(abs(hilbert(sig.data)), floor(sig.cfg.hdr{1}.SamplingFrequency*0.1))';
amp.data(2,:)=(amp.data(1,:)).^2;
amp.data(3,:)= angle(hilbert(sig.data));
amp.label={'amplitude', 'power', 'phase'};
amp.units={'V', 'W'};

