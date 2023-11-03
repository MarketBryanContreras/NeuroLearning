function norm_phase_amp_bins= BC_phase_amp_norm_bins(sig_phi,sig_amp)
%% MS_ModIdx_win: It bins the amp of the fast oscillation in the phases bin
%                 and calculate the mean of it. Indended for theta - gamma modulation as per Tort et al. 2010/2009/2008
%
%    Inputs: 
%    - sig_phi: [struct] filtered LFP data in the csc format (vdmlab codebase) to be used for the phase component.  
%
%    - sig_amp: [struct] filtered LFP data in the csc format (vdmlab codebase) to be used for the amplitude component. 
%
%    
%
%    Outputs: 
%    - norm_phase_amp_bins: [18 x 1] Mean_amplitude_in 18 phase bins 
%
%
%
% BC 2023-08-24   initial version 
%

%% Lets compute the mean amplitude of the fast oscillation in the bins of the slow oscilation phase
phi= angle(hilbert(sig_phi.data(1,:)));
amp= smooth(abs(hilbert(sig_amp.data)), floor(sig_amp.cfg.hdr{1}.SamplingFrequency*0.1))';
%2.reate a composite of the phase_theta, amplitud_SG
phi_bins = -pi:pi/9:pi; 
[~,bins_idx]= histc(phi, phi_bins); %This creates a vector with the corresponding bin of each phase
%3. The phases are binned according to which phase they belong to and their
%mean is calculated
amp_means= zeros(1,length(unique(bins_idx)));
for ii= 1: length(unique(bins_idx));
    amp_means(ii)= nanmean(amp(bins_idx == ii));
end
%4.Each mean amplitude is then divided by the sum over the bins 
norm_phase_amp_bins= amp_means/sum(amp_means);

