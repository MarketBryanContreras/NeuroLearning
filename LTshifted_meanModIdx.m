function [shift_mean,shift_std]=LTshifted_meanModIdx(sig_phi,sig_amp)
%% Lets compute the mean amplitude of the fast oscillation in the bins of the slow oscilation phase
phi= angle(hilbert(sig_phi.data(1,:)));
amp= smooth(abs(hilbert(sig_amp.data)), floor(sig_amp.cfg.hdr{1}.SamplingFrequency*0.1))';
nshuffles=500;

%2.reate a composite of the phase_theta, amplitud_SG
phi_bins = -pi:pi/9:pi; 
 
shiftedModIdx=[];
l=length(phi)-1;
for xx=nshuffles:-1:1
    shiftedPhi=circshift(phi,randi(l));
    [~,bins_idx]= histc(shiftedPhi, phi_bins);
    amp_means= zeros(1,length(unique(bins_idx)));
    for ii= 1: length(unique(bins_idx))
        amp_means(ii)= nanmean(amp(bins_idx == ii));
    end
    norm_phase_amp_bins= amp_means/sum(amp_means);
    shiftedModIdx(xx)=MS_ModIdx(norm_phase_amp_bins);
end
shift_mean=mean(shiftedModIdx);
shift_std=std(shiftedModIdx);

