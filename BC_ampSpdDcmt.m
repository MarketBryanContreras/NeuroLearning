function dcmtAmp=BC_ampSpdDcmt(filteredSig, pos)
%% The function BC_ampSpdDcmt:
%            This fucntion decimates the LFP signal and interpolates the tracking position obtained from DLC to match the number of sampling points.
%            It povides the following information in the order of rows: 1)Amplitude, 2)Power, 3)Phase, 4)Speed 5)Position in the x axis
%   Inputs:
%            -filteredSig[struct]: Filtered signal recorded with the vandermer lab struct
%            -pos[struct]: Position strcuture from the DLCtoTSD function from Eric C.
% 
%   Outputs:
%            -dcmtAmp[type]: Structure with the amp, power, phi and speed
% 
% 
%  First version BC 01-Aug-2024 
 %% Obataining the amp, pwr and phase in a unidemnsional tsd so it can be used by decimate TSD
    amp=filteredSig;
    pwr=filteredSig;
    phi=filteredSig;
    amp.data = smooth(abs(hilbert(filteredSig.data)), floor(filteredSig.cfg.hdr{1}.SamplingFrequency*0.1))';
    pwr.data= amp.data.^2;
    phi.data = angle(hilbert(filteredSig.data));

    %% decimate to match sampling frequency
    cfg_deci = [];
    cfg_deci.decimateFactor = 20;

    %Decimating the amplitude, pwr and phi
    dcmtAmp = decimate_tsd(cfg_deci,amp);
    dcmtPwr= decimate_tsd(cfg_deci,pwr);
    dcmtPhi= decimate_tsd(cfg_deci,phi);
    
    %Interpolate the spd to match the downsampled LFP amp
    spd_interp = interp1(pos.tvec, pos.data(5,:), dcmtAmp.tvec, 'next', 'extrap');
    xpos_interp= interp1(pos.tvec, pos.data(3,:), dcmtAmp.tvec, 'next', 'extrap');
    
    %NORMALIZING THE AMP, PWR AND Phi VALUES BY THE 98TH PERCENTILE
    Amp98=dcmtAmp.data./(prctile(dcmtAmp.data,98)); %Dividing the values by the 98th percentile
    Amp98(Amp98>=1)=1; %Any values bigger than 1 are assigned to 1
    Pwr98=dcmtPwr.data./(prctile(dcmtPwr.data,98)); %Dividing the values by the 98th percentile
    Pwr98(Pwr98>=1)=1; %Any values bigger than 1 are assigned to 1
    Phi98=dcmtPhi.data./(prctile(dcmtPhi.data,98)); %Dividing the values by the 98th percentile
    Phi98(Phi98>=1)=1; %Any values bigger than 1 are assigned to 1

    %Putting pwr, amp and phi into a single strcuture
    dcmtAmp.data(1,:)=Amp98; 
    dcmtAmp.data(2,:)=Pwr98;
    dcmtAmp.data(3,:)=Phi98;
    dcmtAmp.data(4,:)=spd_interp;
    dcmtAmp.data(5,:)=xpos_interp;
    dcmtAmp.label= {'Amp','Pwr','Phi','Spd','Xpos'};
    dcmtAmp.units= {'V','W','Rad','cm/s', 'cm'};