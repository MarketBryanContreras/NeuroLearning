function SpdPwrTbl=BC_SpdPwrTblGenerator(out_archt,inSufix)
%% The function BC_SpdPwrTblGenerator:
%            [This function generates a table with the speed and their respective theta, sg, and fg amplitude, pwr and phase for following analysis. Used in SUPER_BC_stats]
%   Inputs:
%            -out_archt[struct]: out strcuture from BC_pipeline
%            -inSufix[string]: string that selects the suffic of the tables to gather from all the subjects, either inhb, noInhb or running
% 
%   Outputs:
%             -SpdPwrTbl[table]: description 
% 
% 
%  First version BC 30-Jul-2024 
%% 
if nargin<2
    inSufix='inhb';
end
sufix={'inhb', 'noInhb','running'};
if ~ismember(inSufix,sufix);
    error('Please select an suffix from inhb, noInhb or running')
else
    archtSubject=[];
    %Inhibition
    archtInbTtaPhi=[];
    archtInbTtaAmp=[];
    archtInbTtaPwr=[];
    archtInbSGPhi=[];
    archtInbSGAmp=[];
    archtInbSGPwr=[];
    archtInbFGPhi=[];
    archtInbFGAmp=[];
    archtInbFGPwr=[];
    archtInbFGPhi=[];
    archtInbSpeed=[];
    archtInbEpoch=[];
    archt_list=fieldnames(out_archt);

    for iSub = 1:length(archt_list)


        archtInbSpeed=[archtInbSpeed ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).spd];
        archtSubject= [archtSubject repmat(iSub,1,length(out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).tta_phi))];

        archtInbTtaAmp=[archtInbTtaAmp ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).tta_amp];
        archtInbTtaPwr=[archtInbTtaPwr ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).tta_pwr];
        archtInbTtaPhi=[archtInbTtaPhi ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).tta_phi];

        archtInbSGAmp=[archtInbSGAmp ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).sg_amp];
        archtInbSGPwr=[archtInbSGPwr ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).sg_pwr];
        archtInbSGPhi=[archtInbSGPhi ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).sg_phi];

        archtInbFGAmp=[archtInbFGAmp ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).fg_amp];
        archtInbFGPwr=[archtInbFGPwr ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).fg_pwr];
        archtInbFGPhi=[archtInbFGPhi ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).fg_phi];

        archtInbEpoch=[archtInbEpoch ;out_archt.(archt_list{iSub}).D4.(['speedPwrTbl_' inSufix]).epoch];

        if iSub==length(archt_list)
            SpdPwrTbl= table(archtSubject', archtInbEpoch,archtInbSpeed,archtInbTtaAmp,archtInbTtaPwr,archtInbTtaPhi, archtInbSGAmp,archtInbSGPwr,archtInbSGPhi,archtInbFGAmp,archtInbFGPwr,archtInbFGPhi,'VariableNames', {'Subject', 'Epoch','Spd','TtaAmp','TtaPwr','TtaPhi','SGAmp','SGPwr','SGPhi','FGAmp','FGPwr','FGPhi'});
        end
    end
end