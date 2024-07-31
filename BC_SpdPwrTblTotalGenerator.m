function SpdPwrTbl=BC_SpdPwrTblTotalGenerator(struct,inSufix)
%% The function BC_SpdPwrTblGenerator:
%            [This function generates a table with the speed and their respective theta, sg, and fg amplitude, pwr and phase for following analysis]
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
    TtaPhi=[];
    TtaAmp=[];
    TtaPwr=[];
    SGPhi=[];
    SGAmp=[];
    SGPwr=[];
    FGPhi=[];
    FGAmp=[];
    FGPwr=[];
    Epoch=[];
    subject_list=fieldnames(struct);

    for iSub = 1:length(subject_list)


        
        archtSubject= [archtSubject repmat(iSub,1,length(struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).tta_phi))];

        TtaAmp=[TtaAmp ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).tta_amp];
        TtaPwr=[TtaPwr ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).tta_pwr];
        TtaPhi=[TtaPhi ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total'] ).tta_phi];

        SGAmp=[SGAmp ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).sg_amp];
        SGPwr=[SGPwr ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).sg_pwr];
        SGPhi=[SGPhi ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).sg_phi];

        FGAmp=[FGAmp ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).fg_amp];
        FGPwr=[FGPwr ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).fg_pwr];
        FGPhi=[FGPhi ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).fg_phi];

        Epoch=[Epoch ;struct.(subject_list{iSub}).D4.(['speedPwrTbl_' inSufix 'Total']).epoch];

        if iSub==length(subject_list)
            SpdPwrTbl= table(archtSubject', Epoch,TtaAmp,TtaPwr,TtaPhi, SGAmp,SGPwr,SGPhi,FGAmp,FGPwr,FGPhi,'VariableNames', {'Subject', 'Epoch','TtaAmp','TtaPwr','TtaPhi','SGAmp','SGPwr','SGPhi','FGAmp','FGPwr','FGPhi'});
        end
    end
end