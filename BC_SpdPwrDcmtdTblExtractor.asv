function bigTbl=BC_SpdPwrDcmtdTblExtractor(out,inPrefix)
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
    inPrefix='Inhb';
end
prefix={'Inhb', 'NoInhb','Running'};
if ~ismember(inPrefix,prefix);
    error('Please select an prefix from Inhb, NoInhb, or Running')
else
    SubjectList=[];
    bigTbl=[];
    
    %Inhibition
   
    subject_list=fieldnames(out);
    

    for iSub = 1:length(subject_list)
        bigTbl= [bigTbl; out.(subject_list{iSub}).D4.([inPrefix 'DcmtSpeedAmpTbl'])];
        SubjectList= [SubjectList; repmat(iSub,height(out.(subject_list{iSub}).D4.([inPrefix 'DcmtSpeedAmpTbl'])),1)];

        if iSub==length(subject_list)
            bigTbl=addvars(bigTbl,SubjectList,'NewVariableNames','Subject', 'Before', 1);
        end
    end
end