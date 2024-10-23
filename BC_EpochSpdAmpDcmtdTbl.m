function  outTbl=BC_EpochSpdAmpDcmtdTbl(theta_decimated, SG_decimated,FG_decimated,iv_inhb,ii)
%% The function BC_EpochSpdAmpTbl:
%            [Provide a general description of the function here]
%   Inputs:
%            -input1[type]: description 
% 
% 
%   Outputs:
%             -output1[type]: description 
% 
% 
%  First version BC 05-Aug-2024 
%%
if ~isstruct(theta_decimated) && ~isstruct(SG_decimated) && ~isstruct(FG_decimated)

    TtaAmp=NaN;
    TtaPwr=NaN;
    TtaPhi=NaN;

    SGAmp=NaN;
    SGPwr=NaN;
    SGPhi=NaN;

    FGAmp=NaN;
    FGPwr=NaN;
    FGPhi=NaN;
    Spd=NaN;
    epoch=ii;

else

    %The following code create the vector for the creation of the table
    thisThetAmp=restrict(theta_decimated, iv_inhb.tstart(ii), iv_inhb.tend(ii));
    thisSgAmp=restrict(SG_decimated, iv_inhb.tstart(ii), iv_inhb.tend(ii));
    thisFgAmp=restrict(FG_decimated, iv_inhb.tstart(ii), iv_inhb.tend(ii));
    %Create a table with  epoch#, condition, amplitude, velocity
    TtaAmp=thisThetAmp.data(1,:);
    TtaPwr=thisThetAmp.data(2,:);
    TtaPhi=thisThetAmp.data(3,:);

    SGAmp=thisSgAmp.data(1,:);
    SGPwr=thisSgAmp.data(2,:);
    SGPhi=thisSgAmp.data(3,:);

    FGAmp=thisFgAmp.data(1,:);
    FGPwr=thisFgAmp.data(2,:);
    FGPhi=thisFgAmp.data(3,:);
    Spd=thisThetAmp.data(4,:);
    epoch=repmat(ii,1,length(thisSgAmp.data));
    %Return such table as the fucntion
end
outTbl=table(epoch',TtaAmp',TtaPwr',TtaPhi',SGAmp',SGPwr',SGPhi',FGAmp',FGPwr',FGPhi',Spd', 'VariableNames', {'Epoch', 'ttaAmp','ttaPwr', 'ttaPhi','sgAmp','sgPwr','sgPhi','fgAmp','fgPwr','fgPhi','Spd'});
