function info=infoSession(parts)
%% The function infoSession:
%            [Provide a general description of the function here]
%   Inputs:
%            -input1[type]: description 
% 
% 
%   Outputs:
%             -output1[type]: description 
% 
% 
%  First version BC 06-Dec-2024 
%% 
    parts= split(parts,filesep);
    parts=parts{end};
    parts= split(parts,'_');
    info.subject=parts{1};
    info.date=[ parts{2} '_' parts{3} '_' parts{4}];
    info.session=parts{5};
    