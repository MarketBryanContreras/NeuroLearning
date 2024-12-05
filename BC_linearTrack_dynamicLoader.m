function [data_dir,inter_dir, mouse_group]=BC_linearTrack_dynamicLoader(group)
%% The function Dynamic loader: Goes to the dyrectroy of experiemntal orcontrol groups to continue with the the pipeline
%            
%   Inputs:
%            -group[string]: experimental for archt or control for yfp 
% 
% 
%   Outputs:
%             none
% 
% 
%  First version BC 29-Jul-2024 
validInputs = {'control', 'experimental'};
if nargin == 0
    group= 'experimental';
elseif ischar(group) || isstring(group)
    group=lower(group);
    if ~ismember(group, validInputs)
        error('Invalid input. The only acceptable options are "control" or "experimental".');
    else
        sys=computer;
        if contains(sys,'PCWIN')
            % Dynamic based on computer user windows.
            if group=="experimental"
                data_dir = [getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\Archt'];
            else
                data_dir = [getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\eyfp'];
            end
            inter_dir=[getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter']
        elseif contains(sys,'MAC')
            % Dynamic based on computer user for mac
            if group=="experimental"
                data_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'raw' filesep 'Archt'];
            else
                data_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'raw' filesep 'eyfp'];
            end
            inter_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'inter'];
        else disp("This system is not supported on the dynamic loader yet, pleae add it or contact BC for it")
        end
        cd(data_dir)
        if group=="experimental"
            mouse_group=1;
        else
            mouse_group=2;
        end
    end
    
end
