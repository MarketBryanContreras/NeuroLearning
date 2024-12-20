function c_vector=BC_color_genertor(c_name, alpha)
%% The function BC_color_genertor:
%            Provides an array corrsponding to the color specified in color
%            name ccording to MAtlab format. If alpha(transparency) is
%            provided it will return an array with the corresponding alpha
%            value as well. 
%            
%   Inputs:
%            -c_name[string]: Name of the color to call 
%            -alpha[double]: Value from 0(totally translucid) to 1 (completly opaque)
%                            that assigns how transparent the color will be.
% 
% 
%   Outputs:
%             -output1[type]: Array in MAtlab format with the desired color for the OLM project 
% 
% 
%  First version BC 05-Nov-2023 
%% Assigning color pallete
colors= {'Archt_green', 'Oxford_blue', 'Powder_blue', 'Burnt_orange', 'Red_crayola', 'Web_orange','Swamp_green','Torment_blue','indigo','turkey_red'};
color_arrays = {
    [0.530 0.820 0.645],
    [0.039 0.137 0.259],
    [0.682 0.772 0.921],
    [0.758 0.348 0.249],
    [0.937 0.176 0.337],
    [1.00 0.678 0.020],
    [0.345 0.529 0.419],
    [0.474 0.537 0.639],
    [0.392 0.016 0.527],
    [0.694 0.059 0]};
% Convert color names to lower case for case-insensitive comparison
c_name_lower = lower(c_name);

% Find the index of the matching color (case insensitive)
matching_index = find(strcmpi(lower(colors), c_name_lower));

if isempty(matching_index)
    % If no exact match is found, return list of available colors
    color_list = strjoin(colors, ', ');
    error(['Color not found. Available colors are: ', color_list]);
else
    % Return the numeric array for the matching color
    if nargin==1
        c_vector = [color_arrays{matching_index}];
    else nargin==2
        if alpha <0 || alpha>1
            error('Alpha can only be a value between 0 and 1');
        else
            c_vector = [color_arrays{matching_index}, alpha];
        end
    end
end

    
