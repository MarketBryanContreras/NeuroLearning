function interaction_int=BC_object_interaction_intervals(idx_in, min_frames)
%% The function BC_object_interaction_counting:
%            [Provide a general description of the function here]
%   Inputs:
%            -idx_in[1xn Double]: Vector containinf the idx of the frames where the mouse whas in the region of interest 
%            -min_frames[Double]: Minimun number of frames that the mouse has to be present in the region for be considered an interval  
% 
%   Outputs:
%             -interaction_int[Double]: Matrix of the frames (beginning,end [nIntx2]) where
%             the mouse interacted with the objet for more than the input
%             minimum number of frames 
% 
% 
%  First version BC 25-Jan-2024 
%% 
differences=diff([0,idx_in]);
%The smoothing below may cause that some point are out of the specified radious
differences(differences<=29)=1; %smoothing for frames when the mouse returned no longer than 29 frames( 1 seg) 
count=0;
final_intervals=[];
for ii=length(differences):-1:1
    if differences(ii)==1
        if count==0
            temp_start=idx_in(ii);
            count=count+1;
        else
            count=count+1;
        end
    else
        if count>=min_frames
            final_intervals=[final_intervals temp_start idx_in(ii)];
            count=0;
        else
            count=0;
        end
    end
end
final_intervals= fliplr(final_intervals);
interaction_int=(reshape(final_intervals, [2,length(final_intervals)/2]))';


