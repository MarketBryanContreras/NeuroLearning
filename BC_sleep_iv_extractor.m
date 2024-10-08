function [iv_awake, iv_sws, iv_rem]=BC_sleep_iv_extractor(hypno)
%% The function BC_sleep_iv_extractor:
%            [Provide a general description of the function here]
%   Inputs:
%            -hypno[struct]: Structure containing the data and labels for sleep states obtained from "dSub_Sleep_screener" from the CEH2 repository 
% 
% 
%   Outputs:
%             -iv_awake[struct]: Intervals of time where the mouse was awake
%             -iv_sws[struct]: Intervals of time where the mouse was in slow wave sleep
%             -iv_rem[struct]: Intervals of time where the mouse was in rem
% 
% 
%  First version BC 21-May-2024 
%% Initialize and test for errors

%% Identify when data changes from 1,2,or 3
%hypno.data=[1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]%Test data
%jps=find(diff([0 hypno.data])); % All the jumps in the data are stored here
%% Format data to work with this function, it has to be a vertical vector
%Add a zero to the end to mark a change in state no matter what
sz=size(hypno.data,1);
if sz>1
    hypno.data=[(hypno.data)' 0];
else
    hypno.data=[(hypno.data) 0];
end
%% Collecting the time intervals for 1 values
jmps_1=diff(hypno.data==1); % Binarize data where 1(wake) was presnet and calculate difference to mark beginnings(1) and ends of states (-1)
% Adjust the 1st value to compensate for the diff function
if hypno.data(1)==1
    jmps_1=[1 jmps_1]; %add a 1 (begiining of state) in case is is the start state
else
    jmps_1=[0 jmps_1];% add a zero in case it is not
end
wake_end_idx=(find(jmps_1==-1)-1);% find all end of states idx, substract 1 to cover only states where the sleep state is happening
wake_end=hypno.tvec(wake_end_idx); %Get the time-stamps for these values
wake_start_idx=find(jmps_1==1); %find all the start of state idx
wake_start=hypno.tvec(wake_start_idx); %Get the time stamps for these idx
iv_awake=iv(wake_start,wake_end); %Convert to interval format
%Same for SWS
jmps_2=diff(hypno.data==2); 
if hypno.data(1)==2
    jmps_2=[1 jmps_2];
else
    jmps_2=[0 jmps_2];
end
sws_end_idx=(find(jmps_2==-1)-1);
sws_end=hypno.tvec(sws_end_idx);
sws_start_idx=find(jmps_2==1);
sws_start=hypno.tvec(sws_start_idx);
iv_sws=iv(sws_start,sws_end);
%Same for REM
jmps_3=diff(hypno.data==3); % All the jumps from 1 to something else
if hypno.data(1)==3
    jmps_3=[1 jmps_3];
else
    jmps_3=[0 jmps_3];
end
rem_end_idx=(find(jmps_3==-1)-1);
rem_end=hypno.tvec(rem_end_idx);
rem_start_idx=find(jmps_3==1);
rem_start=hypno.tvec(rem_start_idx);
iv_rem=iv(rem_start,rem_end);

%returning hypno to the original dimentions, 
sz2=length(hypno.data);
hypno.data= hypno.data(:,1:sz2-1); 
hypno.data=[(hypno.data)'];



