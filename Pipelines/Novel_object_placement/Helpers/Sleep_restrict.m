function [csc_s,emg_s]=Sleep_restrict(evts,csc,emg,info)
%% The function Sleep_restrict:
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
%Extractinc the time stamps of the recording
    start_OF = evts.t{find(contains(evts.label, 'Starting Recording'))}(1); %Choose the time stamp of the first "Starting recording" which corresponds to the start of the OF
    start_sleep = evts.t{find(contains(evts.label, 'Starting Recording'))}(2); %Choose the time stamp of the *second* "Starting recording" which corresponds to the start of the sleep recording

    end_OF = evts.t{find(contains(evts.label, 'Stopping Recording'))}(1);%Choose the time stamp of the first "Stopping Recording" which corresponds to the end of the OF
    end_sleep = evts.t{find(contains(evts.label, 'Stopping Recording'))}(2);%Choose the time stamp of the second "Stopping Recording" which corresponds to the end of the sleep recording

    %Add the third recoding seesion in case it is a recording from day 2
    if info.session=="D2"
        start_NOPR = evts.t{find(contains(evts.label, 'Starting Recording'))}(3);
        end_NOPR = evts.t{find(contains(evts.label, 'Stopping Recording'))}(3);
    end

    %Printing the duration of OF and sleeping
    fprintf('<strong>OF duration: %.2f mins</strong>\n', (end_OF - start_OF)/60);
    fprintf('<strong>Sleep duration: %.2f mins = %.2f hrs and %.2f min </strong>\n', (end_sleep - start_sleep)/60,floor(((end_sleep - start_sleep)/60)/60), ((end_sleep - start_sleep)/60) -(60*(floor(((end_sleep - start_sleep)/60)/60))));

    if info.session=="D2"
        fprintf('<strong>NOPR duration: %.2f mins</strong>\n', ((end_NOPR - start_NOPR)/60));
    end

    %Restrict the data to just the sleep recording.
    csc_s = restrict(csc, start_sleep, end_sleep);
    emg_s = restrict(emg, start_sleep, end_sleep);

    %Correcting times
    csc_s.tvec= csc_s.tvec-csc_s.tvec(1);
    emg_s.tvec=emg_s.tvec- emg_s.tvec(1);