%%
%SPECTRA_preprocessing.
%
%Preprocessing settings for Frontal theta in PST task. 
%Adapted by Dr. Daniel Robles, Baker Lab, RU. 
%Aug 2023 
%

%% 
%NOTES
%Fixation/trial start = S1

%Stim type 
%A (80%) || B (20%), hence S12 = AB, S21 = BA
%C (70%) || %D (30%), hence S34 = CD, S43 = DC 

%response type
% S4 = left || S5 = right 
% 
% Feedback type 
% S6 = CORRECT left || S7 = INCORRECT left 
% S8 = Correct right || S9 = INCORRECT right 

%S77 = no response detected
%%

%Misc stuff 
%%%006 onward has 'normal' trigger
%%%001:004 has odd triggers

% function [NewMarkers] = PST_update_feedback_marker (Markers)
% 
% Markers_copy = Markers;
% 
% for i = 1:size(Markers,2)
%     
%     if strcmp(Markers(1,i).Description, 'S  6') == 1 ||  strcmp(Markers(1,i).Description, 'S  7') == 1 || strcmp(Markers(1,i).Description, 'S  8') == 1 ||  strcmp(Markers(1,i).Description, 'S  9') == 1
%        if strcmp(Markers(1,i-2).Description , 'S 21') == 1 || strcmp(Markers(1,i-2).Description , 'S 12') == 1 
%         Markers(1,i).Description = [Markers(1,i).Description,'_ab'];
%        
%        elseif strcmp(Markers(1,i-2).Description , 'S 34') == 1 || strcmp(Markers(1,i-2).Description , 'S 43') == 1 
%         Markers(1,i).Description = [Markers(1,i).Description,'_cd'];
% 
%        end
%     end
% 
%     if strcmp(Markers(1,i).Description, 'S201') == 1
%        Markers(1,1:i) = Markers_copy(1,1:i);
%     end
% 
% end
% 
% [NewMarkers] = Markers;

ccc
parts = {'0017'};
session = {'S03'};
tmstarget = {'T03'};
%conditions = {'positive';'negative'};
filepath = ['C:\Users\Kathryn\Documents\PST_DR\raw_files'];% C:\Users\Kathryn\Documents\PST_DR\
exp = 'PST';
i_sess = 1:length(session); 
i_targ = 1:length(tmstarget);

%%%old settings for selection cards%%%
% % selection_cards_far = {'1 2 3 4', '5', '21 22 23 24 25'};
% % selection_cards_near = {'1', '2 3 4 5', '21 22 23 24 25'};


eeg_thresh_1 = [-500,500];
eeg_thresh_2 = [-200,200];

if ~exist([filepath 'segments\']) %when preprocessing files, puts everything ina folder called segments, if not there, makes it
    mkdir([filepath 'segments\']);
end

high_pass = 0.1;
low_pass = 40;

epoch_size = [-2.5  2.5];
baseline = [-200    0];
eeg_thresh_1 = [-1000,1000];
eeg_thresh_2 = [-500,500];

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% i_part = 10;
% i_cond = 1;

%%
for i_part = 1:length(parts)
    
    for i_sess = 1:length(session)
        for i_targ = 1:length (tmstarget)
            
            events = {'1' '2' '3' '4'};
            %         selection_cards = {'1' '2'};
            selection_cards = {'1' '2'};
            
            %%%load data%%%
            disp(['Processing data for participant ' parts{i_part}]);
            disp(['Loading file: ' exp session{i_sess} tmstarget{i_targ} parts{i_part} '.vhdr']);
            %         filename = [parts{i_part} '_' exp '_' conditions{i_cond} '.vhdr'];
            %         EEG = pop_loadbv(filepath, filename, [], []); % loading in EEG data
            %         savename = [parts{i_part} '_' exp '_' conditions{i_cond}];
            %exp-sessopm-tmstarg-part
            filename = [exp session{i_sess} tmstarget{i_targ} '_' parts{i_part} '.vhdr'];
            EEG = pop_loadbv(filepath, filename, [], []); % loading in EEG data
            savename = [exp session{i_sess} tmstarget{i_targ} parts{i_part}];
            % get electrode locations
            EEG=pop_chanedit(EEG, 'load',{'C:\Users\Kathryn\Documents\PST_DR\stms.ced' 'filetype' 'autodetect'});
            
            % arithmetically rereference to linked mastoid
            for x=1:EEG.nbchan-2
                EEG.data(x,:) = (EEG.data(x,:)-((EEG.data(EEG.nbchan-2,:))*.5));
            end
            
            %Filter the data with low pass of 30
            EEG = pop_eegfilt( EEG, high_pass, 0, [], 0);  %high pass filter
            EEG = pop_eegfilt( EEG, 0, low_pass, [], 0);  %low pass filter
            
            
            all_events = length(EEG.event);
            
            
            for i_event = 2:all_events %%% Why start from 2
                if strcmp(EEG.event(i_event).type, 'boundary')
                    continue
                else
                    EEG.event(i_event).type = num2str(str2num((EEG.event(i_event).type(2:end))));
                end
            end
            
            %epoch
            
            EEG = pop_epoch( EEG, events, epoch_size, 'newname',  sprintf('%s epochs' , exp), 'epochinfo', 'yes'); %Changed from [-.2 1] to [-1 2]. DR
            %S1 is standard, S2 = target....
            EEG = pop_rmbase( EEG, baseline);
            
            %    Artifact rejection, trials with range >500 uV
            EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],eeg_thresh_1(1),eeg_thresh_1(2),EEG.xmin,EEG.xmax,0,1);
            
            %   EMCP occular correction
            temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
            EEG = gratton_emcp(EEG,selection_cards,{'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
            EEG.emcp.table %this prints out the regression coefficients
            EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
            
            %    Artifact rejection, trials with range >250 uV
            EEG = pop_rmbase( EEG, baseline); %baseline again since this changed it
            EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],eeg_thresh_2(1),eeg_thresh_2(2),EEG.xmin,EEG.xmax,0,1);
            
            EEG_Copy = EEG;
            
            EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{1}),'renametype','Targets','deleteevents','on','deleteepochs','on','invertepochs','off');
            EEG = pop_editset(EEG, 'setname',[savename '_fft_Targets']);
            EEG = pop_saveset(EEG, 'filename',[savename '_fft_Targets'],'filepath',[filepath 'segments_fft_JK\']);
            
        end
    end
end

%%

