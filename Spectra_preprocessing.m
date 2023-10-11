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

... I need positive ab, negative ab,
    ... positive cd, negative cd
    %%

close all
clear all

parts = {'0017'};
session = {'S03'};
tmstarget = {'T03'};
%conditions = {'positive';'negative'};
filepath = 'C:\Users\Kathryn\Documents\PST_DR\raw_files\';% C:\Users\Kathryn\Documents\PST_DR\
exp = 'PST';
% i_sess = 1:length(session);
% i_targ = 1:length(tmstarget);

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
ref = [8, 16];

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% i_part = 10;
% i_cond = 1;

%%
for i_part = 1:length(parts)

    for i_sess = 1:length(session)
        for i_targ = 1:length (tmstarget)

            % events = {'12' '21' '34' '43'};
            % selection_cards = {'12' '21' '34' '43'};
            events = {'6'  '7'  '8'  '9'};
            selection_cards = {'6'  '7'  '8'  '9'};

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
            %EEG=pop_chanedit(EEG, 'load',{'C:\Users\Kathryn\Documents\PST_DR\stms.ced' 'filetype' 'autodetect'});

            % arithmetically rereference to linked mastoid
            % for x=1:EEG.nbchan-2
            %     EEG.data(x,:) = (EEG.data(x,:)-((EEG.data(EEG.nbchan-2,:))*.5));
            % end


            EEG = pop_reref(EEG, ref);
            %Filter the data with low pass of 30
            EEG = pop_eegfilt( EEG, high_pass, 0, [], 0);  %high pass filter
            EEG = pop_eegfilt( EEG, 0, low_pass, [], 0);  %low pass filter


            all_events = length(EEG.event);


            for i_event = 2:all_events 
                if strcmp(EEG.event(i_event).type, 'boundary')
                    continue
                else
                    EEG.event(i_event).type = num2str(str2num((EEG.event(i_event).type(2:end))));
                end
            end

            %epoch

            EEG = pop_epoch( EEG, events, epoch_size, 'newname',  sprintf('%s epochs' , exp), 'epochinfo', 'yes');
            EEG = pop_rmbase( EEG, baseline);

            %    Artifact rejection, trials with range >500 uV
            EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],eeg_thresh_1(1),eeg_thresh_1(2),EEG.xmin,EEG.xmax,0,1);

            %   EMCP occular correction
            temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
            %EEG = gratton_emcp(EEG,selection_cards,{'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
            EEG = gratton_emcp(EEG,selection_cards,{'VEOG'}); %this assumes the eye channels are called this
            EEG.emcp.table %this prints out the regression coefficients
            EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data

            %    Artifact rejection, trials with range >250 uV
            EEG = pop_rmbase( EEG, baseline); %baseline again since this changed it
            EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],eeg_thresh_2(1),eeg_thresh_2(2),EEG.xmin,EEG.xmax,0,1);

            EEG_Copy = EEG;
            %
            % EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{1}),'renametype','Targets','deleteevents','on','deleteepochs','on','invertepochs','off');
            % EEG = pop_editset(EEG, 'setname',[savename '_fft_Targets']);
            % EEG = pop_saveset(EEG, 'filename',[savename '_fft_Targets'],'filepath',[filepath 'segments\']);

            EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{1}),'renametype','Pos_left','deleteevents','on','deleteepochs','on','invertepochs','off');
            EEG = pop_editset(EEG, 'setname',[savename '_Pos_left']);
            EEG = pop_saveset(EEG, 'filename',[savename '_Pos_left'],'filepath',[filepath 'segments\']);

            EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{2}),'renametype','Neg_left','deleteevents','on','deleteepochs','on','invertepochs','off');
            EEG = pop_editset(EEG, 'setname',[savename '_Neg_left']);
            EEG = pop_saveset(EEG, 'filename',[savename '_Neg_left'],'filepath',[filepath 'segments\']);

            EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{3}),'renametype','Pos_right','deleteevents','on','deleteepochs','on','invertepochs','off');
            EEG = pop_editset(EEG, 'setname',[savename '_Pos_right']);
            EEG = pop_saveset(EEG, 'filename',[savename '_Pos_right'],'filepath',[filepath 'segments\']);

            EEG = pop_selectevent(EEG_Copy, 'type',str2num(events{4}),'renametype','Neg_right','deleteevents','on','deleteepochs','on','invertepochs','off');
            EEG = pop_editset(EEG, 'setname',[savename '_Neg_right']);
            EEG = pop_saveset(EEG, 'filename',[savename '_Neg_right'],'filepath',[filepath 'segments\']);

        end
    end
end





