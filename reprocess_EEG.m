% reprocess_EEG.m

function reprocess_EEG(filenamestem)
% Load preprocessed data
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',[filenamestem,'_filt+resamp200+music.set'],'filepath','/Users/jthompson/data/EEG-Thesis/');

% % Remove extra channels 
% EEG = pop_select( EEG,'channel',{'Fp1' 'AF7' 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'Po3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'AFz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'Po4' 'O2'});
% 
% % Resample to 200 Hz 
% EEG = pop_resample( EEG, 200);
% % save resampled
% EEG.setname=[filenamestem,'_resamp200'];
% EEG = pop_saveset( EEG, 'filename',[filenamestem,'_resamp200.set'],'filepath','/Users/jthompson/data/EEG-Thesis/');
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% 
% % Butterworth IIR filter with cutoff freqs at .1 Hz and 55Hz (less than
% % 60Hz to avoid noise from electical devices, could go down to 30Hz)
% EEG  = pop_basicfilter( EEG,  1:64 , 'Cutoff', [ 0.1 55], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2 );
% 
% % Save filtered and resampled
% EEG.setname=[filenamestem,'_resamp200_filt'];
% EEG = pop_saveset( EEG, 'filename',[filenamestem,'_resamp200+filt.set'],'filepath','/Users/jthompson/data/EEG-Thesis/');
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% 
% % Reject non music segments
% num_events = size(EEG.event, 2);
% regions = ones((num_events/2)+1,2);
% regions(end,2) = size(EEG.data, 2);
% j=1;
% for i=1:2:num_events-1
%     EEG.event(i).duration = EEG.event(i+1).latency - EEG.event(i).latency;
%     regions(j,2) = EEG.event(i).latency-2; %end of section to remove, beginning of music
%     regions(j+1,1) = EEG.event(i+1).latency+2; %begining of next section to remove, end of music
%     j=j+1;
% end
% % reject regions
% EEG = eeg_eegrej( EEG, regions );
% 
% % Save as EEGLAB data set which generates a .set file and a .fdt file
% %EEG = eeg_checkset( EEG );
% EEG.setname=[filenamestem,'_filt+resamp200+music'];
% EEG = pop_saveset( EEG, 'filename',[filenamestem,'_filt+resamp200+music.set'],'filepath','/Users/jthompson/data/EEG-Thesis/');
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Run ICA
% Only include first 64 channels 
%EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Import sensor location
EEG=pop_chanedit(EEG, 'lookup','/Users/jthompson/matlab/eeglab12_0_1_0b/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp');

% dipole fitting
EEG = pop_dipfit_settings( EEG, 'hdmfile','/Users/jthompson/matlab/eeglab12_0_1_0b/plugins/dipfit2.2/standard_BEM/standard_vol.mat','coordformat','MNI','mrifile','/Users/jthompson/matlab/eeglab12_0_1_0b/plugins/dipfit2.2/standard_BEM/standard_mri.mat','chanfile','/Users/jthompson/matlab/eeglab12_0_1_0b/plugins/dipfit2.2/standard_BEM/elec/standard_1005.elc','coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485] ,'chansel',[1:64] );
EEG = pop_multifit(EEG, [1:64] ,'threshold',20,'dipplot','on','plotopt',{'normlen' 'on'});

% Save ICA + dipole
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%EEG = eeg_checkset( EEG );
EEG.setname=[filenamestem,'_filt+resamp200+music+ICA'];
EEG = pop_saveset( EEG, 'filename',[filenamestem,'_filt+resamp200+music+ICA.set'],'filepath','/Users/jthompson/data/EEG-Thesis/');

% clear dataset
ALLEEG = pop_delset(ALLEEG, 1);

% Plot 2D ICA spatial maps
% EEG = eeg_checkset( EEG );
% f=figure;pop_topoplot(EEG,0, [1:64] ,[filenamestem,'_filt+resamp256+music'],[8 8] ,0,'electrodes','off');
% eeglab redraw;
% saveas(f, [filenamestem, '_ICA.fig'],'fig')
% saveas(f, [filenamestem, '_ICA.png'],'png')
% 
% % ICA Spect plots
% % 8.6 Hz
% screen_size = get(0, 'ScreenSize');
% EEG = eeg_checkset( EEG );
% f=figure; pop_spectopo(EEG, 0, [0      3656719.0601], 'EEG' , 'freq', [8.6], 'plotchan', 0, 'percent', 70, 'icacomps', [1:64], 'nicamaps', 5, 'freqrange',[1 35],'electrodes','off');
% set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% saveas(f, [filenamestem, '_ICAspect8.6.fig'],'fig')
% saveas(f, [filenamestem, '_ICAspect8.6.png'],'png')
% % 10.5 Hz
% EEG = eeg_checkset( EEG );
% f=figure; pop_spectopo(EEG, 0, [0      3656719.0601], 'EEG' , 'freq', [10.5], 'plotchan', 0, 'percent', 70, 'icacomps', [1:64], 'nicamaps', 5, 'freqrange',[1 35],'electrodes','off');
% set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% saveas(f, [filenamestem, '_ICAspect10.5.fig'],'fig')
% saveas(f, [filenamestem, '_ICAspect10.5.png'],'png')
% % 21.4 Hz
% EEG = eeg_checkset( EEG );
% f=figure; pop_spectopo(EEG, 0, [0      3656719.0601], 'EEG' , 'freq', [21.4], 'plotchan', 0, 'percent', 70, 'icacomps', [1:64], 'nicamaps', 5, 'freqrange',[1 35],'electrodes','off');
% set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% saveas(f, [filenamestem, '_ICAspect21.4.fig'],'fig')
% saveas(f, [filenamestem, '_ICAspect21.4.png'],'png')
