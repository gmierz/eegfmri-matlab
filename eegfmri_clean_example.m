%% Example using eegfmri_clean

% Open the data
PATH_TO_DATA = ''
DATASET = ''

cd(PATH_TO_DATA)
EEG = pop_loadbv('.', DATASET);

% Assuming you have the script in the PATH, clean the EEG data
cleaneeg = eegfmri_clean(EEG, 'test_data.set');
