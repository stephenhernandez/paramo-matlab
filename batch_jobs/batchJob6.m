clear; close all;
try
    SAG1_triggeredArrayProcessing();
    clearvars -except *Main lfc hfc secDur ampThresh thresh Fs Fs2
    save('SAG1_infrasoundArrayCatalog_4Sec1Hz.mat');
catch
    clearvars -except *Main lfc hfc secDur ampThresh thresh Fs Fs2
    save('SAG1_infrasoundArrayCatalog_4Sec1Hz.mat');
end