function digitizerGain = readDigitizerGain(n)
% see also: "readSensorPolesAndZeros"
if n == 1 %Q330S
    digitizerGain = 419430;
elseif n == 2 %Q330SR
    digitizerGain = 419430;
elseif n == 3 %Smart 24 JICA
    digitizerGain = 609014;
elseif n == 4 %Q330HRS
    digitizerGain = 1677720;
elseif n == 5 %CD24-C2913 (SRAM)
    digitizerGain = 1037703;
elseif n == 6 %RT130
    digitizerGain = 629121;
elseif n == 7 %DM24
    digitizerGain = 348432;
elseif n == 8 %CD24-C2917 (TOMA)
    digitizerGain = 1006711;
elseif n == 9 %CD24-C2919 (BRRN)
    digitizerGain = 1033058;
elseif n == 10 %KEPHREN
    digitizerGain = 1865670;
elseif n == 11 %CD24-C1039 (SUCR)
    digitizerGain = 1041305;
elseif n == 12 %PAYG.2010.180.H10
    digitizerGain = 1677720; %<-- this is just a guess
elseif n == 13 %PVIL
    digitizerGain = 314660; %<-- this is another guess
else
    digitizerGain = 1;
end