function [zz,pp,constant,knetwk,kstnm,khole,kcmpnm,tStart,tEnd,stla,stlo,stel,stde,...
    stdi,staz,Fs,instType,instGain,sensitivity,a0] = read_sac_pole_zero(pole_zero_file_name)

% function to read a sac pole zero file
% usage:
%   [zeros, poles, constant] = read_sac_pole_zero(pole_zero_file_name)
%   creates two output vectors (poles, zeros) and one scalar (constant)
% file is the format as output from rdseed if the option
% "output polezero file" is accepted
% format:
% ZEROS nzeros (only listed if not at 0+0i - this version only allows zeros at the origin!)
% zero_r zero_i
% zero_r zero_i
% ...    ...
% POLES npoles
% pole_r pole_i
% pole_r pole_i
% ...     ...
% CONSTANT constant
% example:
%    ZEROS 3
%    POLES 5
%   -0.0370  0.0370
%   -0.0370  -0.0370
%   -118.7520  423.4880
%   -118.7520  -423.4880
%   -251.3270  0.0000
%   CONSTANT 3.127953e+16

%
% Written by: RWP 5.30.2011
% added functionality to allow #, *, or % as comments in sac pole zero file
% Modified by Stephen Hernandez, 03-09-2015
%

%% First open the file
pz_fid = fopen(pole_zero_file_name);

%% define defaults
zz = 0;
pp = 0;
constant = 0;
knetwk = '';
kstnm = knetwk;
khole = knetwk;
kcmpnm = knetwk;

%% initialize some flags and various variables
buff = 0;
nzeros = 0;
npoles = 0;
zero_scan_flag = -1;
pole_scan_flag = -1;
pole_count_flag = 1;
zero_count_flag = 1;

%% loop over the entire file
while (buff ~= -1)
    % read the next line in the file
    buff = fgets(pz_fid);
    % check to make sure it isn't the end of the file
    if (buff ~= -1)
        %assume first line is ZEROS nzeros
        %tmp = textscan(buff, '%s" "%d');
        tmp = textscan(buff, '%s %d');
        tmp2 = textscan(buff, '%s %s');
        tmp = tmp{1};
        tmp2 = tmp2{2};
        if (strcmp(tmp,'*'))
            if strcmp(tmp2{1},'NETWORK')
                C = textscan(buff,'%s %s %s %s');
                knetwk = string(char(C{4}));
            elseif strcmp(tmp2{1},'STATION')
                C = textscan(buff,'%s %s %s %s');
                kstnm = string(C{4});
            elseif strcmp(tmp2{1},'LOCATION')
                C = textscan(buff,'%s %s %s %s');
                khole = string(C{4});
            elseif strcmp(tmp2{1},'CHANNEL')
                C = textscan(buff,'%s %s %s %s');
                kcmpnm = string(C{4});
            elseif strcmp(tmp2{1},'START')
                C = textscan(buff,'%s %s %s %s');
                timeStr = char(C{4});
                yyyy_ = str2double(timeStr(1:4));
                month_ = str2double(timeStr(6:7));
                day_ = str2double(timeStr(9:10));
                hour_ = str2double(timeStr(12:13));
                min_ = str2double(timeStr(15:16));
                sec_ = str2double(timeStr(18:19));
                tStart = datetime([yyyy_ month_ day_ hour_ min_ sec_]);
            elseif strcmp(tmp2{1},'END')
                C = textscan(buff,'%s %s %s %s');
                timeStr = char(C{4});
                yyyy_ = str2double(timeStr(1:4));
                month_ = str2double(timeStr(6:7));
                day_ = str2double(timeStr(9:10));
                hour_ = str2double(timeStr(12:13));
                min_ = str2double(timeStr(15:16));
                sec_ = str2double(timeStr(18:19));
                tEnd = datetime([yyyy_ month_ day_ hour_ min_ sec_]);
            elseif strcmp(tmp2{1},'LATITUDE')
                C = textscan(buff,'%s %s %s %f');
                    stla = C{4};
            elseif strcmp(tmp2{1},'LONGITUDE')
                C = textscan(buff,'%s %s %s %f');
                stlo = C{4};
            elseif strcmp(tmp2{1},'ELEVATION')
                C = textscan(buff,'%s %s %s %f');
                stel = C{4};
            elseif strcmp(tmp2{1},'DEPTH')
                C = textscan(buff,'%s %s %s %f');
                stde = C{4};
            elseif strcmp(tmp2{1},'DIP')
                C = textscan(buff,'%s %s %s %f');
                stdi = C{4};
            elseif strcmp(tmp2{1},'AZIMUTH')
                C = textscan(buff,'%s %s %s %f');
                staz = C{4};
            elseif strcmp(tmp2{1},'SAMPLE') %sample rate
                C = textscan(buff,'%s %s %s %s %d');
                Fs = C{5};
            elseif strcmp(tmp2{1},'INSTTYPE')
                C = textscan(buff,...
                    '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'); %make up a max number of white-delimited comments...
                C = string(C);

                instType = C(4:end);
                badI = instType == "";
                instType(badI) = [];
                lType = length(instType);
                if lType > 1
                    str1 = instType(1);
                    for j = 2:lType
                        str1 = append(str1," ",instType(j));
                    end
                    instType = str1;
                elseif lType < 1
                    instType = "";
                end
            elseif strcmp(tmp2{1},'INSTGAIN')
                C = textscan(buff,'%s %s %s %f %s');
                instGain = C{4};
            elseif strcmp(tmp2{1},'SENSITIVITY')
                C = textscan(buff,'%s %s %s %f %s');
                sensitivity = C{4};
                sensitivity = sensitivity(1);
            elseif strcmp(tmp2{1},'A0')
                C = textscan(buff,'%s %s %s %f');
                a0 = C{4};
            end
        elseif (strcmp(tmp,'#'))
        elseif (strcmp(tmp,'%'))
        elseif strcmp(tmp, 'ZEROS')
            zero_scan_flag = 1;
            pole_scan_flag = 0;
            tmp = sscanf(buff, '%s %d');
            nzeros = tmp(6);
        elseif strcmp(tmp, 'POLES')
            pole_scan_flag = 1;
            zero_scan_flag = 0;
            tmp = sscanf(buff,'%s %d');
            npoles = tmp(6);
        elseif strcmp(tmp, 'CONSTANT')
            pole_scan_flag = 0;
            tmp = sscanf(buff, '%s %e');
            constant = tmp(9);
        elseif (zero_scan_flag == 1)
            if (zero_count_flag < nzeros+1)
                tmp = sscanf(buff,'%f %f');
                zero_count_flag = zero_count_flag + 1;
                if (zero_count_flag > 1)
                    j=zero_count_flag-1;
                    zz(j) = tmp(1) + tmp(2)*1j;
                end
            end
        elseif (pole_scan_flag == 1)
            if (pole_count_flag < npoles+1)
                tmp = sscanf(buff, '%f %f');
                pole_count_flag = pole_count_flag + 1;
                if (pole_count_flag > 1)
                    j=pole_count_flag-1;
                    pp(j) = tmp(1) + tmp(2)*1j;
                end
            end
        end
    end
end

%% fill in missing poles and zeros
if (length(zz) < nzeros)
    for j=length(zz)+1:nzeros
        zz(j) = 0+0i;
    end
end

if (length(pp) < npoles)
    for j=length(pp)+1:npoles
        pp(j) = 0+0i;
    end
end
fclose(pz_fid);

if instGain == 0
    instGain = 1;
end

if ~nzeros
    zz = [];
end