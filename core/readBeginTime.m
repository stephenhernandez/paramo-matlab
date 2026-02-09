function dn = readBeginTime(beginTime,isBigEndian)
%
% d = readBeginTime(beginTime,isBigEndian)
%
% reads "btime" structure and returns
% D = [YEAR,DAYOFYEAR,HOUR,MINUTE,SECONDS,Seconds0001]
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019

%%
Year					= typecastArray(beginTime(1:2,:),"uint16");    %21:22
DayOfYear				= typecastArray(beginTime(3:4,:),"uint16");    %23:24
Hours					= typecastArray(beginTime(5,:),"uint8");       %25
Minutes					= typecastArray(beginTime(6,:),"uint8");       %26
Seconds					= typecastArray(beginTime(7,:),"uint8");       %27
Seconds0001				= typecastArray(beginTime(9:10,:),"uint16");   %29:30

%%
if isBigEndian
    Year = swapbytes(Year);
    DayOfYear = swapbytes(DayOfYear);
    Hours = swapbytes(Hours);
    Minutes = swapbytes(Minutes);
    Seconds = uint16(swapbytes(Seconds)); %swap, then cast
    Seconds0001 = swapbytes(Seconds0001);
end

%%
mm = size(Year);
mm = ones(mm);
dn = datenum([double(Year),mm,double(DayOfYear),double(Hours),...
    double(Minutes),double(Seconds) + double(Seconds0001)*1e-4]);
