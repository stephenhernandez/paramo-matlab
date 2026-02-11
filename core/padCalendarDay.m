function dayStr = padCalendarDay(day)
% DEPRECATED
% day = padCalendarDay(day)
%
% padCalendarDay num2str from doy with leading zeros if necessary
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
dayStr = sprintf("%03d",day);