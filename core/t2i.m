function i = t2i(t,ref,delta)
%
% t2i time-to-index
%
% i = t2i(t,ref,delta)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
i = 1+floor(seconds(t - ref)/delta);
if i < 1
    i = 1;
end