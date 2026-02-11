function newvector = typecastArray(theArray,newType)
%
% newArray = typecastArray(theArray,newType)
%
% the only thing that this code seems to be doing is typecasting but by
% having converted the input array into a vector first
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019

%%
newvector = typecast(theArray(:),newType);