function mw = m02mw(m0,dir)
if nargin < 2; dir = true; end
mw = (2/3)*(log10(m0)-9.1);
%mw = round(mw*10)/10;

if ~dir
    disp('mw2m0');
    mw = 1.5*m0 + 9.1;
    mw = 10.^mw;
end