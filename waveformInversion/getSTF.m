function stf = getSTF(th,Fs)
% Symmetry is assumed
%th = th - 1;
if th == 0
    stf = 1;
elseif th > 0
    td = 2*th*Fs + 1;
    stf = window(@tukeywin,td,1); %
    %stf = window(@triang,td);
    stf = stf - min(stf);
    %stf = window(@tukeywin,td,1);
    stf = stf/trapz(stf);
else
    disp('half duration (th) must be an integer greater than 0')
    return
end
end