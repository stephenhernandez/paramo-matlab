% april29_2024
clear; close all;
cd ~/data/tmp
load tmp_sangay_28APR2024
raw_shifts = raw_shifts/12;

%%
load tais_tmp_waveforms.mat;
%load pkyu_tmp_waveforms.mat;
%load bulb_tmp_waveforms.mat;
%load puyo_tmp_waveforms.mat;
C = C(:,1:3);
goodI = sum(isnat(pull(C,'ref')),2) == 0;
C = C(goodI,:);
nShifts = sum(goodI);
for i = 1:3
    C_ = C(:,i);
    Z = pull(C_);
    Z = apply_shifts(Z,raw_shifts);
    nShifts = size(raw_shifts);
    for j = 1:nShifts
        newFs = 1./C_(j).delta;
        newRef = C_(j).ref;
        oldRef = newRef;
        newRef = newRef - seconds(raw_shifts(j)/newFs);
        C_(j) = dealHeader(C_(j),Z(:,j),newFs,newRef);
    end
    C(:,i) = C_;
end

%%
npts = 1800;
Z = NaN([npts*3,nShifts]);
for j = 1:nShifts
    C_ = C(j,:);
    Z_ = pull(C_);
    Z_ = Z_(1:npts,:);
    Z_ = Z_';
    Z_ = Z_(:);
    Z(:,j) = Z_;
end
Z(~isfinite(Z)) = 0;
[U,E,~] = svd(normalizeWaveforms(Z),'econ');

%%
% U = populateWaveforms([20,4]);
% for i = 1:Nsub
% U_ = U(i,1);
% U_ = dealHeader(U_,UsubCO1V(:,i),20,datetime(1970,01,01));
% U_.kstnm = "CO1V";
% U_.kcmpnm = "HHZHHNHHE";
% U_.knetwk = "EC";
% U(i,1) = U_;
% end
% U(1)
% for i = 1:Nsub
% U_ = U(i,2);
% U_ = dealHeader(U_,UsubBREF(:,i),20,datetime(1970,01,01));
% U_.kstnm = "BREF";
% U_.kcmpnm = "BHZBHNBHE";
% U_.knetwk = "EC";
% U(i,2) = U_;
% end
% U(1)
% U(1,2)
% for i = 1:Nsub
% U_ = U(i,3);
% U_ = dealHeader(U_,UsubBVC2(:,i),20,datetime(1970,01,01));
% U_.kstnm = "BVC2";
% U_.kcmpnm = "BHZBHNBHE";
% U_.knetwk = "EC";
% U(i,3) = U_;
% end
% U(1,3)
% for i = 1:Nsub
% U_ = U(i,4);
% U_ = dealHeader(U_,UsubBTAM(:,i),20,datetime(1970,01,01));
% U_.kstnm = "BTAM";
% U_.kcmpnm = "BHZBHNBHE";
% U_.knetwk = "EC";
% U(i,4) = U_;
% end
