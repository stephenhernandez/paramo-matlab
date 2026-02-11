function [poles,zeros,sensitivity,normalizationFactor] = readSensorPolesZeros(n)
poles = [];
zeros = [];
sensitivity = [];
normalizationFactor = [];

%% Convention is to convert to displacement (meters), so i added a zero to everything...
if n == 1 %Trillium Compact 120, 749
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0,-434.1+1j*0.0];
    poles = [-0.036910+1j*0.03712,-0.036910-1j*0.03712,-371.2+1j*0.0,-373.9+1j*475.5,...
        -373.9-1j*475.5,-588.4+1j*1508.0,-588.4-1j*1508.0];
    sensitivity = 749.1;
    normalizationFactor = 8.184E+11;
elseif n == 2 % Trillium Compact 120P, 1201
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0,-90.0+1j*0.0,-160.7+1j*0.0,-3108.0+1j*0.0];
    poles = [-0.038520+1j*0.03658,-0.038520-1j*0.03658,-178.0+1j*0.0,...
        -135.0+1j*160.0,-135.0-1j*160.0,-671.0+1j*1154.0,-671.0-1j*1154.0];
    sensitivity = 1201.0;
    normalizationFactor = 308000;
elseif n == 3 % Guralp 40T
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0];
    poles = [-0.011780+1j*0.01178,-0.011780-1j*0.01178,-160.0+1j*0.0,...
        -80.0+1j*0.0,-180.0+1j*0.0];
    sensitivity = 800.0;
    normalizationFactor = 2304000;
elseif n == 4 % Guralp 40T2 (6T! not a "40T2")
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0];
    poles = [-0.023560+1j*0.02356,-62.3816-1j*135.392,-350.0+1j*0.0,-75.0+1j*0.0];
    sensitivity = 1200.0;
    normalizationFactor = 585860000;
elseif n == 5 % Guralp 40Tcuic
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0];
    poles = [-0.023560+1j*0.02356,-0.023650-1j*0.02365,-160.0+1j*0.0,...
        -80.0+1j*0.0,-180.0+1j*0.0];
    sensitivity = 3184.0;
    normalizationFactor = 2304000;
elseif n == 6 % CMG-3ESP
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0];
    poles = [-0.037-1j*0.037000,-0.037000+1j*0.037000,-1130.9733+1j*0.0,...
        -1005.3396+1j*0.0,-502.6548+1j*0.0];
    sensitivity = 2000.0;
    normalizationFactor = 5.71508E+08;
elseif n == 7 % Trillium Compact 120, 754
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0,-3.92e2+1j*0.0,-1.96e3+1j*0.0,-1.49e3+1j*1.74e3,-1.49e3-1j*1.74e3];
    poles = [-0.03691+1j*0.03702,-0.03691-1j*0.03702,-3.43e2+1j*0.0,...
        -370+1j*467.0,-370.0-1j*467.0,-836.0+1j*1522.0,-836.0-1j*1522.0,...
        -4900+1j*4700,-4900-1j*4700,-6900,-1.5e4];
    sensitivity = 754.0;
    normalizationFactor = 4.344928E+17;
elseif n == 8 % Guralp CMG-6T, 30 sec (1200 V/m/s)
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0];
    poles = [-0.02356+1j*0.02356,-0.02356-1j*0.02356,-62.3816+1j*135.392,...
        -62.3816-1j*135.392,-350,-75];
    sensitivity = 1200.0;
    normalizationFactor = 585860000;
elseif n == 9
    zeros = [0.0+1j*0.0, 0.0+1j*0.0, 0.0+1j*0.0, -1.515000e+01+1j*0.0,...
        -1.766000e+02+1j*0.0, -4.631000e+02-1j*4.305000e+02,-4.631000e+02+1j*4.305000e+02];
    poles = [-3.700000e-02-1j*3.700000e-02, -3.700000e-02+1j*3.700000e-02, -1.564000e+01+1j*0.0,...
        -9.734000e+01-1j*4.007000e+02, -9.734000e+01+1j*4.007000e+02, -3.748000e+02+1j*0.0,...
        -5.203000e+02+1j*0.0, -1.053000e+04-1j*1.005000e+04, -1.053000e+04+1j*1.005000e+04,...
        -1.330000e+04+1j*0.0,-2.550970e+02+1j*0.0];
    sensitivity = 1491.5; %2502320000 = 1491.5 * 1677720;
    normalizationFactor = 3.5432e+17;
elseif n == 10 % wood-anderson response
    %[zeros,poles,sensitivity] = read_sac_pole_zero('~/regions/wa2.pz');
    zeros = [0+1j*0,0.0+1j*0.0];
    poles = [-5.49779+1j*5.60886, -5.49779 - 1j*5.60886];
    sensitivity = 2080;
    normalizationFactor = 1;
elseif n == 11 % anatahan sts2
    zeros = [0+1j*0,0.0+1j*0.0,0.0+1j*0.0];
    poles = [-0.03701+1j*0.03701,-0.03701-1j*0.03701,-251.3+1j*0,...
        -131+1j*467.3,-131-1j*467.3];
    normalizationFactor = 4.655677e+16;
    sensitivity = 1;
elseif n == 12
    zeros = [0+1j*0,0.0+1j*0.0,0.0+1j*0.0,-113.349];
    poles = [-0.0444288+1j*0.0444288,-0.0444288-1j*0.0444288,...
        -394.141+1j*394.141,-394.141-1j*394.141,-106.814,...
        -2544.69+1j*1232.45,-2544.69-1j*1232.45,-39478.4];
    normalizationFactor = 9.24047e+16;
    sensitivity = 6.237;
elseif n == 13 %L4C??
    zeros = [0.0+1j*0.0,0.0+1j*0.0,0.0+1j*0.0];
    poles = [-4.44+1j*4.44,-4.44-1j*4.44];
    normalizationFactor = 1;
    sensitivity = 1;
elseif n == 14 %geotech ks-54000 from PAYG, post-2013
    zeros = [0.0+1j*0.0, 0.0+1j*0.0, 0.0+1j*0.0];
    poles = [-5.943130e+01 + 1j*0.000000e+00; ...
        -2.271210e+01 + 1j*2.710650e+01; ...
        -2.271210e+01 - 1j*2.710650e+01; ...
        -4.800400e-03 + 1j*0.000000e+00; ...
        -7.384400e-02 + 1j*0.000000e+00];
    sensitivity = 3.337050e09;
    normalizationFactor = 8.627050e04;
elseif n == 15 %geotech ks-54000 from PAYG, 1999-2009
    zeros = [0.0+1j*0.0, 0.0+1j*0.0, 0.0+1j*0.0];
    poles = [-5.943130e+01 + 1j*0.000000e+00; ...
        -2.271210e+01 + 1j*2.710650e+01; ...
        -2.271210e+01 - 1j*2.710650e+01; ...
        -4.800400e-03 + 1j*0.000000e+00; ...
        -7.384400e-02 + 1j*0.000000e+00];
    sensitivity = 8.561000e08;
    normalizationFactor = 8.608300e04;
elseif n == 16 %smartsolo
    poles = [-2.221110e+01-1j*2.221780e+01;...
        -2.221110e+01+1j*2.221780e+01];
    zeros = [0.000000E+00+1j*0.000000E+00;...
        0.000000E+00+1j*0.000000E+00];
    sensitivity = 1;
    normalizationFactor = 1;
end
