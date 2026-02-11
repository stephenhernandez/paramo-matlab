function sampleHeader(fname,kstnm,slat,slon,Fs,modelType,tOrig,eqmag,elat,elon,...
    prof,fuente,tRef,duration,Ntot,npoles,lfc,hfc,eventID,cmpTmp)
fName = fullfile('~/Desktop',fname);
fNameID = fopen(fName,'w');

%%
header1 = "Red Nacional de Acelerógrafos (RENAC) Ecuador - IG-EPN";
header2 = "Ladrón de Guevara E11-253, Aptdo. 2759 Quito - Ecuador";
header3 = "Teléfonos: (593-2)2225655 ; (593-2)2225627 Fax: (593-2)2567847";
header4 = "#";
header5 = "Datos de la estación";
header6 = sprintf("Nombre de la estación: %s",kstnm);
header7 = sprintf("Coordenadas de la estación: %3.2f LAT.N, %3.2f LON.E",slat,slon);
header8 = "Tipo de suelo: ";
header9 = "#";
header10 = "Datos del acelerógrafo";
header11 = sprintf("Modelo del acelerógrafo: %s",modelType);
header12 = sprintf("Frecuencia de muestreo: %d",Fs);
header12b = "#";
header13 = "Datos del Sismo";
header14 = sprintf("Fecha/Hora del Sismo (yyyy-mm-dd HH:MM:SS): %s",datestr(tOrig,'yyyy-mm-dd HH:MM:SS.FFF'));
header15 = sprintf("Magnitud del Sismo: %2.1f",eqmag);
header16 = sprintf("Coordenadas del Sismo: %4.3f Lat. %4.3f Lon",elat,elon);
header17 = sprintf("Profundidad del Sismo: %4.1f km",prof);
header18 = sprintf("Fuente de los datos epicentrales: %s",fuente);
header19a = "#";
header19 = "Datos de este Registro";


%
fprintf(fNameID,'#%s\n',header1);
fprintf(fNameID,'#%s\n',header2);
fprintf(fNameID,'#%s\n',header3);
fprintf(fNameID,'%s\n',header4);
fprintf(fNameID,'#%s\n',header5);
fprintf(fNameID,'#%s\n',header6);
fprintf(fNameID,'#%s\n',header7);
fprintf(fNameID,'#%s\n',header8);
fprintf(fNameID,'%s\n',header9);
fprintf(fNameID,'#%s\n',header10);
fprintf(fNameID,'#%s\n',header11);
fprintf(fNameID,'#%s\n',header12);
fprintf(fNameID,'%s\n',header12b);
fprintf(fNameID,'#%s\n',header13);
fprintf(fNameID,'#Codigo Evento: %s\n',eventID);
fprintf(fNameID,'#%s\n',header14);
fprintf(fNameID,'#%s\n',header15);
fprintf(fNameID,'#%s\n',header16);
fprintf(fNameID,'#%s\n',header17);
fprintf(fNameID,'#%s\n',header18);
fprintf(fNameID,'%s\n',header19a);
fprintf(fNameID,'#%s\n',header19);
fprintf(fNameID,'#Componente: %s\n',cmpTmp);
fprintf(fNameID,'#Hora de la primera muestra: %s (UTC)\n',datestr(tRef,'yyyy-mm-dd HH:MM:SS.FFF'));
fprintf(fNameID,'#Exactidud del tiempo: %4.3f [seg.]\n',1/Fs);
fprintf(fNameID,'#Duracion del registro: %g [seg.]\n',duration);
fprintf(fNameID,'#Número total de muestras: %d\n',Ntot);
fprintf(fNameID,'#Unidades: cm/s^2\n');
fprintf(fNameID,'#Tipo de filtro: bandpass butterworth, %d poles, lfc: %g [Hz.], hfc: %g [Hz.]\n',npoles,lfc,hfc);
fprintf(fNameID,'%s\n',"#");
fprintf(fNameID,'%s\n',"#####################################################################################");

% ALJ1, ACUE, AZOG, AMCR, APUY, AMIL, AGRD, GYKA, GYPS, GYGU

%
% fprintf(fNameID,'#Tiempo de referencia del archivo: %s\n',datestr(cutStart,'yyyy-mm-dd HH:MM:SS.FFF'));
% fprintf(fNameID,'#Estacion: %s\n',kstnms(i));
% fprintf(fNameID,'#Componente: %s\n',cmpTmp);
% fprintf(fNameID,'#Frecuencia de muestreo (Hz): %d\n',finalFs);
% fprintf(fNameID,'#Latitud, longitud de la Estacion: %3.2f, %3.2f\n',S1.stla,S1.stlo);
% fprintf(fNameID,'#Tiempo de Origen del Evento: %s\n',datestr(origt,'yyyy-mm-dd HH:MM:SS.FFF'));
% fprintf(fNameID,'#Latitud, longitud del Evento: %3.2f, %3.2f\n',origlat,origlon);
% fprintf(fNameID,'#Magnitud del Evento: %2.1f\n',origmag);
% fprintf(fNameID,'#Unidades: cm/s^2\n');

%
fclose(fNameID);

