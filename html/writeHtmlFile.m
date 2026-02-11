function writeHtmlFile(SNCL,fName)
% cd ~/public_html/helis/v2/
% fName = fullfile('~','public_html','helis','v2',strcat(SNCL,'_1Hz8Hz.html'));
fNameID = fopen(fName,'w');

fprintf(fNameID,"<!DOCTYPE html>\n");
fprintf(fNameID,"<html>\n");
fprintf(fNameID,"<head>\n");
fprintf(fNameID,"<meta charset='UTF-8'>\n");
fprintf(fNameID,"<meta http-equiv='X-UA-Compatible' content='IE=edge'>\n");
fprintf(fNameID,"<meta name='viewport' content='width=device-width,initial-scale=1.0'>\n");
fprintf(fNameID,"<link rel='stylesheet' href='css/style.css'>\n");
fprintf(fNameID,"<title>Helicorder %s</title>\n",SNCL);
fprintf(fNameID,"</head>\n");
fprintf(fNameID,"<body>\n");
fprintf(fNameID,"<header>\n");
fprintf(fNameID,"<h1>\n");
fprintf(fNameID,"%s 1 - 8 Hz.\n",SNCL);
fprintf(fNameID,"</h1>\n");
fprintf(fNameID,"</header>\n");
fprintf(fNameID,"<main>\n");
fprintf(fNameID,"<p>\n");
fprintf(fNameID,"<b>\n");
fprintf(fNameID,"*El uso de este sitio es para uso exclusivo de los cientificos del Instituto Geofisico de la Escuela Politecnica Nacional y/o sus colegas invitados.\n");
fprintf(fNameID,"</b>\n");
fprintf(fNameID,"</p>\n");
jpgName = fullfile('..','assets','jpg',strcat(SNCL,'_1Hz8Hz.jpg'));
fprintf(fNameID,"<a href='%s' target='_blank'>\n",jpgName);
fprintf(fNameID,"<img src='%s' alt='%s, Filtrado'>\n",jpgName,SNCL);
fprintf(fNameID,"</a>\n");
fprintf(fNameID,"</main>\n");
fprintf(fNameID,"<footer>\n");
fprintf(fNameID,"<div class='author'>\n");
fprintf(fNameID,"by: stephen hernandez\n");
fprintf(fNameID,"</div>\n");
fprintf(fNameID,"</footer>\n");
fprintf(fNameID,"</body>\n");
fprintf(fNameID,"</html>\n");
fclose(fNameID);