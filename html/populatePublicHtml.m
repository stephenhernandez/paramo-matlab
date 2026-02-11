function populatePublicHtml(masterDirectory)
if nargin < 1; masterDirectory = 'ecuadorian-seismicity'; end

%%
regions = ["ecuador-continental","costa","galapagos","cotopaxi",...
    "pisayambo","guagua","tungurahua","quito","reventador","cuicocha",...
    "pululahua","esmeraldas","sierra-negra","cerro-azul","punta-galera",...
    "cayambe","chimborazo","quilotoa","sangay","fernandina","jama",...
    "puna","plata","chiles","antisana","alcedo","wolf","manta",...
    "bahia-elizabeth","darwin","outer-rise","balao","pedernales"];
mainIndex(masterDirectory,regions,'index.html');

%%
[currentYear,currentMonth,currentDay] = datevec(today);
years = (2011:currentYear)';
yyyyStr = string(num2str(years));
months = (01:12)';
monthNumStr = [];

for i = 1:length(months)
    if months(i) < 10
        monthNumStr = [monthNumStr; ['0',num2str(months(i))]];
    else
        monthNumStr = [monthNumStr; num2str(months(i))];
    end
end
monthStr = ["january","february","march","april","may","june","july","august",...
    "september","october","november","december"];

days1 = (1:28)'; days1str = [];
days2 = (1:29)'; days2str = [];
days3 = (1:30)'; days3str = [];
days4 = (1:31)'; days4str = [];

for i = 1:length(days1)
    if days1(i) < 10
        days1str = [days1str; ['0',num2str(days1(i))]];
    else
        days1str = [days1str; num2str(days1(i))];
    end
end
for i = 1:length(days2)
    if days2(i) < 10
        days2str = [days2str; ['0',num2str(days2(i))]];
    else
        days2str = [days2str; num2str(days2(i))];
    end
end
for i = 1:length(days3)
    if days3(i) < 10
        days3str = [days3str; ['0',num2str(days3(i))]];
    else
        days3str = [days3str; num2str(days3(i))];
    end
end
for i = 1:length(days4)
    if days4(i) < 10
        days4str = [days4str; ['0',num2str(days4(i))]];
    else
        days4str = [days4str; num2str(days4(i))];
    end
end

%%
monthNumStr = string(monthNumStr);
days1str = string(days1str);
days2str = string(days2str);
days3str = string(days3str);
days4str = string(days4str);

monthSearchList1 = [1 3 5 7 8 10 12];
monthSearchList2 = [4 6 9 11];
monthSearchList3 = 2;

for i = 1:length(regions)
    newMasterDir = strcat(masterDirectory,"/",regions(i));
    mainIndex(newMasterDir,yyyyStr,'index.html',strcat("region:",regions(i)));
    
    for j = 1:length(years)
        if years(j) == currentYear
            lia = ismember(months,(1:currentMonth)');
            monthNumStrNow = monthNumStr(lia);
            monthsNow = months(lia);
            m1I = monthSearchList1 <= currentMonth;
            monthSearchList1Now = monthSearchList1(m1I);
            m2I = monthSearchList2 <= currentMonth;
            monthSearchList2Now = monthSearchList2(m2I);
            m3I = monthSearchList3 <= currentMonth;
            monthSearchList3Now = monthSearchList3(m3I);
        else
            monthsNow = months;
            monthNumStrNow = monthNumStr;
            monthSearchList1Now = monthSearchList1;
            monthSearchList2Now = monthSearchList2;
            monthSearchList3Now = monthSearchList3;
        end
        newMasterDir = strcat(masterDirectory,"/",regions(i),"/",yyyyStr(j));
        mainIndex(newMasterDir,monthNumStrNow,'index.html',strcat("year:",yyyyStr(j)));
        
        for k = 1:length(monthsNow)
            if ismember(k,monthSearchList1Now)
                if k == currentMonth && years(j) == currentYear
                    days4strNow = days4str(days4 <= currentDay);
                else
                    days4strNow = days4str;
                end
                newMasterDir = strcat(masterDirectory,"/",regions(i),"/",yyyyStr(j),"/",monthNumStrNow(k));
                mainIndex(newMasterDir,days4strNow,'index.html',strcat("month:",monthStr(k)));
            elseif ismember(k,monthSearchList2Now)
                if k == currentMonth && years(j) == currentYear
                    days3strNow = days3str(days3 <= currentDay);
                else
                    days3strNow = days3str;
                end
                newMasterDir = strcat(masterDirectory,"/",regions(i),"/",yyyyStr(j),"/",monthNumStrNow(k));
                mainIndex(newMasterDir,days3strNow,'index.html',strcat("month:",monthStr(k)));
            elseif ismember(k,monthSearchList3Now) % month is february
                leapYearCheck = years(j);
                if leapyear(leapYearCheck)
                    if k == currentMonth && years(j) == currentYear
                        days2strNow = days2str(days2 <= currentDay);
                    else
                        days2strNow = days2str;
                    end
                    newMasterDir = strcat(masterDirectory,"/",regions(i),"/",yyyyStr(j),"/",monthNumStrNow(k));
                    mainIndex(newMasterDir,days2strNow,'index.html',strcat("month:",monthStr(k)));
                else
                    if k == currentMonth && years(j) == currentYear
                        days1strNow = days1str(days1 <= currentDay);
                    else
                        days1strNow = days1str;
                    end
                    newMasterDir = strcat(masterDirectory,"/",regions(i),"/",yyyyStr(j),"/",monthNumStrNow(k));
                    mainIndex(newMasterDir,days1strNow,'index.html',strcat("month:",monthStr(k)));
                end
            end
        end
    end
    
end

end

function mainIndex(masterDirectory,pages,outFile,masterName)
if nargin < 3; outFile = 'index.html'; end
if nargin < 4; masterName = 'ecuadorian tectonic and volcanic activity'; end

%%
masterDirectory = char(strcat("~/public_html/",masterDirectory));
if ~exist(masterDirectory,'dir')
    mkdir(masterDirectory);
end

cd(masterDirectory);
lPages = length(pages);

indexString = "<!DOCTYPE html>";
indexString(2) = "<html>";
indexString(3) = "<meta charset='utf-8'>";
indexString(4) = strcat("<head> <title>",masterName,"</title></head>");
indexString(5) = " ";
indexString(6) = "<style type=""text/css"">";
indexString(7) = "body";
indexString(8) = "{";
indexString(9) = "		margin: 0 auto;";
indexString(10) = "		width: 960px;";
indexString(11) = "		background: none repeat scroll 0 0 #fff;";
indexString(12) = "		padding: 0;";
indexString(13) = "		font-family: ""HelveticaNeue-Light"", ""Helvetica Neue Light"", ""Helvetica Neue"", Helvetica, sans-serif;";
indexString(14) = "  		font-weight: 300;";
indexString(15) = "	}";
indexString(16) = "	";
indexString(17) = "	header";
indexString(18) = "	{";
indexString(19) = "		background: #182F51;";
indexString(20) = "		float: center;";
indexString(21) = " }	";
indexString(22) = "	";
indexString(23) = "	";
indexString(24) = "	h1";
indexString(25) = "	{";
indexString(26) = "		color: #fff;";
indexString(27) = "		float: center;";
indexString(28) = "		font-size: 38px;";
indexString(29) = "	}";
indexString(30) = "	";
indexString(31) = "	footer";
indexString(32) = "	{";
indexString(33) = "		float: left;";
indexString(34) = "		width: 940px;";
indexString(35) = "		display: block;";
indexString(36) = "		font-size: 13px;";
indexString(37) = "		background: #cbcbcb;";
indexString(38) = "		padding: 5px 10px;";
indexString(39) = "		margin-top:10px;";
indexString(40) = "	}";
indexString(41) = "	";
indexString(42) = "	.author";
indexString(43) = "	{";
indexString(44) = "		float: right;";
indexString(45) = "		width: 50%;";
indexString(46) = "		text-align: right;";
indexString(47) = "	}";
indexString(48) = "	";
indexString(49) = "	/* Style the tab */";
indexString(50) = "	.tab ";
indexString(51) = "	{";
indexString(52) = "  	overflow: hidden;";
indexString(53) = "  	border: 1px solid #ccc;";
indexString(54) = "  	background-color: #f1f1f1;";
indexString(55) = "	}";
indexString(56) = "	";
indexString(57) = "	/* Style the buttons inside the tab */";
indexString(58) = "	.tab button ";
indexString(59) = "	{";
indexString(60) = "  	background-color: inherit;";
indexString(61) = "  	float: center;";
indexString(62) = "  	border: none;";
indexString(63) = "  	outline: none;";
indexString(64) = "  	cursor: pointer;";
indexString(65) = "  	padding: 14px 14px;";
indexString(66) = "  	transition: 0.3s;";
indexString(67) = "  	font-size: 36px;";
indexString(68) = "	}";
indexString(69) = "	";
indexString(70) = "	/* Change background color of buttons on hover */";
indexString(71) = " 	.tab button:hover ";
indexString(72) = " 	{";
indexString(73) = "  	background-color: #ddd;";
indexString(74) = " 	}";
indexString(75) = "</style>";
indexString(76) = " ";
indexString(77) = "<body>";
indexString(78) = "<header>";
indexString(79) = "	<h1 align=""center"">";
indexString(80) = strcat("			",masterName);
indexString(81) = "	</h1>";
indexString(82) = "</header>";
indexString(83) = " ";
indexString(84) = "<p style=""font-size:18px;"" align=""justify"">";
indexString(85) = "this is a generic repository for a variety of images and tables related to ecuadorian";
indexString(86) = "tectonic and volcanic activity. ideally, this repository will serve as a staging ";
indexString(87) = "ground for new ideas that could eventually/potentially be ";
indexString(88) = "incorporated into other internal systems (e.g., 'sam'), but at a future date.";
indexString(89) = "</p>";
indexString(90) = " ";
indexString(91) = "<p align=""center"">";
indexString(92) = "<b>";
indexString(93) = "*use of this site is for the exclusive use of scientists at the instituto geof&iacutesico de la ";
indexString(94) = "escuela polit&eacutecnica nacional and/or their invited colleagues.";
indexString(95) = "</b>";
indexString(96) = "</p>";
indexString(97) = " ";
indexString(98) = "<div class=""tab"" align=""center"">";

n = length(indexString);

for i = 1:lPages
    n = n+1;
    indexString(n) = strcat("<button onclick=""window.location.href = '",pages(i),"';"">",pages(i),"</button>");
    createDir = char(strcat(masterDirectory,"/",pages(i)));

    if ~exist(createDir,'dir')
        mkdir(createDir);
    end
end
indexString(n+1) = "</div>";
indexString(n+2) = "";
indexString(n+3) = "<footer>";
indexString(n+4) = "	<div class='author'>";
indexString(n+5) = "		by: stephen hernandez";
indexString(n+6) = "	</div>";
indexString(n+7) = "</footer>";
indexString(n+8) = " ";
indexString(n+9) = "</body>";
indexString(n+10) = "</html>";

PWD = pwd;
disp(PWD);
fileID = fopen(outFile,'w');
fprintf(fileID,'%s\n',indexString);
fclose(fileID);
end
