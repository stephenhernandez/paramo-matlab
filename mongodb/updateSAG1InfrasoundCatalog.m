function updateSAG1InfrasoundCatalog(S,SAG1InfrasoundArray,allMySNCLs,updatedSNCLs,...
    writeToDB,conn,dbName1,dbName2)

lia = ismember(SAG1InfrasoundArray,updatedSNCLs);
if ~sum(lia)
    fprintf("necessary SAG1 sncls havent been updated, no change, return.\n");
    return;
else
    fprintf("\n");
end

[lia,locb] = ismember(SAG1InfrasoundArray,allMySNCLs);
sumlia = sum(lia);
if ~sumlia
    return;
end

%%
S = S(locb(lia));
[vI,nUsed,t,amp,baz,vel,meanCC,medCC] = ...
    process_sag1_infrasound_array_v1(S);

lT = length(t);
if ~lT
    return;
end

if ~writeToDB
    return;
end

%%
epochRef = datetime(1970,01,01);
nDocs = count(conn,dbName1);
if ~nDocs
    tLast = max(t);
    tLastStr = tLast;
    t = milliseconds(t - epochRef);
    tLast = max(t);
    insert(conn,dbName2,table(tLast,tLastStr));
    insert(conn,dbName1,table(t,amp,baz,vel,meanCC,medCC,nUsed,vI));
    return;
end

%%
t2 = max(t);
t = milliseconds(t - epochRef);
tLast = find(conn,dbName2);
tLast = tLast.tLast;
tI = t > tLast;
if ~sum(tI)
    fprintf(1,"no new events, doing nothing\n");
    return;
end

%%
remove(conn,dbName2,"{}");
tLast = max(t);
tLastStr = t2;
disp(tLastStr);
insert(conn,dbName2,table(tLast,tLastStr));

%%
t = t(tI);
amp = amp(tI);
baz = baz(tI);
vel = vel(tI);
meanCC = meanCC(tI);
medCC = medCC(tI);
nUsed = nUsed(tI);
vI = vI(tI);
insert(conn,dbName1,table(t,amp,baz,vel,meanCC,medCC,nUsed,vI));
