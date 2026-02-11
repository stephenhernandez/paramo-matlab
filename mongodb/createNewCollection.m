function collectionNames = createNewCollection(conn,dbName,collectionNames,verboseFlag)
lia = ismember(dbName,collectionNames);
sumlia = sum(lia);
if sumlia
    return;
end

if verboseFlag
    fprintf('%s as collection does not exist, creating....\n',dbName);
end

try
    createCollection(conn,dbName);
    collectionNames = conn.CollectionNames;
catch ME
    warning(ME);
    return;
end