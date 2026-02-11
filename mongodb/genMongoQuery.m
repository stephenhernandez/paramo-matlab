function mongoQuery = genMongoQuery(field,operation,value,fmtStr)
mongoQuery = sprintf('{""%s"": {""$%s"": %s}}',field,operation,num2str(value,fmtStr));