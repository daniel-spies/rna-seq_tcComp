function result = extractResultsMatlab(file)
    inFile = load(file);
    outFile = strrep(file,'.mat','.txt');
    
    if ischar(inFile.result{1})
        printCell(outFile,inFile.result);
    else
        candidates = getCandidates(inFile.result);
        printCell(outFile,candidates);
    end
end

function printCell(fileName,data)
    fileID = fopen(fileName,'w');
    [nrows,ncols] = size(data);
    for row = 1:nrows
        fprintf(fileID,'%s\t%5.4e\n',data{row,:});
    end    
    fclose(fileID);
end

function candidates = getCandidates(result,labels)
    candidates = {};
    for idx =1:size(result,1)
       if ~isempty(result{idx})
           candidates{end+1,1} = result{idx}.label;
           candidates{end,2} = result{idx}.bayes_factor;
       end
    end
end