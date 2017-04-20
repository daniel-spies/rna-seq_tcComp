function result =combine_results_split(folder,outFile,nBlock)
    % extract bayes factor and gene from each block and merge into one cell
    filesIn = dir(strcat(folder,'/*.mat'));
    result = {};
    time = 0;
    for file = 1:length(filesIn)
       data = load(strcat(folder,'/',filesIn(file).name));
       data.result = data.result(~cellfun(@isempty,data.result)); %% filter out possible empty cells
       time = time + str2double(regexp(data.time,':','split'));
       ext = cell(size(data.result,1),2);
       for idx=1:size(data.result,1)
          ext{idx,1} = data.result{idx,1}.label;
          ext{idx,2} = data.result{idx,1}.bayes_factor;
       end
       result = [result;ext];
    end
    % convertig back to seconds to sum up and convert back to time format
    outTime = time(1) * 24 * 3600 + time(2) * 3600 + time(3) * 60 + time(4);
    outTime = datestr(outTime/24/3600/str2num(nBlock),'mm:DD:HH:MM:SS');

    % write list of DEG genes and their gene_label
    printCell(strcat(outFile,'.txt'),result);
    save(strcat(outFile,'.mat'),'result','outTime');
end

function printCell(fileName,data)
    fileID = fopen(fileName,'w');
    [nrows,ncols] = size(data);
    for row = 1:nrows
        fprintf(fileID,'%s\t%5.4e\n',data{row,:});
    end    
    fclose(fileID);
end