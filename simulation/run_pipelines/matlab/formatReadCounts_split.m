function [] = formatReadCounts_filtered(folder,split)
    fnames = dir(strcat(folder,'/*.mat'));
    fnames = {fnames.name}';
    FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'once'));
    files = fnames(FIND('.*(TP|wn).mat'));
    for idx=1:numel(files)
        param = regexp(files{idx},'\d*M_(?<repl>\d*)rep_(?<TP>\d*)TP(_[0-9\.]*wn)?.mat','names');
        file = strcat(folder,'/',files{idx});
        load(file);
        treat = cell(length(treatment),1);
        contr = cell(length(control),1);
        for index=1:length(control) contr(index)=mat2cell(control(:,:,index),str2num(param.repl),str2num(param.TP));end
        for index=1:length(treatment) treat(index)=mat2cell(double(treatment(:,:,index)),str2num(param.repl),str2num(param.TP));end

        % split data sets into 'split' blocks in order not to waste resources on brutus as the parfor loop seems to quite CPU inefficient ( only 8% used .. )
        block_size = round(linspace(0,size(treat,1),str2num(split)+1));
        for blockIdx=1:str2num(split)
            contr_i = contr(block_size(blockIdx)+1:block_size(blockIdx+1));
            treat_i = treat(block_size(blockIdx)+1:block_size(blockIdx+1));
            labels_i = labels(block_size(blockIdx)+1:block_size(blockIdx+1));
            
            [pathstr,name,ext] = fileparts(file);
            outDir=strcat(folder,'/',name);
            [s,mess,messid] = mkdir(outDir);
            outFile=strcat(outDir,'/',strrep(files{idx},'.mat',strcat('_block_',num2str(blockIdx),'_formated.mat')));
            save(outFile,'contr_i','treat_i','labels_i','conSize','treSize');
        end
    end
    %delete(strcat(folder,'/',files.name));
    return;
end
