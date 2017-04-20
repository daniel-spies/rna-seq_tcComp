function [] = formatReadCounts(folder)
    d = dir(folder);
    isub = [d(:).isdir];
    nameFolds = {d(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..','raw'})) = []; %exclude folders

    for subfolder=1:numel(nameFolds)
        fnames = dir(strcat(folder,'/',nameFolds{subfolder},'/*.mat'));
        fnames = {fnames.name}';
        FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'once'));
        files = fnames(FIND('.*(TP|wn).mat'));

        for idx=1:numel(files)
            param = regexp(files{idx},'\d*M_(?<repl>\d*)rep_(?<TP>\d*)TP(_[0-9\.]*wn)?.mat','names');
            load(strcat(folder,'/',nameFolds{subfolder},'/',files{idx}));
            treat = cell(length(treatment),1);
            contr = cell(length(control),1);
            for index=1:length(control) contr(index)=mat2cell(control(:,:,index),str2num(param.repl),str2num(param.TP));end
            for index=1:length(treatment) treat(index)=mat2cell(double(treatment(:,:,index)),str2num(param.repl),str2num(param.TP));end
            outFile=strcat(folder,'/',nameFolds{subfolder},'/',strrep(files{idx},'.mat','_formated.mat'));
            save(outFile,'contr','treat','labels','conSize','treSize');
        end
    end
    return;
end
