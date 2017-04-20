function result = DyNB_caller(file,nCPU,walltime)
    tic;
    % read in data
    [pathstr,name,ext] = fileparts(file);
    data = load(file);
    param = regexp(name,'\d*M_(?<repl>\d*)rep_(?<TP>\d*)TP(_[0-9\.]*wn)?_block_\d*_formated','names');
    TP=str2num(param.TP);
    data.t = [1:TP]'/TP;
    data.t_star = [0:TP]'/TP;
    data.t_indices = arrayfun(@(x) find(x == data.t_star),data.t);
    cluster = parcluster('BrutusLSF24h');
    %cluster.SubmitArguments = ['-W ',walltime,' -R "select[model==Opteron8380]"'];
    cluster.SubmitArguments = ['-W ',walltime];
    pool = parpool(cluster,str2num(nCPU));
    addAttachedFiles(pool,'~/dynb');

    % check if time points are included in the t_star vector
    if ~issame(data.t,data.t_star(data.t_indices))
        disp('ERROR: t should be included in t_star!')
        return
    end

    % get the coefficients from the dispersion estimation
    coefficients = estimate_dispersion(data.contr_i,data.treat_i,data.conSize,data.treSize);
    result = cell(numel(data.treat_i),1);

    parfor idx=1:length(data.treat_i)
        % let us check that the dimensions of the input variables make sense
        if size(data.contr_i{idx},2) ~= length(data.t) || size(data.treat_i{idx},2) ~= length(data.t)
            disp('ERROR: size(data1,2) and size(data2,2) should be equal to length(t)')
            continue
        end
        if ~issame(size(data.contr_i{idx}),size(data.conSize)) || ~issame(size(data.treat_i{idx}),size(data.treSize))
            disp('ERROR: size(data1) should be equal to size(data1_sizefactors) and size(data2) should be equal to size(data2_sizefactors)!')
            continue
        end

        if (sum(data.contr_i{idx}(:) == 0) == numel(data.contr_i{idx}))|| (sum(data.treat_i{idx}(:) == 0) == numel(data.treat_i{idx}))
            disp('WARNING: not enough reads!')
            continue
        end

        result{idx} = DyNB(double(data.contr_i{idx}),double(data.treat_i{idx}),data.t,data.t_star,data.t_indices,data.conSize,data.treSize,coefficients);
        result{idx}.label = data.labels_i{idx};
    end
    delete(pool);
    time=datestr(toc/24/3600, 'DD:HH:MM:SS');
    outFile=strrep(name,'_formated','_dynb_result.mat');
    outDir=strrep(pathstr,'data','result');
    [s,mess,messid] = mkdir(outDir);
    save(strcat(outDir,'/',outFile),'result','time');
    return
end

function coefficients = estimate_dispersion(data1,data2,data1_sizefactors,data2_sizefactors)
    % estimate dispersion
    x = [];
    y = [];
    for idx=1:length(data1)
            % NOTICE: the timepoint t=0 is shared here
            data = [reshape(double(data1{idx})./data1_sizefactors,[numel(data1{idx}) 1]); reshape(double(data2{idx}(:,2:end))./data2_sizefactors(:,2:end),[numel(data2{idx}(:,2:end)) 1])];
            q = mean(data)';
        % shouldn't this be data1./data1_sizefactor?
            w = mean(var([double(data1{idx})./data1_sizefactors double(data2{idx}(:,2:end))./data2_sizefactors(:,2:end)]))';
            x = [x; q];
            y = [y; w];
    end
    z = x*.1/length(data).*sum(sum(1./[data1_sizefactors data2_sizefactors(:,2:end)]));
    tmp = y;
    X = [ones(size(x(x > 0 & tmp > 0 & ~isinf(tmp)))) log10(x(x > 0 & tmp > 0 & ~isinf(tmp))) log10(x(x > 0 & tmp > 0 & ~isinf(tmp))).^2];
    Y = log10(tmp(x > 0 & tmp > 0 & ~isinf(tmp)));
    coefficients = (X'*X)\X'*Y;
end