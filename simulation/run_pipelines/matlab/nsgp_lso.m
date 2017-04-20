function result = nsgp_hmc(file,nCPU,walltime)
    tic;
    [pathstr,name,ext] = fileparts(file);
    data = load(file);

    result = cell(length(data.contr_i),1);
    method='grad';
    
    param = regexp(name,'\d*M_(?<repl>\d*)rep_(?<TP>\d*)TP(_[0-9\.]*wn)?_block_\d*_formated','names');
    TP = linspace(0,str2num(param.TP),str2num(param.TP))';
    replicates=str2num(param.repl);
    labels = repelem(TP,replicates);

    cluster = parcluster('BrutusLSF24h');
    %cluster.SubmitArguments = '-W 8:00 -R "select[model==Opteron8380]"';
    cluster.SubmitArguments = ['-W ',walltime];
    pool = parpool(cluster,str2num(nCPU));
    addAttachedFiles(pool,'~/adaptivegp/code');

    parfor i = 1:length(data.contr_i)
        % make sure data is double format as otherwise later on it will throw an error in the logmvnpdf function
        contr = double(data.contr_i{i}(:));
        treat = double(data.treat_i{i}(:));
        MLL =[];

        if (sum(contr) == 0)
            contr = rand(length(contr),1);
        elseif (length(unique(contr))<2)
            contr = contr.*((0.8 + rand(length(contr),1)*(1.2-0.8)));
        end

        if (sum(treat) == 0)
            treat = rand(length(treat),1);
        elseif (length(unique(treat))<2)
            treat = treat.*((0.8 + rand(length(treat),1)*(1.2-0.8)));
        end

        [gpmodel,samples,mll,mse,nmse,nlpd] = nsgp(labels,contr, 'lso', method);
        MLL(2)=max(mll);
        [gpmodel,samples,mll,mse,nmse,nlpd] = nsgp(labels,treat, 'lso', method);
        MLL(3)=max(mll);
        [gpmodel,samples,mll,mse,nmse,nlpd] = nsgp([labels;labels],[contr;treat], 'lso', method);
        MLL(1)=max(mll);

        out = struct('bayes_factor',prod(MLL(2:3))/MLL(1),'label',data.labels_i{i});
        result{i} = out;
    end
    delete(pool);
    time=datestr(toc/24/3600, 'DD:HH:MM:SS');
    outFile=strrep(name,'_formated','_nsgp_lso_result.mat');
    outDir=strrep(pathstr,'data','result');
    [s,mess,messid] = mkdir(outDir);
    save(strcat(outDir,'/',outFile),'result','time');
    return
end