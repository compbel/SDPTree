function [obj_FME,T_fme,nRear] = runFME(fmeFolder,jobID,D,tree0)
    
    auxFolder = [fmeFolder filesep 'auxFiles'];
    if ~isfolder(auxFolder)
        mkdir(auxFolder);
    end

    n = length(D);
    if n <= 3
        obj_FME = NaN;
        return;
    end
    leafList = "T"  + (1:n);
    outfileP = [auxFolder filesep 'test' int2str(jobID) '.phy'];
    outfileN = [auxFolder filesep 'test' int2str(jobID) '.nwk'];
    infile_fme = [auxFolder filesep 'test' int2str(jobID) '_res.nwk'];
    statfile_fme = [auxFolder filesep 'test' int2str(jobID) '_stat.txt'];
    if exist(outfileP,"file")
        delete(outfileP);
    end
    if exist(outfileN,"file")
        delete(outfileN);
    end
    if exist(infile_fme,"file")
        delete(infile_fme);
    end
    if exist(statfile_fme,"file")
        delete(statfile_fme);
    end
    writePhylip(D, outfileP);


    if ~isempty(tree0)
        writeNewick(tree0,outfileN);
        command = [ fmeFolder '/fastme -i ' auxFolder '/test' int2str(jobID) ...
            '.phy -u ' auxFolder '/test' int2str(jobID) '.nwk ... ' ...
        '-o ' auxFolder '/test' int2str(jobID) '_res.nwk -I ' auxFolder '/test' int2str(jobID) '_stat.txt -s'];
    else
        command = [fmeFolder '/fastme -i ' auxFolder '/test' int2str(jobID) '.phy -m UNJ ... ' ...
        '-o ' auxFolder '/test' int2str(jobID) '_res.nwk -I ' auxFolder '/test' int2str(jobID) '_stat.txt -s'];
    end

    [status, out] = system(command); 
    out

    phyloFME = phytreeread(infile_fme);
    T_fme = phytree2digraph(phyloFME,leafList);
    T_fme.Edges.Weight = [];
    S = string(T_fme.Nodes.Name);
    T_fme.Nodes.Name = regexprep(S, '(?i)^branch\s*(\d+)$', 'B$1');
    S = string(T_fme.Nodes.Name);
    T_fme.Nodes.Name = regexprep(S, '^T(\d+)$', '$1'); 
    Tunr_fme = rooted2unrooted(T_fme);

    leaves_unr = (degree(Tunr_fme)==1);
    D_FME = distances(Tunr_fme);
    D_FME = D_FME(leaves_unr,leaves_unr);
    obj_FME = sum(sum(D./(2.^(D_FME))));

    txt = fileread(statfile_fme);
    nRear = str2double(regexp(txt, 'Performed\s+(\d+)\s+SPR\(s\)\.', 'tokens', 'once'));

    if nRear > 0
        pat = sprintf('SPR\\s+%d:\\s+new tree length is\\s+([0-9]*\\.?[0-9]+)', nRear);
        obj_FME = str2double(regexp(txt, pat, 'tokens', 'once'));
    end

    if exist(outfileP,"file")
        delete(outfileP);
    end
    if exist(outfileN,"file")
        delete(outfileN);
    end
    if exist(infile_fme,"file")
        delete(infile_fme);
    end
    if exist(statfile_fme,"file")
        delete(statfile_fme);
    end


