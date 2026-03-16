clear
addpath '/opt/gurobi1002/linux64/matlab'
addpath(genpath('/home/skumsp/YALMIP'))          
addpath('/home/skumsp/mosek/10.2/toolbox/r2017a','-end');

fmeFolder = './fastme';

ns = 10:5:50;
m = 1000;
p = 0.05;
nTests = 20;

% simulType = 'points';
% simulType = 'tree_seq';
% simulType = 'RIM';
simulType = 'RDSM';
% simulType = 'phylobench';


if strcmp(simulType,'RIM') || strcmp(simulType,'RDSM')
    nTests = 10;
end

algParam = struct;
algParam.unrooted = 1;

resFold = ['resSim_' simulType];
if ~isfolder(resFold)
    mkdir(resFold);
end

crashDir = 'crash_dumps';
if ~isfolder(crashDir)
    mkdir(crashDir);
end

for n = ns
    algParam.maxDepth = ceil(n/2);
    if strcmp(simulType,'phylobench')
        fileType = fullfile(['phylobench' filesep 'nexus-' int2str(n)], '*.nex');
        files = dir(fileType);
        nTests = length(files);
    end
    K = n-1;
    for i = 1:nTests
        try
            resFile = [resFold filesep 'res_' int2str(n) '_' int2str(i) '.mat'];
            if isfile(resFile)
                continue;
            end
            if strcmp(simulType,'tree_seq')
                [tree_true,~] = simulateYuleTree(n); 
                leafSeq = simulateEvolTree(tree_true,m,p,0);
                D = squareform(pdist(leafSeq,'hamming'));
            end
            if strcmp(simulType,'points')
                D = squareform(pdist(rand(n,m)));
            end
            if strcmp(simulType,'RIM')
                symbols = [num2cell('a':'z')];
                testID = ['RIM' int2str(n) symbols{i}];
                matrFold = ['randMatr' filesep testID];
                matrFile = [matrFold filesep testID '.txt'];
                D = readmatrix(matrFile);
                D = D(:,2:end);
            end
            if strcmp(simulType,'RDSM')
                symbols = [num2cell('a':'z')];
                testID = ['RDSM' int2str(n) symbols{i}];
                matrFold = ['randMatr' filesep testID];
                matrFile = [matrFold filesep testID '.txt'];
                D = readmatrix(matrFile);
                D = D(:,2:end);
            end
            if strcmp(simulType,'phylobench')
                filename = [files(i).folder filesep files(i).name];
                try
                    aln = multialignread(filename);
                catch
                    aln = readNexusAlignment(filename);
                end
                if length(aln)~=n
                    aln = readNexusAlignment(filename);
                end
                D   = seqpdist(aln, 'Method','Jukes-Cantor', 'Alphabet','AA','SquareForm',true);
            end
            
            [T0,~] = simulateCoalesTree(n);         
            [F_FME,T_FME,nRear_FME] = runFME(fmeFolder,i,D,T0);
            [F_FME_NJ,T_FME_NJ,nRear_FME_NJ] = runFME(fmeFolder,i,D,[]);
            [F_NJ, T_NJ] = runNJ(D,algParam);
        
            algParamCurr = algParam;
            algParamCurr.heurRuleCh = 'sep';
            t = tic;
            out = solveSDP(D,algParamCurr);
            tm = toc(t);
            F_SDP = out.F;
            T_SDP = out.T;
            [F_SDP_SPR,T_SDP_SPR,nRearSPR] = runFME(fmeFolder,i,D,T_SDP);
    
            algParamCurr.heurRuleCh = 'prof';
            t = tic;
            out1 = solveSDP(D,algParamCurr);
            tm1 = toc(t);
            F_SDP1 = out1.F;
            T_SDP1 = out1.T;
            [F_SDP_SPR1,T_SDP_SPR1,nRearSPR1] = runFME(fmeFolder,i,D,T_SDP1);
                    
            res = struct();
            res.D = D;
            res.T0 = T0;
            res.F_FME = F_FME;
            res.F_NJ = F_NJ;
            res.F_FME_NJ = F_FME_NJ;
            res.F_SDP = F_SDP;
            res.T_SDP = T_SDP;
            res.F_SDP_SPR = F_SDP_SPR;
            res.T_SDP_SPR = T_SDP_SPR;
            res.tm = tm;
            res.nRearSPR = nRearSPR;
            res.F_SDP1 = F_SDP1;
            res.T_SDP1 = T_SDP1;
            res.F_SDP_SPR1 = F_SDP_SPR1;
            res.T_SDP_SPR1 = T_SDP_SPR1;
            res.tm1 = tm1;
            res.nRearSPR1 = nRearSPR1;
            res.nRear_FME = nRear_FME;
            res.nRear_FME_NJ = nRear_FME_NJ;
            save(resFile,"-fromstruct",res);
        catch 
            ts = getCurrentTask(); 
            wid = 0; 
            if ~isempty(ts) 
                wid = ts.ID; 
            end
            fn = fullfile(crashDir, sprintf('crash_w%d_%s.mat', wid, datestr(now,'yyyymmdd_HHMMSSFFF')));
            es = struct();
            es.n = n;
            es.i = i;
            save(fn,'-fromstruct',es);              
            continue;
        end
    end
end