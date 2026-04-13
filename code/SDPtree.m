function out = SDPtree(D,maxDepth,mergeType,fmeFolder)

jobID = randi(100000);
algParam = struct;
algParam.maxDepth = maxDepth;
algParam.heurRuleCh = mergeType;
algParam.unrooted = 1;
algParam.maxDepthType = 'log';
outSDP = solveSDP(D,algParam);
[F_SDP_SPR,T_SDP_SPR,nRearSPR] = runFME(fmeFolder,jobID,D,outSDP.T);
out.F = F_SDP_SPR;
out.T = T_SDP_SPR;