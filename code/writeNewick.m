function writeNewick(G, filename)
    
    root = find(indegree(G)==0, 1);
    n  = numnodes(G);
    L  = find(outdegree(G)==0);           
    leafID = zeros(n,1); 
    leafID(L) = 1:numel(L);

    function s = build(u)
        ch = successors(G,u);
        if isempty(ch)                 
            s = sprintf('T%d', leafID(u));
        else                              
            s = sprintf('(%s:1,%s:1)', build(ch(1)), build(ch(2)));
        end
    end

    nwk = [build(root) ';'];
    fid = fopen(filename,'w'); fprintf(fid,'%s\n', nwk); fclose(fid);
end