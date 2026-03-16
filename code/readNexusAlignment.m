function aln = readNexusAlignment(filename)
% Read NEXUS alignment MATRIX block (supports interleave).
% Returns struct array aln(i).Header, aln(i).Sequence (like multialignread).

L = regexp(fileread(filename), '\r?\n', 'split');

L = regexprep(L, '\[[^\]]*\]', '');

i0 = find(strcmpi(strtrim(L), 'matrix'), 1, 'first');
if isempty(i0), error('No MATRIX block found.'); end

i1 = i0 + find(contains(strtrim(L(i0+1:end)), ';'), 1, 'first');
if isempty(i1), error('MATRIX block has no terminating ";".'); end

mp = containers.Map('KeyType','char','ValueType','char');
order = {}; 

for k = i0+1 : i1-1
    s = strtrim(L{k});
    if s=="" || startsWith(s,'['), continue; end
    s = regexprep(s, '\s+', ' ');      
    tok = regexp(s, '^(\S+)\s+(\S+)$', 'tokens', 'once'); 
    if isempty(tok), continue; end
    nm = tok{1}; chunk = tok{2};

    if isKey(mp,nm)
        mp(nm) = [mp(nm) chunk];      
    else
        mp(nm) = chunk;
        order{end+1} = nm; %#ok<AGROW>
    end
end

n = numel(order);
aln = repmat(struct('Header','', 'Sequence',''), n, 1);
for i = 1:n
    aln(i).Header   = order{i};
    aln(i).Sequence = mp(order{i});
end
end
