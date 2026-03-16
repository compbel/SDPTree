function writePhylip(D, filename)
% writePhylip Write an n×n matrix to a PHYLIP distance file.
% - Uses taxon names T1..Tn (padded/truncated to 10 chars).

    [n, ~] = size(D);
    names = arrayfun(@(k) sprintf('T%d', k), 1:n, 'UniformOutput', false);
    names10 = cellfun(@(s) pad(s(1:min(end,10)), 10, 'right'), names, 'UniformOutput', false);

    fid = fopen(filename, 'w');
    fprintf(fid, '%d\n', n);

    fmtNum = repmat('%10.6f', 1, n);
    for i = 1:n
        fprintf(fid, '%-10s', names10{i});
        fprintf(fid, fmtNum, D(i, :));
        fprintf(fid, '\n');
    end
    fclose(fid);
end