function dec_perms = GetRepeatedPermutation(options, length)
    dec_perms = [];
    if (length < 1)
        return;
    end
    dec_cell = cell(length, 1);
    [dec_cell{:}] = ndgrid(options);
    dec_perms = cellfun(@(options){options(:)}, dec_cell);
    dec_perms = [dec_perms{:}];
end