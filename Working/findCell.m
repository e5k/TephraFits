% Finds the cell index of a string in a cell array
function idx = findCell(cell2find, str2find)
cell2find   = cellfun(@num2str, cell2find, 'UniformOutput', false);
idx         = find(not(cellfun('isempty', strfind(cell2find, str2find))),1);

% Check that the string to find isn't only a subset of a longer string
if ~strcmpi(cell2find(idx), str2find)
    idx = [];
end
