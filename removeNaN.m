
function [dates, ClPr] = removeNaN(dates, ClPr)

N = length(ClPr(1, :));

% Replace NaNs until common start
index = zeros(N, 1);
for i = 1:N
    index(i) = min(find(~isnan(ClPr(:, i))));
end
startIndex = max(index);
ClPr(1:startIndex-1, :) = [];
dates(1:startIndex-1, :) = [];

% Replacing next NaNs with previous values
[row, col] = find(isnan(ClPr));
for i = 1:length(col)
    ClPr(row(i), col(i)) = ClPr(row(i)-1, col(i));
end

end
