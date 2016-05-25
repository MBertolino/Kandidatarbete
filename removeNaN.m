function [dates, clPr] = removeNaN(dates, clPr)
%{
 removeNaN removes all rows until the first row has no NaN's.
 Then all NaN's are replaced with the previous row's value.

 Input:
 dates (matrix) - dates for each asset
 clPr (matrix) - closing prices for each asset with NaN's

 Output:
 dates (matrix) - dates for each asset w/ out NaN-days
 clPr (matrix) - closing prices for each asset w/ out NaN's

 2016 Iliam Barkino, Mattias Bertolino
%}

N = length(clPr(1, :));

% Replace NaNs until common start
index = zeros(N, 1);
for i = 1:N
    index(i) = min(find(~isnan(clPr(:, i))));
end
startIndex = max(index);
clPr(1:startIndex-1, :) = [];
dates(1:startIndex-1, :) = [];

% Replacing next NaNs with previous values
[row, col] = find(isnan(clPr));
for i = 1:length(col)
    clPr(row(i), col(i)) = clPr(row(i)-1, col(i));
end

end
