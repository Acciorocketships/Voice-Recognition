function output = innerprod(x,y)
if nargin < 2
    y = x;
end
if isrow(x)
    x = x';
else
    if ~iscolumn(x)
        error('The x input is not a row or column vector');
    end
end
if isrow(y)
    y = y';
else
    if ~iscolumn(y)
        error('The y input is not a row or column vector');
    end
end
if length(x) ~= length(y)
    error('The length of x and y must be equal');
end
output = sum(x .* conj(y));
return