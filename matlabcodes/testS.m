% Function: this is used for getting all the original series for test
% Input: n-the number of state
%Output: all the posible of parts series
function [ o ] = testS(n)
    b = [-14;-12;-11;-10;-9;-7;-3;-2;-1;1;2;3;4;5;6;7;8;9;10;11;12;13;14];
    o = b;
    for i = 2:n
        l = size(o,1);
        temp = repmat(b',l,1);
        temp = temp(:);
    end




end

