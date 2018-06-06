% Computes the Walsh function
% 
% `wal` is the sequency ordered Walsh function on [0,1).
% 
% INPUT:
% n - Walsh frequency (number of sign changes on [0,1)) n = 0,1,2,... 
% k - Real numbers in the interval [0,1).
%
% OUTPUT:
% The values of the Walsh function. If n and k are vectors the function returns
% the matrix 
%           [ wal(n(1), k(1)),   ..., wal(n(1), k(end));
%             ...            ,   ..., ...              ; 
%             wal(n(end), k(1)), ..., wal(n(end), k(end))]             
% 
function y = wal(n,k);
    
    if (exist('wal') ~= 3)
        error('fastwht: wal function is not compiled');
    end 
    
end
