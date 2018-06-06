% Computes indices in the sequency ordered Hadamard matrix
% 
% This function computes
% >> W = N*fwht(eye(N));
% >> W(n,k) % This is what the function computes
% 
% INPUT:
% N - Size of Hadamard matrix, N = 2^r for positive integer r.
% n - Row number.
% k - Column number.
%
% OUTPUT:
% Returns the subsampled rows and columns of the sequency ordered Hadamard
% matrix.
% 
function y = had_mat_idx(N, n,k);
    
    if (exist('had_mat_idx') ~= 3)
        error('fastwht: had_mat_idx function is not compiled');
    end 
    
end

