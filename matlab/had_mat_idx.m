% Computes indices in the sequency ordered Hadamard matrix
% 
% This function computes 
% >> W = N*fwht(eye(N));
% >> W(n,k) % This is what the function computes
% 
% Note that in the code snipped above, we compute the Walsh-Hadamard transform 
% N times (one time for each column in eye(N)). The implemented function only 
% computes the Walsh-Hadamart transform of columns of eye(N), needed to 
% evaluate W(n,k) 
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

