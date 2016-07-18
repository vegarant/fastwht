% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright 2016 Vegard Antun
%

eps = 1e-10;

% Test that both transforms create the same result
for N = [32,64,128]

    U1 = zeros(N,N);
    U2 = zeros(N,N);

    for order = {'sequency', 'dyadic', 'hadamard'}

        for i = 1:N

            x = zeros(N,1);
            x(i) = 1;

            y1 = fwht(x,N, order{1});
            y2 = fastwht(x, N, order{1});

            U1(:,i) = y1;
            U2(:,i) = y2;

        end

        zeroNorm = norm(U1-U2);
        assert(abs(zeroNorm) < eps, 'fastwht not equal to fwht 1');
    end
end

% Test with shorter length array
% NOTE: Matlab fwht does not allow the use of N value which is not a power of 2,
%       while it does allow that the length of 'x' is not a power of 2.
%       In fastwht this inconsistency has been removed.

for N = [8, 15, 32];

    for order = {'sequency', 'dyadic', 'hadamard'}

        x = randn(N,1);

        y1 = fwht(x, 16, order{1});
        y2 = fastwht(x, 16, order{1});

        zeroNorm = norm(y1-y2);
        assert(abs(zeroNorm) < eps, 'fastwht not equal to fwht 2');

    end

end

% Test with row vector
for N = [32]

    for order = {'sequency', 'dyadic', 'hadamard'}

        x = randn(1,N);

        y1 = fwht(x, [], order{1});
        y2 = fastwht(x, [], order{1});

        [s1_M, s1_N] = size(y1);
        [s2_M, s2_N] = size(y1);

        assert((s1_M == s2_M & s1_N == s2_N), 'fastwht not equal to fwht 4');

    end

end

% Test with no order
for N = [32]

    x = randn(1,N);

    y1 = fwht(x, N);
    y2 = fastwht(x, N);

    zeroNorm = norm(y1-y2);
    assert(abs(zeroNorm) < eps, 'fastwht not equal to fwht 5');

end

% Test with only x as input argument
for N = [64]

    x = randn(1,N);

    y1 = fwht(x);
    y2 = fastwht(x);

    zeroNorm = norm(y1-y2);
    assert(abs(zeroNorm) < eps, 'fastwht not equal to fwht 5');

end

% Test with complex numbers
N = 64;

x_re = randn(N,1);
x_im = randn(N,1);

x = x_re + x_im*1j;

y1 = fwht(x);
y2 = fastwht(x);

zeroNorm = norm(y1-y2);
assert(abs(zeroNorm) < eps, 'fastwht not equal to fwht complex 6');



% Test two dimensional tensor product expansion of fastwht
% Create hadamard matrices

for M = [16,16]
    for N = [8, 32]
        for order = {'sequency', 'dyadic', 'hadamard'}
            % Create matrix of size M
            U_M = zeros(M,M);
             
            for i=1:M
                ei = zeros(M,1); ei(i) = 1;
                U_M(:,i) = fwht(ei, [], order{1});
            end 
             
            % Create matrix of size N
           U_N = zeros(N,N);
             
            for i=1:N
                ei = zeros(N,1); ei(i) = 1;
                U_N(:,i) = fwht(ei, [], order{1});
            end 

            X  = randn(M,N);
            H1 = U_M*X*U_N';
            H2 = fastwht(X, [M,N], order{1});
            
            zeroNorm = norm(H1-H2, 'fro');
            assert(abs(zeroNorm) < eps, 'fastwht tensor product not working');
        end
    end
end

