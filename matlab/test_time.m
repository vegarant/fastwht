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

%
% This program tests the performance of the different implementations.
%

n_start = 6;
n_end   = 18;

resolution = n_start:n_end;

t_fwht = [];
t_fastwht = [];
t_fft = [];

for k = resolution

    N = 2^k;
    x = randn(N,1);

    f1 = @() fwht(x);
    f2 = @() fastwht(x);
    f3 = @() fft(x);

    t1 = timeit(f1);
    t2 = timeit(f2);
    t3 = timeit(f3);

    t_fwht = [t_fwht, t1];
    t_fastwht = [t_fastwht, t2];
    t_fft = [t_fft, t3];

end

lwidth = 2;
fsize  = 15;

fig = figure();
hold('on');
plot(resolution, t_fwht(resolution-n_start+1), 'linewidth', lwidth );
plot(resolution, t_fastwht(resolution-n_start+1), 'linewidth', lwidth );
plot(resolution, t_fft(resolution-n_start+1), 'linewidth', lwidth );
legend({'fwht', 'fastwht', 'fft'}, 'fontsize', fsize, 'Location', 'northwest');
xlabel('Length of vector 2^{R}', 'fontsize', fsize);
set(gca,'yscale','log');

%saveas(fig, 'compare_performance.png');

