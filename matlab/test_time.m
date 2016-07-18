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

n_start = 10;
fwht_end = 18;
n_end   = 27;


t_fwht = [];
t_fastwht = [];
t_fft = [];

for k = n_start:fwht_end

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

for k = fwht_end+1:n_end

    N = 2^k;
    x = randn(N,1);
    
    f2 = @() fastwht(x);
    f3 = @() fft(x);

    t2 = timeit(f2);
    t3 = timeit(f3);

    t_fastwht = [t_fastwht, t2];
    t_fft = [t_fft, t3];

end

lwidth = 2;
fsize  = 15;

fig = figure();
hold('on');
plot(n_start:fwht_end, t_fwht, 'linewidth', lwidth );
plot(n_start:n_end, t_fastwht, 'linewidth', lwidth );
plot(n_start:n_end, t_fft, 'linewidth', lwidth );
legend({'fwht', 'fastwht', 'fft'}, 'fontsize', fsize, 'Location', 'northwest');
xlabel('Length of vector 2^{R}', 'fontsize', fsize);
ylabel('Time (seconds)', 'fontsize', fsize);
set(gca,'yscale','log');

saveas(fig, 'compare_performance.png');

