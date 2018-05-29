classdef TestWal < matlab.unittest.TestCase
    % TestWal tests the wal function
    %
    % USAGE:
    % >> testCase = TestWal;
    % >> res = run(testCase)
    %
    properties
        eps;
    end  
    
    methods (Test)
        % Constructor
        function obj = TestWal(testCase)
            obj.eps = 1e-8;
        end

        function test_w256(testCase)
            eps = testCase.eps;

            N = 256;
            W = N*fwht(eye(N));
            X = wal(0:N-1, (0:N-1)/N);
            testCase.verifyTrue(norm(W-X,'fro') < eps);

        end
        
        function test_one_value(testCase)
            eps = testCase.eps;

            N = 256;
            W = N*fwht(eye(N));
            for i = 1:30
                x = wal(130, i/N);
                y = W(131, i+1);
                testCase.verifyTrue(abs(x-y) < eps);
            end
        end
        
        function test_n_row_value(testCase)
            eps = testCase.eps;

            N = 512;
            W = N*fwht(eye(N));
            for i = 1:30
                x = wal(0:N-1, i/N);
                y = W(1:N, i+1);
                testCase.verifyTrue(norm(x-y) < eps);
            end
        end
        
      function test_n_col_value(testCase)
            eps = testCase.eps;

            N = 512;
            W = N*fwht(eye(N));
            for i = 1:30
                x = wal(i, (0:N-1)/N);
                y = W(i+1, 1:N);
                testCase.verifyTrue(norm(x-y) < eps);
            end
        end
        
        function testCloseTo1(testCase) 
            eps = testCase.eps;
            N = 256;
            x = linspace(0,1,N);
            x(end) = x(end) - 1e-14;

            W = N*fwht(eye(N));
            X = wal(0:N-1, x);

            testCase.verifyTrue(norm(W-X,'fro') < eps);
        end
        

    end
    methods (Access=private)
    end
end 
