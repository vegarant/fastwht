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
        
        function test128SingleArg(testCase)
            N = 128;
            testEquailityWithMatrixSingleArg(testCase, N);
        end

        function testVectorN128(testCase)

            eps = testCase.eps;

            N = 128;
            k = 50;
            ei = zeros(N,1);
            ei(k+1) = 1;
            ei = N*fwht(ei);
            x = wal(N,0:N-1, k);

            testCase.verifyTrue(norm(ei-x,'fro') < eps);
        end
        
        function testVectorN128kSingleArg(testCase)

            eps = testCase.eps;

            N = 128;
            k = randi(N);
            ei = zeros(N,1);
            ei(k+1) = 1;
            ei = N*fwht(ei)';
            x = wal(N,k, 0:N-1);
            
            testCase.verifyTrue(norm(ei-x,'fro') < eps);
        end

        function testNargGtrTKarg(testCase) 

            eps = testCase.eps;
            N = 256;
            s = 5;
            k = randi(N,[1,s]);
            X = wal(N, 0:N-1, k-1);

            W = zeros(N,s);
            for i = 1:s
                W(k(i),i) = 1;
            end
            W = N*fwht(W);
            testCase.verifyTrue(norm(W-X,'fro') < eps);

        end

        function testKargGtrTNarg(testCase) 

            eps = testCase.eps;
            N = 256;
            s = 5;
            n = randi(N,[1,s]);
            X = wal(N, n-1, 0:N-1);
            
            W = N*fwht(eye(N));
            W = W(n,:);
            testCase.verifyTrue(norm(W-X,'fro') < eps);
            
        end

        function testIn01_w256(testCase)
            eps = testCase.eps;
            
            N = 256;
            W = N*fwht(eye(N));
            X = wal(N,0:N-1, (0:N-1)/N);
            testCase.verifyTrue(norm(W-X,'fro') < eps);
            
        end
        
        function testCloseTo1(testCase) 
            eps = testCase.eps;
            N = 256;
            x = linspace(0,1,N);
            x(end) = x(end) - 1e-14;
            
            W = N*fwht(eye(N));
            X = wal(N,0:N-1, x);
            
            testCase.verifyTrue(norm(W-X,'fro') < eps);
        end
        

    end
    methods (Access=private)
        function testEquailityWithMatrixSingleArg(testCase, N)
            eps = testCase.eps;
            W = N*fwht(eye(N));
            X = zeros(N);
            for i = 1:N
                for j = 1:N
                    X(i,j) = wal(N,i-1,j-1);
                end
            end
            testCase.verifyTrue(norm(W-X,'fro') < eps);
        end
    end
end 
