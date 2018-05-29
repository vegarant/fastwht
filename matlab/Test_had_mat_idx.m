classdef Test_had_mat_idx < matlab.unittest.TestCase
    % Test_had_mat_idx tests the had_mat_idx function
    %
    % USAGE:
    % >> testCase = Test_had_mat_idx;
    % >> res = run(testCase)
    %
    properties
        eps;
    end  
    
    methods (Test)
        % Constructor
        function obj = Test_had_mat_idx(testCase)
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
            ei(k) = 1;
            ei = N*fwht(ei);
            x = had_mat_idx(N,1:N, k);

            testCase.verifyTrue(norm(ei-x,'fro') < eps);
        end
        
        function testVectorN128kSingleArg(testCase)

            eps = testCase.eps;

            N = 128;
            k = randi(N);
            ei = zeros(N,1);
            ei(k) = 1;
            ei = N*fwht(ei)';
            x = had_mat_idx(N,k, 1:N);
            
            testCase.verifyTrue(norm(ei-x,'fro') < eps);
        end

        function testNargGtrTKarg(testCase) 

            eps = testCase.eps;
            N = 256;
            s = 5;
            k = randi(N,[1,s]);
            X = had_mat_idx(N, 1:N, k);

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
            X = had_mat_idx(N, n, 1:N);
            
            W = N*fwht(eye(N));
            W = W(n,:);
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
                    X(i,j) = had_mat_idx(N,i,j);
                end
            end
            testCase.verifyTrue(norm(W-X,'fro') < eps);
        end
    end
end 

