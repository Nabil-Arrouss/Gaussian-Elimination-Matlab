    % Computing Gaussian elimination 
    % Function gaussel1 with two input parameters: coefficient matrix 'MatrixA' and right-side vector 'VectorB'.
    function [x] = gaussel1(MatrixA,VectorB)
    
    % right-side Vector Length 
    N = length(VectorB);
    
    % The solution vector 'X' initialized to 0, it will store solutions later
    X = zeros(N,1);
    
    % Concatenate the coefficient matrix 'MatrixA' and right-side vector 'VectorB' 
    %    to form an augmented matrix 'AugMatrix'
    AugMatrix = [MatrixA VectorB];
    
    % in case of underdetermined LES. Checks if the number of columns in
    % MatrixA is greater than the number of rows in VectorB, if true, it
    % indicates an underdetermined system
    if size(MatrixA, 2) > size(VectorB, 1)
        disp("This is an underdetermined LES");
    end
    
    % If the user want to see the steps
    user = input(" Enter a number: \n 1. Display the matrices during computation. \n 0. Do Not display the matrices \n:");
    
    % The opperation for GE
    for j = 1:N-1
        for i = j+1:N
            % Calculate the multiplier for row operations
            m = AugMatrix(i,j)/AugMatrix(j,j); 
            % Row operation to eliminate the element below the pivot 
            AugMatrix(i,:) = AugMatrix(i,:) - m*AugMatrix(j,:);
    
            % If the user want to display the steps of GE.
            % If row or coloumn swap needed, throw an error and terminate the program.
            if AugMatrix(i,j+1) == 0
                error(" GE can't be solved without row or column swap!") ;
            end
        end
            % Display the augmented matrix if the user chooses
            if user
                  disp(AugMatrix);
            end
    end
    
    % Assigning the solutions to the solution vector 'X'
    X(N) = AugMatrix(N,N+1)/AugMatrix(N,N);
    
    % Calculating what are the solutions are equal to. Back-substitution (technique to solve system of linear equations)
    for k = N-1:-1:1
        X(k) = (AugMatrix(k,N+1)- AugMatrix(k,k+1:N)*X(k+1:N))/AugMatrix(k,k);
            % Display the augmented matrix during back-substitution if the user chooses
            if user
                disp(AugMatrix);
            end
    end
    
    x=X;
    end


% TEST:
% simple 2x2 system
% A = [2, 3; 4, 1];
% B = [5; 6];
% x = gaussel1(A, B);
% disp(x);

% 3x3 system
% A = [3, -2, 5; 2, 6, -4; 1, -3, -3];
% B = [21; -38; 19];
% x = gaussel1(A, B);
% disp(x);

% Testing with a singular matrix (no solution)
% A = [2, 4, 6; 1, 2, 3; -1, -2, -3];
% B = [18; 9; -9];
% x = gaussel1(A, B);
% disp(x);

% underdetermined system 
% A = [1, 2; 3, 4; 5, 6];
% B = [7; 8; 9];
% x = gaussel1(A, B);
% disp(x);