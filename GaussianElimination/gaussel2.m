%Computing Gaussian elimination with steps of partial and whole pivoting method.

function [x, L, U] = gaussel2(A, b, useWholePivoting)
    % gaussel2 is for complete pivoting and partial pivoting.
    
    [m, n] = size(A);  % Get the dimensions of the coefficient matrix A.
    x = zeros(n, 1);   % Initialize the solution vector x.
    L = eye(n);        % Initialize the lower triangular matrix L. (identity-matrix)
    U = A;             % Initialize the upper triangular matrix U.

    % Check if the user wants to use whole pivoting.
    if useWholePivoting
        % Use whole pivoting method.
        fprintf('Using whole pivoting method.\n');
        [x, L, U] = wp(A, b);
        return;
    
    else
        fprintf('Using partial pivoting method.\n');
        % Use partial pivoting method.
        for k = 1:n - 1
            % Performing Partial-pivoting
            for p = k + 1:n
                % Switching the biggest pivot element in the rows
                [~, maxRow] = max(abs(U(k:n, k)));
                maxRow = maxRow + k - 1;
            end

            for i = k + 1:n
                l = U(i, k) / U(k, k);  % Compute the multiplier l(i, k).
                L(i, k) = l;           % Store the multiplier in the lower triangular matrix.
                for j = k + 1:n
                    U(i, j) = U(i, j) - l * U(k, j);  % Update elements in the upper triangular matrix.
                end
                b(i) = b(i) - l * b(k);  % Update the right-hand side vector b.
            end

            % Display the matrix A after each step
            display(A);
        end
            % Performing Gaussian elimination
        for k = 1:n - 1
            for i = k + 1:n
                U(i, k) = 0;  % Zero out elements below the diagonal.
            end
            % Display the matrix A after zeroing out elements
            display(A);
        end

        % Solve the system using back-substitution
        x(n) = b(n) / U(n, n);
        for i = n - 1:-1:1
            sum = 0;
            for j = i + 1:n
                sum = sum + U(i, j) * x(j);
            end
            x(i) = (b(i) - sum) / U(i, i);
            display(A);
        end    

                % Check for singularity
                % Checks if a zero pivot is encountered, which would indicate a singular matrix. 
                % If so, and if whole pivoting is allowed, switches to the wp function
                if U(maxRow, k) == 0
                      if useWholePivoting
                         fprintf('Switching to whole pivoting method due to a zero pivot element.\n');
                         [x, L, U] = wp(A, b);
                      end     
                      return;
                end
    end        
end

% nested function for gaussian elimination with whole pivoting
function [x, L, U] = wp(A, b)
    [n, ~] = size(A);  % Get the number of rows (assumed to be a square matrix).
    p = 1:n;           % Initialize the pivot row permutation vector.
    q = 1:n;           % Initialize the pivot column permutation vector;

    for k = 1:n - 1
        % Find the maximum element in the current submatrix
        submatrix = A(k:n, k:n);
        [~, maxElementIdx] = max(abs(submatrix(:)));
        [maxRowIdx, maxColIdx] = ind2sub(size(submatrix), maxElementIdx);
        maxRow = maxRowIdx + k - 1;
        maxCol = maxColIdx + k - 1;

        % Check for singularity to ensure that system has a unique solution
        if A(maxRow, maxCol) == 0

            fprintf('Switching to whole pivoting method due to a zero pivot element.\n');
            [x, L, U] = wp(A, b);  % Switch to Gaussian elimination with whole pivoting
                      
            return;
        end

        % Swap rows in A
        A([k, maxRow], k:n) = A([maxRow, k], k:n);

        % Swap columns in A
        A(k:n, [k, maxCol]) = A(k:n, [maxCol, k]);

        % Update pivot permutation vectors
        p([k, maxRow]) = p([maxRow, k]);
        q([k, maxCol]) = q([maxCol, k]);

        if A(k, k) == 0
            break;
        end

        % Perform Gaussian elimination with whole pivoting
        A(k + 1:n, k) = A(k + 1:n, k) / A(k, k);
        i = k + 1:n;
        A(i, i) = A(i, i) - A(i, k) * A(k, i);

        % Display the matrix A after each step
        display(A);
    end

    L = tril(A, -1) + eye(n);  % Extract the lower triangular matrix L.
    U = triu(A);              % Extract the upper triangular matrix U.

    % Solve the system using back-substitution
    y = L \ b;
    x = U \ y;
end


% TEST
%  Solvable System with Partial Pivoting
%A1 = [3, -2, 5; 2, 6, -4; 1, -3, -3];
% b1 = [21; -38; 19];
% [x1, L1, U1] = gaussel2(A1, b1, false);
% disp('Solution using Partial Pivoting:');
% disp(x1);

% Solvable System with Whole Pivoting
% A2 = [3, -2, 5; 2, 6, -4; 1, -3, -3];
% b2 = [21; -38; 19];
% [x2, L2, U2] = gaussel2(A2, b2, true);
% disp('Solution using Whole Pivoting:');
% disp(x2);

