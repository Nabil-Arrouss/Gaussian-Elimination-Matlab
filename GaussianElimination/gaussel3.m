% Gaussian Elimination for computing inverse of an square matrix

function [invM, dete] = gaussel3(MatrixA)

% Determine the length of the square matrix.
s = length(MatrixA);

% Initialize the identity matrix with the same dimensions as MatrixA.
I = eye(s);

% deteTemp used to initialize the determinant to 1 and then multiplied by the diagonal elements of the matrix.
deteTemp = 1;

% Ask the user if they want to display matrices during computation.
user = input("Enter '1' to display matrices or '0' to skip displaying:");

% Performing Gaussian elimination for the lower triangle part of the matrix.
for c = 1:s
    for r = c+1:s
        %  Check for zero pivot (MatrixA(c, c)) and throw an error if the matrix
        % can't be solved without row swap.
        if MatrixA(c, c) == 0
            error("Matrix cannot be solved without row swap.");
        end

        % Calculate the coefficient for row operations.
        coef = MatrixA(r, c) / MatrixA(c, c);

        % Apply the row operations to both the input matrix and the identity matrix.
        % This maintains the equivalence between the operations performed on
        % MatrixA and the transformations required to compute its inverse.
        MatrixA(r, :) = MatrixA(r, :) - coef * MatrixA(c, :);
        I(r, :) = I(r, :) - coef * I(c, :);

        % Display matrices if the user wants.
        if user
            display([MatrixA, I]);
        end
    end
end

% Computing the determinant by multiplying the diagonal elements.
% The determinant of a triangular matrix is the product of its diagonal elements.
for c = 1:s
    deteTemp = deteTemp * MatrixA(c, c);
end

% Performing Gaussian elimination for the upper triangle part of the matrix.
for r = 1:s
    for c = r+1:s
        % Check for zero pivot and throw an error if the matrix
        % can't be solved without row swap.
        if MatrixA(c, c) == 0
            error("Matrix cannot be solved without row swap.");
        end

        % Calculate the coefficient for row operations.
        coef = MatrixA(r, c) / MatrixA(c, c);

        % Apply the row operations to both the input matrix and the identity matrix.
        MatrixA(r, :) = MatrixA(r, :) - coef * MatrixA(c, :);
        I(r, :) = I(r, :) - coef * I(c, :);

        % Display matrices if the user wants.
        if user
            display([MatrixA, I]);
        end
    end
end

% Turn the diagonal to ones to get the identity matrix.
% This step normalizes the matrix, completing the process of converting MatrixA into an identity matrix
% and I into the inverse of the original MatrixA.
for c = 1:s
    I(c, :) = I(c, :) / MatrixA(c, c);
    MatrixA(c, :) = MatrixA(c, :) / MatrixA(c, c);

    % Display matrices if the user wants.
    if user
        display([MatrixA, I]);
    end
end

% The inverse matrix is stored in 'invM'.
invM = I;
% The determinant is stored in 'dete'.
dete = deteTemp;
end

% TEST:

%  test 1: Regular 3x3 Matrix 
% MatrixA1 = [4, 7, 2; 3, 5, 1; 2, 4, 3];
% [invM1, dete1] = gaussel3(MatrixA1);
% disp('Inverse of Matrix A1:');
% disp(invM1);
% disp('Determinant of Matrix A1:');
% disp(dete1);


% test 2: Identity Matrix
% MatrixA4 = eye(3);
% [invM4, dete4] = gaussel3(MatrixA4);
% disp('Inverse of Matrix A4 (Should be Identity Matrix):');
% disp(invM4);
% disp('Determinant of Matrix A4 (Should be 1):');
% disp(dete4);

% test 3: Diagonal Matrix
% MatrixA5 = diag([3, 5, 7]);
% [invM5, dete5] = gaussel3(MatrixA5);
% disp('Inverse of Matrix A5:');
% disp(invM5);
% disp('Determinant of Matrix A5:');
% disp(dete5);

