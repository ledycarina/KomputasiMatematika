function X = jacobi(A, b, X0, max_iter, tol)
% Input:
%   A        - Matriks koefisien (n x n)
%   b        - Vektor hasil (n x 1)
%   X0       - Penyelesaian Awal
%   max_iter - Maksimum iterasi (default: 1000)
%   tol      - Toleransi error (default: 1e-5)
% Output:
%   X        - Solusi sistem

% Penanganan input default
if nargin < 5
    tol = 1e-5;
end
if nargin < 4
    max_iter = 1000;
end
if nargin < 3
    X0 = zeros(size(b));
end

n = length(b);
X = X0;

fprintf('\n=== Iterasi Jacobi ===\n');
fprintf('Iter\tx1\t\tx2\t\tx3\t\tx4\t\tError\n');

for k = 1:max_iter
    X_old = X;
    
    % Perhitungan tiap komponen
    X(1) = (b(1) - (A(1,2)*X_old(2) + A(1,3)*X_old(3))) / A(1,1);
    X(2) = (b(2) - (A(2,1)*X_old(1) + A(2,3)*X_old(3) + A(2,4)*X_old(4))) / A(2,2);
    X(3) = (b(3) - (A(3,1)*X_old(1) + A(3,2)*X_old(2) + A(3,4)*X_old(4))) / A(3,3);
    X(4) = (b(4) - (A(4,2)*X_old(2) + A(4,3)*X_old(3))) / A(4,4);
    
    % Hitung error
    err = norm(X - X_old, inf);
    
    % Tampilkan progres
    fprintf('%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
            k, X(1), X(2), X(3), X(4), err);
    
    % Cek konvergensi
    if err < tol
        fprintf('\nKonvergen pada iterasi ke-%d\n', k);
        break;
    end
end
end 