function [X, iter] = gauss_seidel_diktat(A, b, X, tol, max_iter)
% Untuk menyelesaikan Contoh 4.1
% Format pemanggilan:
% [X, iter] = gauss_seidel_diktat(A, b, X0, tol, max_iter)
% 
% Output:
% X - Solusi sistem persamaan
% iter - Jumlah iterasi yang dilakukan

if nargin < 3
    X = zeros(size(b)); % Nilai awal default
end
if nargin < 4
    tol = 1e-6; % Toleransi default
end
if nargin < 5
    max_iter = 1000; % Iterasi maksimum default
end

n = size(A, 1);
X_old = X; % Simpan nilai sebelumnya
err = tol + 1; % Inisialisasi error
iter = 0;

fprintf('\nMemulai iterasi Gauss-Seidel untuk Contoh 4.1\n');
fprintf('Iter\tx1\t\tx2\t\tx3\t\tx4\t\tError\n');

% Implementasi persis seperti diktat
while (iter <= max_iter) && (err > tol)
    % Persamaan pertama (i=1)
    X(1) = (b(1) - A(1,2:n)*X(2:n)) / A(1,1);
    
    % Persamaan tengah (i=2 sampai n-1)
    for i = 2:n-1
        X(i) = (b(i) - A(i,1:i-1)*X(1:i-1) - A(i,i+1:n)*X_old(i+1:n)) / A(i,i);
    end
    
    % Persamaan terakhir (i=n)
    X(n) = (b(n) - A(n,1:n-1)*X(1:n-1)) / A(n,n);
    
    % Hitung error
    err = norm(X - X_old);
    
    % Tampilkan progres iterasi
    fprintf('%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
            iter, X(1), X(2), X(3), X(4), err);
    
    X_old = X;
    iter = iter + 1;
end

if iter > max_iter && err > tol
    warning('Iterasi maksimum tercapai sebelum konvergen');
else
    fprintf('\nIterasi berhenti pada iterasi ke-%d\n', iter-1);
    fprintf('Dengan error akhir: %.8f\n', err);
end
end