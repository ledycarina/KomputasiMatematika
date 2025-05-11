% Penyelesaian contoh 4.1 dengan metode Jacobi
% Definisi sistem persamaan
A = [10  -1   2   0;
     -1  11  -1   3;
      2  -1  10  -1;
      0   3  -1   8];
  
b = [6; 25; -11; 15];

% Cek dominasi diagonal
diag_dom = all(abs(diag(A)) > sum(abs(A),2) - abs(diag(A)));
if ~diag_dom
    error('Matriks tidak dominan diagonal - metode mungkin tidak konvergen');
end

% Penyelesaian Awal
X0 = zeros(size(b));

% Panggil fungsi Jacobi
solusi = jacobi(A, b, X0);

%Tampilkan hasil
fprintf('\n=== Solusi Akhir ===\n');
fprintf('x1 = %.4f\nx2 = %.4f\nx3 = %.4f\nx4 = %.4f\n', ...
        solusi(1), solusi(2), solusi(3), solusi(4));


% % Contoh 4.1 penyelesaian dengan metode Gauss Seidel
% A = [10  -1   2   0;
%      -1  11  -1   3;
%       2  -1  10  -1;
%       0   3  -1   8];
% 
% b = [6; 25; -11; 15];
% 
% % Tebakan awal
% X0 = zeros(4,1);
% 
% % Parameter
% toleransi = 1e-6;
% iter_maks = 100;
% 
% % Panggil fungsi
% [solusi, iter] = gauss_seidel(A, b, X0, toleransi, iter_maks);
% 
% % Tampilkan hasil akhir
% fprintf('\n=== HASIL AKHIR ===\n');
% fprintf('Iterasi berhenti pada k = %d\n', iter);
% fprintf('Solusi sistem:\n');
% fprintf('x1 = %.4f\n', solusi(1));
% fprintf('x2 = %.4f\n', solusi(2));
% fprintf('x3 = %.4f\n', solusi(3));
% fprintf('x4 = %.4f\n', solusi(4));

% function Contoh10_1_Riemann_Trapesium 
% % Menggunakan Aturan Riemann (kiri/tengah/kanan) dan Trapesium
% 
% f = @(x) 2*x.^3;  % Fungsi integran
% a = 0; b = 1;     % Batas integrasi
% h = 0.1;          % Lebar langkah
% n = (b-a)/h;      % Jumlah subinterval
% 
% % Hitung dengan berbagai metode
% Riemann_kiri = riemann(f, [a b], n, 'kiri');
% Riemann_tengah = riemann(f, [a b], n, 'tengah');
% Riemann_kanan = riemann(f, [a b], n, 'kanan');
% Trapesium = trapesium(f, [a b], n);
% 
% % Tampilkan hasil
% disp(['Hasil Riemann Kiri   : ', num2str(Riemann_kiri)]);
% disp(['Hasil Riemann Tengah : ', num2str(Riemann_tengah)]);
% disp(['Hasil Riemann Kanan  : ', num2str(Riemann_kanan)]);
% disp(['Hasil Trapesium      : ', num2str(Trapesium)]);
% disp(['Nilai Eksak         : 0.5']);  % ?2x? dx = x?/2 |_0^1 = 0.5
% 
% % Fungsi Aturan Riemann
% function In = riemann(f, x, n, metode)
%     h = (x(2)-x(1))/n;
%     xvek = x(1):h:x(2)-h;  % Titik ujung subinterval
%     
%     switch metode
%         case 'kiri'
%             x_titik = xvek;               % Titik ujung kiri
%         case 'tengah'
%             x_titik = xvek + h/2;         % Titik tengah
%         case 'kanan'
%             x_titik = xvek + h;           % Titik ujung kanan
%     end
%     
%     yvek = f(x_titik);
%     In = h * sum(yvek);
% end
% 
% % Fungsi Aturan Trapesium (versi dikoreksi)
% function Tn = trapesium(f, x, n)
%     h = (x(2)-x(1))/n;
%     xvek = x(1):h:x(2);
%     yvek = f(xvek);
%     Tn = h/2 * (yvek(1) + 2*sum(yvek(2:n)) + yvek(n+1));
% end
% 
% end