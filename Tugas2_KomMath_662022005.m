%whileloop : untuk mengulang ketika kondisinya benar / Perulangan while mengeksekusi pernyataan secara berulang selama kondisinya benar.

% a = 15;
% while( a < 35 )
%    fprintf('value of a: %d\n', a);
%    a = a + 1;
% end

% Nested if adalah fungsi yang sepenuhnya terdapat dalam fungsi induk.
% x = 350 ;
% y = 700; 
%    if( x == 350 )
%       if( y == 700 )
%          fprintf('Value of x is 350 and y is 700\n' );
%       end
%    end
%    fprintf('Exact value of x is : %d\n', x );
%    fprintf('Exact value of y is : %d\n', y );
   

%menentukan Kategori Usia
usia = input('Masukkan usia Anda: ');

if usia >= 0
    if usia < 13
        disp('Kategori: Anak-anak');
    elseif usia < 18
        disp('Kategori: Remaja');
    elseif usia < 60
        disp('Kategori: Dewasa');
    else
        disp('Kategori: Lansia');
    end
else
    disp('Usia tidak valid');
end
