%Tugas4

% % Figure 1
% t = 0:0.01:12*pi;
% x = sin(t).*(exp(cos(t))-2*cos(4*t)-(sin(t/12)).^5);
% y = cos(t).*(exp(cos(t))-2*cos(4*t)-(sin(t/12)).^5);
% figure(1)
% plot(x,y,'r','LineWidth',1.5)
% title('Figure Butterfly Curve')
% axis equal

% %Figure 2
% [x,y,z] = meshgrid(-2:0.1:2);
% f = (x.^2 + (9/4)*y.^2 + z.^2 - 1).^3 - x.^2.*z.^3 - (9/80)*y.^2.*z.^3;
% figure(2)
% p = patch(isosurface(x,y,z,f,0));
% set(p,'FaceColor','red','EdgeColor','none');
% view(3); axis equal; camlight; lighting gouraud
% title('3D Heart')

% %Figure 3
% x = linspace(-5,5,100);
% [X,Y] = meshgrid(x);
% psi = exp(-(X.^2+Y.^2)/2).*(X.^2-Y.^2);
% figure(3);
% surf(X,Y,abs(psi).^2); shading interp; title('Quantum Probability Density');

% %Figure 4
% u = linspace(0, 2*pi, 50);
% v = linspace(0, 2*pi, 50);
% [U,V] = meshgrid(u,v);
% X = (2+cos(U/2).*sin(V)-sin(U/2).*sin(2*V)).*cos(U);
% Y = (2+cos(U/2).*sin(V)-sin(U/2).*sin(2*V)).*sin(U);
% Z = sin(U/2).*sin(V) + cos(U/2).*sin(2*V);
% figure(4);
% surf(X,Y,Z); title('Klein Bottle (4D Surface Projection)');

% % Figure 5
% n = 200;
% x = linspace(-2.5, 1.5, n);
% y = linspace(-2, 2, n);
% [X,Y] = meshgrid(x,y);
% C = X + 1i*Y;
% Z = zeros(size(C));
% M = zeros(size(C));
% for k = 1:100
%     Z = Z.^2 + C;
%     M(abs(Z) < 2) = k;
% end
% figure(5);
% imagesc(x,y,M); colormap(jet); axis equal; title('Figure 5');

% % Figure 6
% [x,y] = meshgrid(-10:0.5:10);
% r = sqrt(x.^2 + y.^2);
% z = sin(r)./r;
% figure(6);
% surf(x,y,z); title('Sinc Function Surface');

% % Figure 7
% [x,y,z] = meshgrid(-2:0.4:2);
% r = sqrt(x.^2 + y.^2 + z.^2);
% Bx = x./r.^3; By = y./r.^3; Bz = z./r.^3;
% figure(7);
% quiver3(x,y,z,Bx,By,Bz);

% % Figure 8
% [x,y] = meshgrid(-1:0.1:1);
% z = x.^2 - y.^2;
% figure(8);
% surf(x,y,z); 

% % Figure 9
% c = -0.7 + 0.27i;
% n = 200;
% x = linspace(-1.5, 1.5, n);
% y = linspace(-1.5, 1.5, n);
% [X,Y] = meshgrid(x,y);
% Z = X + 1i*Y;
% M = zeros(size(Z));
% for k = 1:50
%     Z = Z.^2 + c;
%     M(abs(Z) < 2) = k;
% end
% figure(9);
% imagesc(x,y,M); colormap(hot); axis equal; 

% % Figure 10
% u = linspace(0, 2*pi, 50);
% v = linspace(-0.4, 0.4, 20);
% [U,V] = meshgrid(u,v);
% X = (1 + V.*cos(U/2)).*cos(U);
% Y = (1 + V.*cos(U/2)).*sin(U);
% Z = V.*sin(U/2);
% figure(10);
% surf(X,Y,Z);

% % Figure 11
% [x,y] = meshgrid(-2:0.1:2);
% z = -1./sqrt(x.^2 + y.^2 + 0.1);
% figure(11);
% surf(x,y,z); 

% % Figure 12
% u = linspace(0, 2*pi, 50);
% v = linspace(-1, 1, 20);
% [U,V] = meshgrid(u,v);
% X = cosh(V).*cos(U);
% Y = cosh(V).*sin(U);
% Z = V;
% figure(12);
% surf(X,Y,Z); 

% % Figure 13
% a = 0.2; b = 0.2; c = 5.7;
% dt = 0.01; T = 100; t = 0:dt:T;
% x = zeros(length(t),1); y = x; z = x;
% x(1) = 1; y(1) = 1; z(1) = 1;
% for i=1:length(t)-1
%     x(i+1) = x(i) + (-y(i)-z(i))*dt;
%     y(i+1) = y(i) + (x(i)+a*y(i))*dt;
%     z(i+1) = z(i) + (b+z(i)*(x(i)-c))*dt;
% end
% figure(13);
% plot3(x,y,z); 


% % Figure 14
% theta = linspace(0,pi,30);
% phi = linspace(0,2*pi,40);
% [Theta,Phi] = meshgrid(theta,phi);
% R = 1;
% X = R.*sin(Theta).*cos(Phi);
% Y = R.*sin(Theta).*sin(Phi);
% Z = R.*cos(Theta);
% W = R.*cos(Theta/2); % 4th dimension projection
% figure(14);
% scatter3(X(:),Y(:),Z(:),50,W(:),'filled'); 

% % Figure 15
% [x,y] = meshgrid(-3:0.1:3);
% r = sqrt(x.^2 + y.^2);
% theta = atan2(y,x);
% psi = exp(-r.^2/2).*exp(1i*2*theta);
% figure(22);
% surf(x,y,abs(psi).^2); shading interp; 


% % Figure 16
% [x,y,z] = meshgrid(-2:0.4:2);
% I = 1; mu0 = 1;
% r = sqrt(x.^2 + y.^2);
% Bx = -mu0*I*y./(2*pi*r.^2);
% By = mu0*I*x./(2*pi*r.^2);
% Bz = zeros(size(x));
% figure(16);
% quiver3(x,y,z,Bx,By,Bz); 


% % Figure 17
% x = linspace(-2,2,50);
% y = linspace(-2,2,50);
% [X,Y] = meshgrid(x,y);
% Z = X + 1i*Y;
% W = log(Z);
% figure(17);
% surf(X,Y,imag(W)); 

% % Figure 18
% x = linspace(-10,10,500);
% t = linspace(0,5,100);
% [X,T] = meshgrid(x,t);
% c = 2; % wave speed
% u = 0.5*c*sech(sqrt(c)/2*(X-c*T)).^2;
% figure(18);
% surf(X,T,u); shading interp; 

% % Figure 19
% r = linspace(2,20,50);
% phi = linspace(0,2*pi,50);
% [R,Phi] = meshgrid(r,phi);
% X = R.*cos(Phi);
% Y = R.*sin(Phi);
% Z = 1./sqrt(R-2); % potential
% figure(19);
% surf(X,Y,Z); 

% % Figure 20
% [x,y] = meshgrid(-5:0.1:5);
% V = 100*(sqrt(x.^2 + y.^2) > 3 & sqrt(x.^2 + y.^2) < 3.2);
% E = 5;
% psi = besselj(0,sqrt(E)*sqrt(x.^2 + y.^2)).*exp(-0.1*(x.^2 + y.^2));
% figure(20);
% surf(x,y,abs(psi).^2 + V); shading interp; 

% % Figure 21
% theta = linspace(0,pi/2,50);
% phi = linspace(0,2*pi,50);
% [Theta,Phi] = meshgrid(theta,phi);
% X = sin(Theta).*cos(Phi);
% Y = sin(Theta).*sin(Phi);
% Z = cos(Theta);
% W = sin(2*Theta); % entanglement measure
% figure(21);
% scatter3(X(:),Y(:),Z(:),50,W(:),'filled'); 

% % Figure 22
% theta = linspace(0,2*pi,50);
% phi = linspace(0,2*pi,50);
% [theta,phi] = meshgrid(theta,phi);
% r = 4*(1-cos(theta)/2);
% x = (r.*cos(theta)+6).*cos(phi);
% y = (r.*cos(theta)+6).*sin(phi);
% z = r.*sin(theta);
% figure(22)
% surf(x,y,z)
% axis equal

% % Figure 23
% sigma = 10; beta = 8/3; rho = 28;
% f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
% [t,a] = ode45(f,[0 100],[1 1 1]);
% figure(23)
% plot3(a(:,1),a(:,2),a(:,3))
% grid on

% % Figure 24
% phi = linspace(-pi,pi,50);
% theta = linspace(-pi/2,pi/2,50);
% [phi,theta] = meshgrid(phi,theta);
% a = 1; b = 1; m = 6; n1 = 1; n2 = 7; n3 = 8;
% r = (abs(cos(m*phi/4)/a).^n2 + abs(sin(m*phi/4)/b).^n3).^(-1/n1);
% x = r.*cos(phi).*cos(theta);
% y = r.*sin(phi).*cos(theta);
% z = r.*sin(theta);
% figure(24)
% surf(x,y,z)
% axis equal

% % Figure 25
% [X,Y] = meshgrid(-10:0.5:10);
% R = sqrt(X.^2 + Y.^2) + eps;
% Z = besselj(0,R);
% figure(25)
% surf(X,Y,Z)

% % Figure 26
% theta = 0:0.01:2*pi;
% k = 5;
% r = cos(k*theta);
% figure(26)
% polar(theta, r, 'r')

% % Figure 27
% [x,y,z] = meshgrid(-5:0.5:5);
% psi = exp(-(x.^2+y.^2+z.^2)/4).*(x+y.*1i);
% figure(27)
% slice(x,y,z,abs(psi).^2,[0],[0],[-3 0 3])
% shading interp

% % Figure 28
% [x,y] = meshgrid(0:0.01:1);
% noise = randn(size(x));
% filter = fspecial('gaussian',[50 50],5);
% cmb = imfilter(noise,filter);
% figure(28)
% imagesc(cmb)
% colormap jet

% % Figure 29
% theta = linspace(0,pi,50);
% phi = linspace(0,2*pi,50);
% [Theta,Phi] = meshgrid(theta,phi);
% Psi = cos(Theta/2).*exp(1i*0) + sin(Theta/2).*exp(1i*Phi);
% figure(29)
% surf(cos(Phi).*sin(Theta), sin(Phi).*sin(Theta), cos(Theta), angle(Psi))
% shading interp

% % Figure 30
% [x,y] = meshgrid(-10:0.2:10);
% r = sqrt(x.^2 + y.^2);
% z = exp(-r.^2/10).*sin(3*r-20*0);
% figure(30)
% surf(x,y,z)

% % Figure 31
% figure(31)
% [u,v] = meshgrid(linspace(0,2*pi,200), linspace(0,pi,200));
% a = 1; b = 0.5; c = 0.3; n = 5;
% % Persamaan parametrik kompleks
% x = (a + b*cos(n*u).*sin(n*v)).*cos(u).*sin(v);
% y = (a + b*cos(n*u).*sin(n*v)).*sin(u).*sin(v);
% z = (c + b*sin(n*u).*sin(n*v)).*cos(v);
% surf(x,y,z,'EdgeColor','none','FaceAlpha',0.9)
% title('FontSize',16,'FontWeight','bold')
% colormap(hsv)
% axis equal off
% light('Position',[1 1 1],'Style','infinite')
% lighting gouraud
% material([0.4 0.6 0.5 30 1.0])
% view(45,30)

% % Figure 31
% figure(32)
% 
% % Hitung nilai r sesuai rumus
% theta_deg = 40; % sudut dalam derajat
% theta_rad = deg2rad(theta_deg); % konversi ke radian
% r = sin(theta_rad)^2 + cos(theta_rad);
% 
% % Hitung x dan y
% theta = 0:0.01:2*pi; % sudut dari 0 sampai 2pi
% x = r * cos(theta);
% y = r * sin(theta);
% 
% % Plot lingkaran dengan radius r
% plot(x, y, 'b', 'LineWidth', 2);
% hold on;
% 
% % Plot titik khusus pada theta = 0
% plot(r*cos(0), r*sin(0), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 
% % Plot garis radius
% plot([0 r*cos(0)], [0 r*sin(0)], 'r--');
% 
% % Pemetaan dengan f(z) = e^2
% e_squared = exp(2); % Nilai e^2
% plot(e_squared*cos(theta), e_squared*sin(theta), 'g--', 'LineWidth', 1.5);
% 
% % Anotasi
% text(0, r+0.2, sprintf('r = sin?(40?) + cos(40?) = %.3f', r), ...
%     'HorizontalAlignment', 'center');
% text(e_squared, 0, sprintf('f(z) = e? ? %.3f', e_squared), ...
%     'VerticalAlignment', 'bottom');
% 
% % Format plot
% title('Visualisasi Rumus E-MATH', 'FontSize', 14);
% xlabel('x'); ylabel('y');
% axis equal;
% grid on;
% legend('Lingkaran r', 'Titik (r,0)', 'Garis radius', 'Lingkaran e^2', ...
%     'Location', 'best');

% % Figure 31
% figure(31)
% clf;
% % Parameter
% a = 0.12;
% zoom_factor = 0.13;
% % Buat grid untuk x, y, z
% [x, y, z] = meshgrid(linspace(-1, 1, 50), linspace(-1, 1, 50), linspace(-1, 1, 50));
% % Hitung komponen pertama rumus
% term1 = 0.5*z;
% component1 = 64*term1.^7 - 112*term1.^5 + 56*term1.^3 - 7*term1;
% part1 = a*0.99*component1;
% % Hitung komponen kedua rumus (polynomial produk)
% part2 = (0.7818314825 - 0.3765101982*y - 0.7818314825*x) ...
%       .* (0.7818314824 - 0.8460107361*y - 0.1930964297*x) ...
%       .* (0.7818314825 - 0.6784479340*y + 0.5410441731*x) ...
%       .* (0.7818314825 + 0.8677674789*x) ...
%       .* (0.7818314824 + 0.6784479330*y + 0.541044172*x) ...
%       .* (0.7818314824 + 0.8460107358*y - 0.193096429*x) ...
%       .* (0.7818314821 + 0.3765101990*y - 0.781831483*x) ...
%       .* z.^2;
% % Gabungkan kedua komponen
% F = part1 - part2;
% % Visualisasi dengan isosurface
% isovalue = 0; % Nilai isosurface yang akan divisualisasikan
% p = patch(isosurface(x, y, z, F, isovalue));
% isonormals(x, y, z, F, p);
% % Set properties visual
% set(p, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
% colormap(jet);
% light('Position', [1 1 1], 'Style', 'infinite');
% lighting gouraud;
% material shiny;
% axis([-zoom_factor zoom_factor -zoom_factor zoom_factor -zoom_factor zoom_factor]);
% view(3);
% title('Visualisasi Rumus Kompleks 3D', 'FontSize', 14);
% xlabel('X'); ylabel('Y'); zlabel('Z');
% colorbar;
% grid on;
% rotate3d on;