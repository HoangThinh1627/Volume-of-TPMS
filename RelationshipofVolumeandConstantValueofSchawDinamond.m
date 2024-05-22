clc;
close all;
clear;
loop = 1;
% plot1 = zeros(3,1);
k = 1;
for c = 0.1:0.05:1
d = pi;   
d2 = pi; 
s = pi/40 ; 

[x, y, z] = meshgrid(-d:s:d2, -d:s:d2, -d:s:d2);

Zmax = max(z(:));
Zmin = min(z(:));

m1 = 5; % Max Number of Unit cells 
n1 = 1; % Min Number of Unit cells

m = m1 / 2;
n = n1 / 2;

k1 = (m - n) / (Zmax - Zmin);
c0 = (k1 * Zmin * Zmin) / 2;
c1 = -(Zmin * k1) + n;

% Uniform scaling function applied only to Z
g = (k1 * z + c1);

u = cos(x).*cos(y).*cos(z) - sin(x).*sin(y).*sin(z);


% Sheet
S = (u + c) .* (u - c);

% Compute isosurfaces and caps Primitive 
[F1, V1] = isosurface(x, y, z, S, 0);
[F2, V2] = isocaps(x, y, z, S, 0, 'below');

% Combine faces and vertices Primitive
F3 = [F1; F2 + size(V1, 1)];
V3 = [V1; V2];

% Create the patch object for the visualization
% P = patch('Vertices', V3, 'Faces', F3, 'FaceColor', 'red', 'EdgeColor', 'none');

% Surface Area Calculation
SA = 0;

for i = 1:size(F3, 1)
    %Primitive
    v1 = V3(F3(i, 1), :);
    v2 = V3(F3(i, 2), :);
    v3 = V3(F3(i, 3), :); 
    edge1 = v2 - v1;
    edge2 = v3 - v1;
    area = 0.5 * norm(cross(edge1, edge2));
    SA = SA + area;    
end

% Volume Calculation Primitive
VF = permute(reshape(V3(F3,:),[size(F3) 3]),[3 1 2]);
Vol = 1/6*sum(dot(cross(VF(:,:,1),VF(:,:,2),1),VF(:,:,3),1));

plotp(1,loop+1) = c;
plotp(2,loop+1) = Vol/((2*d)^3);
plotp(3,loop+1) = SA;
loop = loop + 1;
end
figure;
axis equal;
plot(plotp(2,:), plotp(1,:), '*', 'MarkerSize', 7,'Color', 'red');
axis([0 1 0 5])
ylabel('Constant value, c');
xlabel('Relative Density, \rho');
title('Relationship of Relative Density and Constant value');
legend('Primitive');
hold on
grid on
hold on
n = 3/0.1+1;
m = 5;
for i = 1:2*m
    xsum(i) = sum(plotp(2,:).^(i));
end
a(1,1) = n;
b(1,1) = sum(plotp(1,:));
for j = 2:m+1
    a(1,j) = xsum(j-1);
end
for i = 2:m+1
    for j = 1:m+1
        a(i,j) = xsum(j+i-2);
    end
    b(i,1) = sum(plotp(2,:).^(i-1).*plotp(1,:));
end
p = (a\b)';
for i = 1:m+1
    Pcoef(i) = p(m+2-i);
end
epsilon = 0:0.015:1.2;
ptest = [37.027 -67.980 8.187 57.591 -44.188 11.489 0.840 0];
cfit = polyval(Pcoef,epsilon);
cfit1 = polyval(Pcoef,plotp(2,:));
 plot(plotp(2,:), cfit1,'LineWidth', 2);
% plot(epsilon, voltest,'LineWidth', 2,'Color', 'blue');
% legend('Primitive','Equation for Data','t(p)=37.027p^7+-67.980p^6+8.187p^5+57.591p^4+-44.188p^3+11.489p^2+0.840p+0');
legend('Schwarz Dinamond','Equation for Data');
SStot = sum((plotp(1,:) - mean(cfit1)).^2); % Total sum of squares
SSres = sum((plotp(1,:) - cfit1).^2); % Residual sum of squares
R_squared = 1 - (SSres / SStot); % R-squared value
% Add R^2 annotation
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', ['R^2 = ', num2str(R_squared)], ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white');
ylabel('Constant value, c');
xlabel('Relative Density, \rho');
title('Relationship of Relative Density and Constant value');

