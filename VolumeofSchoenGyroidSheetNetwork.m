clc;
close all;
clear;
loop=1;
for c = 0.1:0.05:1.5
d = pi;   
d2 = pi; 
s = pi/30 ; 

[x, y, z] = meshgrid(-d:s:d2, -d:s:d2, -d:s:d2);

Zmax = max(z(:));
Zmin = min(z(:));

m1 = 1; % Max Number of Unit cells 
n1 = 1; % Min Number of Unit cells

m = m1 / 2;
n = n1 / 2;

k1 = (m - n) / (Zmax - Zmin);
c0 = (k1 * Zmin * Zmin) / 2;
c1 = -(Zmin * k1) + n;

% Uniform scaling function applied only to Z
g = (k1 * z + c1);
u = sin(x).*cos(y) + sin(y).*cos(z) + sin(z).*cos(x);

% Sheet
S = (u + c) .* (u - c);

% Compute isosurfaces and caps
[F1, V1] = isosurface(x, y, z, S, 0);
[F2, V2] = isocaps(x, y, z, S, 0, 'below');


% Combine faces and vertices
F3 = [F1; F2 + size(V1, 1)];
V3 = [V1; V2];

% Create the patch object for the visualization
% P = patch('Vertices', V3, 'Faces', F3, 'FaceColor', 'red', 'EdgeColor', 'none');
% axis equal;


% Surface Area Calculation
SA = 0;
for i = 1:size(F3, 1)
    v1 = V3(F3(i, 1), :);
    v2 = V3(F3(i, 2), :);
    v3 = V3(F3(i, 3), :); 
    edge1 = v2 - v1;
    edge2 = v3 - v1;
    area = 0.5 * norm(cross(edge1, edge2));
    SA = SA + area;
end
k = reshape(V3(F3,:),[size(F3) 3]);
k2= V3(F3,:);

% Volume Calculation
VF = permute(reshape(V3(F3,:),[size(F3) 3]),[3 1 2]);
Vol = 1/6*sum(dot(cross(VF(:,:,1),VF(:,:,2),1),VF(:,:,3),1));

plotg(1,loop+1) = c;
plotg(2,loop+1) = Vol/(2*2*2*pi*pi*pi);
plotg(3,loop+1) = SA;
loop = loop + 1;
end

figure;
axis equal;
plot(plotg(2,:), plotg(1,:), '*', 'MarkerSize', 7,'Color', 'red');
axis([0 1 0 5])
ylabel('Constant value, c');
xlabel('Relative Density, \rho');
title('Relationship of Relative Density and Constant value');
legend('Gyroid','c(\rho)=67.8793(\rho)^8-237.17(\rho)^7+326.72(\rho)^6-221.16(\rho)^5+74.43(\rho)^4-9.23(\rho)^3-0.704(\rho)^2+1.742(\rho)+0');
hold on
grid on
hold on
n = 3/0.1+1;
m = 8;
for i = 1:2*m
    xsum(i) = sum(plotg(2,:).^(i));
end
a(1,1) = n;
b(1,1) = sum(plotg(1,:));
for j = 2:m+1
    a(1,j) = xsum(j-1);
end
for i = 2:m+1
    for j = 1:m+1
        a(i,j) = xsum(j+i-2);
    end
    b(i,1) = sum(plotg(2,:).^(i-1).*plotg(1,:));
end
p = (a\b)';
for i = 1:m+1
    Pcoef(i) = p(m+2-i);
end

epsilon = 0:0.01:2;
cfit1 = polyval(Pcoef,plotg(2,:));
plot(plotg(2,:), cfit1,'LineWidth', 2);
legend('Gyroid','c(\rho)=67.8793(\rho)^8-237.17(\rho)^7+326.72(\rho)^6-221.16(\rho)^5+74.43(\rho)^4-9.23(\rho)^3-0.704(\rho)^2+1.742(\rho)+0');
SStot = sum((plotg(1,:) - mean(cfit1)).^2); % Total sum of squares
SSres = sum((plotg(1,:) - cfit1).^2); % Residual sum of squares
R_squared = 1 - (SSres / SStot); % R-squared value
% Add R^2 annotation
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', ['R^2 = ', num2str(R_squared)], ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white');
ylabel('Constant value, c');
xlabel('Relative Density, \rho');
title('Relationship of Relative Density and Constant value');