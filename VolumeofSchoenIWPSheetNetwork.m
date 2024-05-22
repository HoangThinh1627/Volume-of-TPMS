clc;
close all;
clear;
loop=1;
for c = 0.1:0.05:5
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

u = 2 * (cos(x) .* cos(y) + cos(y) .* cos(z) + cos(z) .* cos(x)) - (cos(2 * x) + cos(2 * y) + cos(2 * z));
% c = 0.3; % Level Set

% Sheet
S = (u + c) .* (u - c);
S1 =  (u);
S2 = (u + c);

% Compute isosurfaces and caps
[F1, V1] = isosurface(x, y, z, S, 0);
[F2, V2] = isocaps(x, y, z, S, 0, 'below');

% SAn = 0;
% area1 = 0;
% for i = 1:size(Fn, 1)
%     v1 = Vn(Fn(i, 1), :);
%     v2 = Vn(Fn(i, 2), :);
%     v3 = Vn(Fn(i, 3), :); 
%     edge1 = v2 - v1;
%     edge2 = v3 - v1;
%     area1 = 0.5 * norm(cross(edge1, edge2));
%     SAn = SAn + area1;
% end


% Combine faces and vertices
F3 = [F1; F2 + size(V1, 1)];
V3 = [V1; V2];

% Create the patch object for the visualization
P = patch('Vertices', V3, 'Faces', F3, 'FaceColor', 'red', 'EdgeColor', 'none');
axis equal;


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
Vol = 1/6*sum(dot(cross(VF(:,:,1),VF(:,:,2),1),VF(:,:,3),1)); % close to 4/3*pi volume of the sphere of raduius 1

ploti(1,loop) = c;
ploti(2,loop) = Vol/(2*2*2*pi*pi*pi);
ploti(3,loop) = SA;
loop = loop + 1;

end


figure;
axis equal;
plot(ploti(2,:), ploti(1,:), '*', 'MarkerSize', 7,'Color', 'red');
axis([0 1 0 5])
ylabel('Constant value, c');
xlabel('Relative Density, \rho');
title('Relationship of Relative Density and Constant value');
hold on
grid on
hold on
n = 3/0.1+1;
m = 8;
for i = 1:2*m
    xsum(i) = sum(ploti(2,:).^(i));
end
a(1,1) = n;
b(1,1) = sum(ploti(1,:));
for j = 2:m+1
    a(1,j) = xsum(j-1);
end
for i = 2:m+1
    for j = 1:m+1
        a(i,j) = xsum(j+i-2);
    end
    b(i,1) = sum(ploti(2,:).^(i-1).*ploti(1,:));
end
p = (a\b)';
for i = 1:m+1
    Pcoef(i) = p(m+2-i);
end
epsilon = 0:0.01:2;
volfit = polyval(Pcoef,epsilon);
% plot(epsilon, volfit,'LineWidth', 2);
legend('Gyroid','c(\rho)=443.04(\rho)^8-1.26(\rho)^7+1.30(\rho)^6-485.16(\rho)^5-69.42(\rho)^4-102.74(\rho)^3-26.85(\rho)^2+6.06(\rho)+0');
epsilon = 0:0.01:2;
cfit1 = polyval(Pcoef,ploti(2,:));
plot(ploti(2,:), cfit1,'LineWidth', 2,'Color', 'red');
legend('Gyroid','c(\rho)=67.8793(\rho)^8-237.17(\rho)^7+326.72(\rho)^6-221.16(\rho)^5+74.43(\rho)^4-9.23(\rho)^3-0.704(\rho)^2+1.742(\rho)+0');
SStot = sum((ploti(1,:) - mean(cfit1)).^2); % Total sum of squares
SSres = sum((ploti(1,:) - cfit1).^2); % Residual sum of squares
R_squared = 1 - (SSres / SStot); % R-squared value
% Add R^2 annotation
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', ['R^2 = ', num2str(R_squared)], ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white');
ylabel('Constant value, c');
xlabel('Relative Density, \rho');
title('Relationship of Relative Density and Constant value');

