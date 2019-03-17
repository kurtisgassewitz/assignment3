%% Assignment 3: Part 2
%
% The Finite Difference method was then used to calculate the the potential
% within a region with a bottle neck inserted. The same approach was taken
% in assignment 2 part 2. The results can be found below. 


close all;
clear;

Sigma = 1;

nx = 50;
ny = 50;

G = sparse (nx*ny, nx*ny);
B = zeros(1, nx*ny);

cMap = zeros (nx, ny);

%Loop to assign Conductivity
for i = 1:nx
    for j = 1:ny
        if ((i>=0.4*nx) && (i<=0.6*nx) && (j<=0.4*ny)) || ((i>=0.4*nx) && (i<=0.6*nx) && (j>=0.6*ny))
            cMap(i,j) = .01;
        else
            cMap(i,j) = Sigma;
        end   
    end
end

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        
        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            B(n)    = 1;
        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
        elseif j == 1
            nxm = j + (i - 2)*ny;
            nxp = j + (i)*ny;
            nyp = j + 1 + (i - 1)*ny;
            
            rxm = (cMap(i,j) + cMap(i - 1, j))/2;
            rxp = (cMap(i,j) + cMap(i + 1, j))/2;
            ryp = (cMap(i,j) + cMap(i, j + 1))/2;
            
            G(n, n) = -(rxm + rxp + ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;
            
        elseif j == ny
            nxm = j + (i - 2)*ny;
            nxp = j + (i)*ny;
            nym = j - 1 + (i - 1)*ny;
            
            rxm = (cMap(i,j) + cMap(i - 1, j))/2;
            rxp = (cMap(i,j) + cMap(i + 1, j))/2;
            rym = (cMap(i,j) + cMap(i, j - 1))/2;
        
            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
            
        else
            
            nxm = j + (i - 2)*ny;
            nxp = j + (i)*ny;
            nym = j - 1 + (i - 1)*ny;
            nyp = j + 1 + (i - 1)*ny;
            
            rxm = (cMap(i,j) + cMap(i - 1, j))/2;
            rxp = (cMap(i,j) + cMap(i + 1, j))/2;
            rym = (cMap(i,j) + cMap(i, j - 1))/2;
            ryp = (cMap(i,j) + cMap(i, j + 1))/2;
        
            G(n, n) = -(rxm + rxp + ryp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
            G(n, nyp) = ryp;
            
        end
        
    end
    
end

V = G\B';

vMap = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        n = j + (i - 1)*ny;
        vMap(i,j) = V(n);
    end
end

vMap_T = vMap';

%% 
% The Finite Difference Method was then applied to solve for the current
% flow through the region. A voltage map was created for the region and can
% be seen in the Figure below. 

%Potential Plot
figure(1)
surf(vMap_T);
title('Voltage');

%% 
% The electric field was then determined for the region. The electric field was found by taking the gradient of the voltage matrix. 
% The electric field of the region can be seen in the Figure below. Note
% that there was an adjustment on the electric field values which would
% account for the true area of the region. 

%Electric Field Plot
[Ex, Ey] = gradient(-vMap_T);
EX = Ex/5E-8;
EY = Ey/5E-8;
figure(2)
quiver(EX, EY);
title ('Electric Field');


