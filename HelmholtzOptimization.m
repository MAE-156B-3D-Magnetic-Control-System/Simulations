clear all
clc
close all

V=12;
mu0=4*pi*10^-7; % permeability of air
Rad=linspace(0.025,0.1,50); % solenoid radii in meters
nTurns=linspace(1,150,50); % current through wires, Amps

% separationDist = linspace(0.025,.5,50); % cycle through separtion distances

z_axial = linspace(-.4,.4,10);

%B = zeros(length(z_axial), length(separationDist));
%B_left = B;
%B_right = B;

Field=zeros(50);
 for i=1:length(nTurns) % varies the number of turns in the coil
     for j=1:length(Rad) % this for loop varies the radius of the solenoid
         I(i,j)=V/Resistance(Rad(j),nTurns(i)); % I=V/Resistance
         if I(i,j)>=20
             I(i,j)=20;
             B=0.5*mu0*nTurns(i)*I(i,j)*(Rad(j)^2)*((Rad(j)^2+(z_axial+(Rad(j)/2)).^2).^(-3/2)+(Rad(j)^2+(z_axial-Rad(j)/2).^2).^(-3/2));
             Field(i,j)=B(5);
             hold on
         else 
             B=0.5*mu0*nTurns(i)*I(i,j)*(Rad(j)^2)*((Rad(j)^2+(z_axial+(Rad(j)/2)).^2).^(-3/2)+(Rad(j)^2+(z_axial-Rad(j)/2).^2).^(-3/2));
             Field(i,j)=B(5);             
             hold on
         end
     end
 end

figure(1)
surf(Rad,nTurns,Field)
colorbar
ylabel('Number of Turns')
xlabel('Radius [m]')
zlabel('Magnetic Field Density [Tesla]]')
title('Optimization of Radius and Number of Turns')



