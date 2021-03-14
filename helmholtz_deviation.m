%% Joseph Martin 1/11/21
%% Tanner Hanson 1/30/21

% Helmholtz deviation simulation
% Ideal helmholtz: B = 0.5*mu0_Fe*n*I*(R^2)*((R^2+(z+(R/2)).^2).^(-3/2)+(R^2+(z-R/2).^2).^(-3/2));

%% Simulating Solenoids

clear all
clc
close all

% nTurns = 50;
V=15;
mu0=4*pi*10^-7;
%I = 5; % current through wires, Amps
% R = 0.055; % solenoid radius, meters

Rad=linspace(0.025,0.1,10); % radii in meters
nTurns=linspace(1,1500,10);

mu0 = 4 * pi * 10^-7; % permeability of air


% separationDist = linspace(0.025,.5,50); % cycle through separtion distances

z_axial = linspace(-.4,.4,10);

%B = zeros(length(z_axial), length(separationDist));
%B_left = B;
%B_right = B;
count=0;
% I=V/Resistance

 for i=1:length(nTurns) % varies the number of turns in the coil
     for j=1:length(Rad) % this for loop varies the radius of the solenoid
         I(i,j)=V/Resistance(Rad(j),nTurns(i));
        % B=0.5*mu0*nTurns(i)*I(i,j)*(Rad(j)^2)*((Rad(j)^2+(z_axial+(Rad(j)/2)).^2).^(-3/2)+(Rad(j)^2+(z_axial-Rad(j)/2).^2).^(-3/2));
       %  Field(i,j)=B(500);
         hold on
     end
 end
 Rad=[Rad;Rad;Rad;Rad;Rad;Rad;Rad;Rad;Rad;Rad]
 nTurns=[nTurns;nTurns;nTurns;nTurns;nTurns;nTurns;nTurns;nTurns;nTurns;nTurns]'
 z_axial=zeros(10);%[z_axial;z_axial;z_axial;z_axial;z_axial;z_axial;z_axial;z_axial;z_axial;z_axial]
B=0.5.*mu0.*nTurns.*I.*(Rad.^2).*((Rad.^2+(z_axial+(Rad./2)).^2).^(-3/2)+(Rad.^2+(z_axial-Rad./2).^2).^(-3/2));
figure(1)
surf(Rad,nTurns,B)%Field)
colorbar
ylabel('Number of Turns')
xlabel('Radius [m]')
zlabel('Magnetic Field Density [Tesla]]')


%% Use this simulation version
% nTurns = 50;
V=12;
mu0=4*pi*10^-7;
%I = 5; % current through wires, Amps
% R = 0.055; % solenoid radius, meters

Rad=linspace(0.025,0.1,50); % radii in meters
nTurns=linspace(1,150,50);

mu0 = 4 * pi * 10^-7; % permeability of air

% separationDist = linspace(0.025,.5,50); % cycle through separtion distances

z_axial = linspace(-.4,.4,10);

%B = zeros(length(z_axial), length(separationDist));
%B_left = B;
%B_right = B;
% I=V/Resistance
Field=zeros(50);
 for i=1:length(nTurns) % varies the number of turns in the coil
     for j=1:length(Rad) % this for loop varies the radius of the solenoid
         I(i,j)=V/Resistance(Rad(j),nTurns(i));
         if I(i,j)>20
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

%B=0.5.*mu0.*nTurns.*I.*(Rad.^2).*((Rad.^2+(z_axial+(Rad./2)).^2).^(-3/2)+(Rad.^2+(z_axial-Rad./2).^2).^(-3/2));
figure(1)
surf(Rad,nTurns,Field)
colorbar
ylabel('Number of Turns')
xlabel('Radius [m]')
zlabel('Magnetic Field Density [Tesla]]')



%%
 Rad=0.055;
 nTurns=100;
 z_axial=0;
 I=V/Resistance(Rad,nTurns)
B1=0.5.*mu0.*nTurns.*I.*(Rad.^2).*((Rad.^2+(z_axial+(Rad./2)).^2).^(-3/2)+(Rad.^2+(z_axial-Rad./2).^2).^(-3/2))

 Rad=0.055;
 nTurns=2*10^8;
 z_axial=0;
 I=V/Resistance(Rad,nTurns)
B2=0.5.*mu0.*nTurns.*I.*(Rad.^2).*((Rad.^2+(z_axial+(Rad./2)).^2).^(-3/2)+(Rad.^2+(z_axial-Rad./2).^2).^(-3/2))


%%

clear all
clc
close all

% nTurns = 50;
V=15;
mu0=4*pi*10^-7;
%I = 5; % current through wires, Amps
% R = 0.055; % solenoid radius, meters

Rad=linspace(0.030,0.1,10); % radii in meters
nTurns=linspace(1,150,10);

mu0 = 4 * pi * 10^-7; % permeability of air
%mu0_Fe = 6.3e-3; % permeability of iron

% separationDist = linspace(0.025,.5,50); % cycle through separtion distances
z_axial = linspace(-.4,.4,10);
z_axial=0;
%B = zeros(length(z_axial), length(separationDist));
%B_left = B;
%B_right = B;
count=0;
% I=V/Resistance

 for i=1:length(nTurns) % varies the number of turns in the coil
     for j=1:length(Rad) % this for loop varies the radius of the solenoid
         I(i,j)=V/Resistance(Rad(j),nTurns(i));
         B(i,j)=0.5*mu0*nTurns(i)*I(i,j)*(Rad(j)^2)*((Rad(j)^2+(z_axial+(Rad(j)/2)).^2).^(-3/2)+(Rad(j)^2+(z_axial-Rad(j)/2).^2).^(-3/2));
%nTurns(i)
%Rad(j)
         Field(i,j)=B(i,j);
         hold on
     end
 end

figure(1)
surf(Rad,nTurns,B)
colorbar
ylabel('Number of Turns')
xlabel('Radius [m]')
zlabel('Magnetic Field Density [Tesla]]')





% nTurns = 50;
V=15;
mu0=4*pi*10^-7;
%I = 5; % current through wires, Amps
% R = 0.055; % solenoid radius, meters

Rad=linspace(0.030,0.1,10); % radii in meters
nTurns=linspace(1,150,10);

mu0 = 4 * pi * 10^-7; % permeability of air
%mu0_Fe = 6.3e-3; % permeability of iron

% separationDist = linspace(0.025,.5,50); % cycle through separtion distances
separationDist=ones(10)*0.055;
%z_axial = linspace(-.4,.4,10);
z_axial=0;
%B = zeros(length(z_axial), length(separationDist));
%B_left = B;
%B_right = B;
count=0;
% I=V/Resistance

 for i=1:length(nTurns) % varies the number of turns in the coil
     for j=1:length(Rad) % this for loop varies the radius of the solenoid
         I(i,j)=V/Resistance(Rad(j),nTurns(i));
         B(i,j)=0.5*mu0*nTurns(i)*I(i,j)*(separationDist(j)^2)*((separationDist(j)^2+(z_axial+(separationDist(j)/2)).^2).^(-3/2)+(separationDist(j)^2+(z_axial-separationDist(j)/2).^2).^(-3/2));
%nTurns(i)
%Rad(j)
         Field(i,j)=B(i,j);
         hold on
     end
 end

figure(2)
surf(Rad,nTurns,B)
colorbar
ylabel('Number of Turns')
xlabel('Radius [m]')
zlabel('Magnetic Field Density [Tesla]]')


















% B4=0.5*mu0*n3*I3*(R3^2)*((R3^2+(z+(D3/2)).^2).^(-3/2)+(R3^2+(z-D3/2).^2).^(-3/2));
% 
% zplot=linspace(-.3,.3,50);
% figure(3)
% plot(zplot,B4)
% title(['Field Density vs Distance in Helmholtz Coil at Distance 55 mm'])
% xlabel('Distance from Center of Helmholtz Coil (mm)')
% ylabel('Field Density (Tesla)')
% 
% B4(25)

%%
% 
% %%
% n3=50;
% I3=5;
% mu0_Fe=6.3*10^-3;
% R3=0.055;
% mu0=4*pi*10^-7;
% D3=.1;
% 
% z=linspace(-.2,.2,50);
% D3=linspace(0.03,.3,10);
% for n=1:length(D3)
% %B4=0.5*mu0_Fe*n3*I3*(R3^2)*((R3^2+(z+(R3/2)).^2).^(-3/2)+(R3^2+(z-R3/2).^2).^(-3/2));
% B4=0.5*mu0_Fe*n3*I3*(R3^2)*((R3^2+(z+(D3(n)/2)).^2).^(-3/2)+(R3^2+(z-D3(n)/2).^2).^(-3/2));
% 
% zplot=linspace(-.3,.3,50);
% figure(3)
% plot(zplot,B4)
% title(['Field Density vs Distance in Helmholtz Coil at Distance ' num2str(D3(n))])
% xlabel('Distance from Center of Helmholtz Coil (mm)')
% ylabel('Field Density (Tesla)')
% pause(.5)
% end
% 
% 
% 
% %%
% 
% for i=1:length(separationDist)
% 
%     % Calculate Field Intensity, B
%     B(:,i) = 0.5 * mu0_Fe * nTurns * I * R^2 * ( ...
%         (R^2 + (z_axial + separationDist(i)/2).^2).^(-3/2) + ... % left coil
%         (R^2+(z_axial-separationDist(i)/2).^2).^(-3/2)... % right coil
%         ); 
%     
%     B_left(:,i) = 0.5 * mu0_Fe * nTurns * I * R^2 * ( ...
%         (R^2 + (z_axial + separationDist(i)/2).^2).^(-3/2));
%     
%     B_right(:,i) = 0.5 * mu0_Fe * nTurns * I * R^2 * ( ...
%         (R^2+(z_axial-separationDist(i)/2).^2).^(-3/2)); 
% end
% %%
% % plot
% h = figure();
% 
% zplot = z_axial*1000;
% filename = 'Helmholtz Deviation.gif';
% 
% for i=1:length(separationDist)
%     
%     coilPair = plot(zplot,B(:,i), 'linewidth', 3, 'color', 'black');
%     hold on
%     leftCoil = plot(zplot,B_left(:,i), 'linewidth', 2, 'color', 'blue');
%     rightCoil = plot(zplot,B_right(:,i), 'linewidth', 2, 'color', 'red');
%     title(['Axial field intensity for separation distance: ' num2str(round(1000*(separationDist(i)))) 'mm'], 'fontsize', 18)
%     axis([min(zplot) max(zplot) min(min(B)) max(max(B))]);
%     xlabel('Position, relative to center (mm)');
%     ylabel('Field Density (Tesla)');
%     
%     leftBox = xline(-25);
%     rightBox = xline(25);
%     legend([coilPair, leftCoil, rightCoil], 'Total field dist', 'Left coil Field', 'Right coil field')
%     
%     text(min(zplot)-.05*min(zplot), max(max(B))-.1*max(max(B)), ...
%         ['Coil radius: ' num2str(R*1000) 'mm\newline' ...
%         'Current: ' num2str(I) 'A\newline' ...
%         'Iron core: yes'], 'fontsize', 14)
%     
%     hold off
%     grid on
%     
%     drawnow;
%     
%     f = getframe(h); 
% 
%     im = frame2im(f);
%     [ind,cmap] = rgb2ind(im,256);
% 
%     if i == 1
%         imwrite(ind,cmap,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(ind,cmap,filename,'gif','DelayTime',0.2, 'WriteMode','append');
%     end
%     
% end