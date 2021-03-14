clear all
clc
close all

%% Define parameters for the helmholtz coil
mu0_Cu=0.999994; %relative permeability of copper
mu0=4*pi*10^-7; %permeability of free space
mu0=4*pi*10^-7;
%n %number of coils
R=0.055; %radius of helmholtz coils
I=5; %current from a laptop power source %need to see what the charger was rated at

%% Helmholtz coil equation
% Equations from https://en.wikipedia.org/wiki/Helmholtz_coil
n=linspace(10,1000,1000); %number of coils
x=linspace(0,0.05,50);


n2=50;
B1=0.5*mu0.*n*I*(R^2)*((R^2+(0+(R/2)).^2).^(-3/2)+(R^2+(0-R/2).^2).^(-3/2));
%B1=(4/5)^(3/2)*((mu0*n*I)/R); %determines field density at the center of the Helholtz coil
%B2=(4/5)^(3/2)*((mu0*n2*I)/R);
%B2=(mu0*n2*I*R^2)./(2*(R^2+x.^2).^(3/2)); %idealizes the system as a having n loops
B2=0.5*mu0*n2*I*(R^2)*((R^2+(x+(R/2)).^2).^(-3/2)+(R^2+(x-R/2).^2).^(-3/2));

figure(1)
sgtitle('Helmholtz Coil Simulations','FontSize')%, 30)
subplot(3,1,1)
plot(n,B1)
hold on
% yline(1)
title('Field Density vs Number of Coils at Center of Helmholtz Coil')
xlabel('Number of Coils')
ylabel('Field Density (Tesla)')
% ax = gca;
% ax.FontSize = 18;

xplot=linspace(0,50,50);
%figure(2)
subplot(3,1,3) %this plot might be garbage
plot(xplot,B2)
hold on
%yline(1)
title('Field Density Inside of Helmholtz Coils')
xlabel('Distance from Center of Helmholtz Coil (mm)')
ylabel('Field Density (Tesla)')
xlim([0,50])
legend('N=50 Coil Turns and 5 Amps')
% ax = gca;
% ax.FontSize = 18;

I3=linspace(0.1,20,100);
B3=(4/5)^(2/3)*((mu0*n2*I3)/R);
subplot(3,1,2)
plot(I3,B3)
hold on
%yline(1)
title('Field Density at Center of Coil vs Current at Center of Helmholtz Coil')
xlabel('Electrical Current (amps)')
ylabel('Field Density (Tesla)')
% ax = gca;
% ax.FontSize = 18;

%% Iron core parameters
% R2=0.0508; %radius of iron core
% L=0.0254; %length of iron core
% V=pi*R^2*L; %Volume of iron core
% mu0_Fe=1150; %magnetic permeability of iron
% M=1; %need an actual number here! M is supposed to be the magnetization of the magnets
% x2=linspace(0.01,0.05,500);
% 
% 
% % our iron core is 99.9% iron; Darren got it from China
% %% Iron Core Calculations
% % Equation from https://en.wikipedia.org/wiki/Force_between_magnets
% % Uses equation for the force between two cylindrical magnets
% F=((pi*mu0_Fe)/4)*(M^2)*(R2^4)*((1./(x2.^2))+(1./((x2+2*L).^2))-(2./((x2+L).^2)));
% 
% figure(1)
% subplot(2,2,4)
% plot(x2,F)
% title('Magnetic Force vs Distance Between Magnets')
% xlabel('Distance (meters)')
% ylabel('Magnetic Force (Newtons)')


%% Simulating Solenoids
n3=50;
I3=5;
mu0_Fe=6.3*10^-3;
R3=0.055;
mu0=4*pi*10^-7;
%D3=.1;

z=linspace(-.2,.2,50);
%D3=linspace(0.03,.3,10);
D3=0.055;

B4=0.5*mu0*n3*I3*(R3^2)*((R3^2+(z+(D3/2)).^2).^(-3/2)+(R3^2+(z-D3/2).^2).^(-3/2));

zplot=linspace(-.3,.3,50);
figure(3)
plot(zplot,B4)
title(['Field Density vs Distance in Helmholtz Coil at Distance 55 mm'])
xlabel('Distance from Center of Helmholtz Coil (mm)')
ylabel('Field Density (Tesla)')

Test1=B4(25)
%%
n3=50;
I3=5;
mu0_Fe=6.3*10^-3;
R3=0.08;
mu0=4*pi*10^-7;
%D3=.1;

z=linspace(-.2,.2,50);
%D3=linspace(0.03,.3,10);
D3=0.055;

B4=0.5*mu0*n3*I3*(R3^2)*((R3^2+(z+(D3/2)).^2).^(-3/2)+(R3^2+(z-D3/2).^2).^(-3/2));

zplot=linspace(-.3,.3,50);
figure(3)
plot(zplot,B4)
title(['Field Density vs Distance in Helmholtz Coil at Distance 55 mm'])
xlabel('Distance from Center of Helmholtz Coil (mm)')
ylabel('Field Density (Tesla)')

Test2=B4(25)
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
%% 
% %% Magnetic Field Strengh Between two Helmholtz Coils
% % Equation from https://virtuelle-experimente.de/en/b-feld/b-feld/helmholtzspulenpaar.php
% 
% I3=50; % get this current using the follwoing link
% % https://www.amazon.com/Supply%EF%BC%88SMPS%EF%BC%89-Monitoring-Industrial-Transformer-220VAC-DC12V/dp/B0781V1D7Q/ref=sr_1_4?dchild=1&keywords=50+amp+power+supply&qid=1607408251&sr=8-4
% % https://www.amazon.com/PowerMax-PM4-100A-Converter-Battery/dp/B01ER3LH5W?th=1
% N1=linspace(100,1000,100); %number of windings
% R3=0.1; %Radius of helmholtz coils and distance between the two coils
% 
% B=mu0*((8*I3*N1)/(sqrt(125)*R3));
% 
% figure(3)
% plot(N1,B)
% title('Field Density Between Two Helmholtz Coils')
% xlabel('Number of Coils')
% ylabel('Magnetic Field Density (Teslas')


%% 


B5=0.5*mu0*I3*R3^2./(2*(R3^2+z.^2).^3/2);









