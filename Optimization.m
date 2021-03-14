clear all
clc
close all

%% Section 1: Helmholtz Coil/Solenoit Optimization
V=12; %Voltage of power supply; can be adjusted to 14 V
mu0=4*pi*10^-7;% permeability of air

Rad=linspace(0.02,0.1,50); % radii in meters
nTurns=linspace(1,150,50); %number of coil turns

separationDist = 0.1328;% 0.055;%%separation distance between th two coils
z_axial = 0; %linspace(-.4,.4,10);%linspace(-.4,.4,41); %was .-4 to .4 before linspace(-.4,.4,10);%

Field=zeros(50);

MaxCurrent=20;
 for i=1:length(nTurns) % varies the number of turns in the coil
     for j=1:length(Rad) % this for loop varies the radius of the solenoid
         I(i,j)=V/Resistance(Rad(j),nTurns(i));% I=V/Resistance
         if I(i,j)>MaxCurrent
             I(i,j)=MaxCurrent;
             B=0.5*mu0*nTurns(i)*I(i,j)*(Rad(j)^2)*((Rad(j)^2+(z_axial+(separationDist/2)).^2).^(-3/2)+(Rad(j)^2+(z_axial-separationDist/2).^2).^(-3/2));
             Field(i,j)=B;
             hold on
         else 
             B=0.5*mu0*nTurns(i)*I(i,j)*(Rad(j)^2)*((Rad(j)^2+(z_axial+(separationDist/2)).^2).^(-3/2)+(Rad(j)^2+(z_axial-separationDist/2).^2).^(-3/2));
             Field(i,j)=B;             
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
title('             Optimization of Radius and Number of Turns when Separation Distance=132mm') 

%%
Field=zeros(50);
 for i=1:length(nTurns) % varies the number of turns in the coil
     for j=1:length(Rad) % this for loop varies the radius of the solenoid
         I(i,j)=V/Resistance(Rad(j),nTurns(i)); % I=V/Resistance
         if I(i,j)>MaxCurrent
             I(i,j)=MaxCurrent;
             B=0.5*mu0*nTurns(i)*I(i,j)*(Rad(j)^2)*((Rad(j)^2+(z_axial+(Rad(j)/2)).^2).^(-3/2)+(Rad(j)^2+(z_axial-Rad(j)/2).^2).^(-3/2));
             Field(i,j)=B;
             hold on
         else 
             B=0.5*mu0*nTurns(i)*I(i,j)*(Rad(j)^2)*((Rad(j)^2+(z_axial+(Rad(j)/2)).^2).^(-3/2)+(Rad(j)^2+(z_axial-Rad(j)/2).^2).^(-3/2));
             Field(i,j)=B;             
             hold on
         end
     end
 end

figure(2)
surf(Rad,nTurns,Field)
colorbar
ylabel('Number of Turns')
xlabel('Radius [m]')
zlabel('Magnetic Field Density [Tesla]]')
title('Optimization of Radius and Number of Turns for True Helmholtz Coil') %separation distance=radius

%% Section 2: C Shape Magnet Calculation
% https://www.youtube.com/watch?v=4UFKl9fULkA&fbclid=IwAR1hIaDQsot-R8Fp2-mJGmyal9paKnfbvzH6iq5BSfTiHpJCsAZ-LJ9RbjA
Voltage=12; %can modify power supply to produce 14V
Radius=0.005;
NumberOfTurns=75;
Current=Voltage/Resistance(Radius,NumberOfTurns);
MagnetomotiveForce=Current*NumberOfTurns; %units of amp*turns
mu0=4*pi*10^-7;% permeability of air
mu0Fe=5000; %5000 for 99.8% pure Iron; 200000 for 99.95% pure iron; values taken from https://en.wikipedia.org/wiki/Permeability_(electromagnetism) 

% R=L/(mu0*A)where R is magnetic reluctance in this case
LGap=0.0207;
LCore=pi*(0.0762)-LGap; %circumference of system
AreaCore=0.01^2;
RCore=LCore/(mu0Fe*mu0*AreaCore);
AGap=(0.0254+LGap)^2;
RGap=LGap/(mu0*AGap);
RTot=RCore+RGap;

% Phi=MagnetomotiveFore/RTot; where Phi is field flux
Phi=MagnetomotiveForce/RTot;

% FluxDensity=Phi/A; %want this maximized
FluxDensity=Phi/AGap;



%% Section 3: C Shape Magnet Optimization
Gauge=[10,12,14,16,18,20,22,24]; %Gauge of copper wire
NumberOfTurns=linspace(1,150,50); %number of coil turns
Voltage=12; %can modify power supply to produce 14V
mu0=4*pi*10^-7;% permeability of air
mu0Fe=5000; %5000 for 99.8% pure Iron; 200000 for 99.95% pure iron; values taken from https://en.wikipedia.org/wiki/Permeability_(electromagnetism) 

% R=L/(mu0*A)where R is magnetic reluctance in this case
LGap=0.0207;
LCore=pi*(0.0762)-LGap; %circumference of system
AreaCore=0.01^2;
RCore=LCore/(mu0Fe*mu0*AreaCore);
AGap=(0.0254+LGap)^2;
RGap=LGap/(mu0*AGap);
RTot=RCore+RGap;

for i=1:length(Gauge)
    for j=1:length(NumberOfTurns)
        Current(i,j)=Voltage/Resistance2(Gauge(i),NumberOfTurns(j));
        if Current(i,j)>=MaxCurrent
            Current(i,j)=MaxCurrent;
            MagnetomotiveForce(i,j)=Current(i,j)*NumberOfTurns(j); %units of amp*turns
            Phi(i,j)=MagnetomotiveForce(i,j)/RTot; % Phi=MagnetomotiveFore/RTot; where Phi is field flux
            FluxDensity(i,j)=Phi(i,j)/AGap;
        else
            MagnetomotiveForce(i,j)=Current(i,j)*NumberOfTurns(j); %units of amp*turns
            Phi(i,j)=MagnetomotiveForce(i,j)/RTot; % Phi=MagnetomotiveFore/RTot; where Phi is field flux
            FluxDensity(i,j)=Phi(i,j)/AGap;
        end
    end
end

figure(3)
surf(NumberOfTurns,Gauge,FluxDensity)
colorbar
ylabel('WireGauge')
xlabel('Number of Turns')
zlabel('Magnetic Field Density [Tesla]]')
title('  Optimization of Wire Gauge and Number of Turns for C-Shaped Magnet')

%% Section 4: Magnetic Force Calculation
% Equation for solenoid magnetic force from https://sciencing.com/calculate-magnetic-force-solenoid-6310220.html
n=150; %n represents number of coil turns
current=12; %current running through the solenoid in amps
mu0=4*pi*10^-7; % permeability of air
SolenoidRadius=0.055; % solenoid radius in meters
SolenoidArea=pi*SolenoidRadius^2;
gapLength= 0.06604; %SolenoidRadius/2;

SolenoidMagneticForce=((n*current)^2)*mu0*(SolenoidArea/(2*gapLength^2));





% Example Verification
% n=1000; %n represents number of coil turns
% current=10; %current running through the solenoid in amps
% mu0=4*pi*10^-7; % permeability of air
% SolenoidRadius=0.055; % solenoid radius in meters
% SolenoidArea=0.5;
% gapLength=1.5;
% 
% SolenoidMagneticForce=((n*current)^2)*mu0*(SolenoidArea/(2*gapLength^2));



























