function ResistanceFromCoil = Resistance(Radius,nTurns)
%RESISTANCE Computes the resistance in solenoid based on the number of
%wrappings used. This function is written with values for 16 AWG copper
%wire from http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/wirega.html
%   Detailed explanation goes here
OhmsPerMeter=13.2/1000*3; %13.2 ohms per 1000 meters
Circumference=Radius*2*pi;
WireUsed=Circumference*nTurns;
ResistanceFromCoil=WireUsed*OhmsPerMeter;
end

