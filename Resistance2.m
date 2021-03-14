function ResistanceFromCoil = Resistance2(Gauge,nTurns)
%RESISTANCE Computes the resistance in solenoid based on the number of
%wrappings used. This function is written with values for 16 AWG copper
%wire from http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/wirega.html
%   Detailed explanation goes here
Side=0.0254; %1 inch sides
switch(Gauge)
    case(24)
        OhmsPerMeter=84.2/1000; %ohms per 1000 meters
        Perimeter=4*Side;
        WireUsed=Perimeter*nTurns;
        ResistanceFromCoil=WireUsed*OhmsPerMeter;
    case(22)
        OhmsPerMeter=52.7/1000; %ohms per 1000 meters
        Perimeter=4*Side;
        WireUsed=Perimeter*nTurns;
        ResistanceFromCoil=WireUsed*OhmsPerMeter; 
    case(20)
        OhmsPerMeter=33.2/1000; %ohms per 1000 meters
        Perimeter=4*Side;
        WireUsed=Perimeter*nTurns;
        ResistanceFromCoil=WireUsed*OhmsPerMeter;
    case(18)
        OhmsPerMeter=20.9/1000; %ohms per 1000 meters
        Perimeter=4*Side;
        WireUsed=Perimeter*nTurns;
        ResistanceFromCoil=WireUsed*OhmsPerMeter;
    case(16)
        OhmsPerMeter=13.2/1000; %ohms per 1000 meters
        Perimeter=4*Side;
        WireUsed=Perimeter*nTurns;
        ResistanceFromCoil=WireUsed*OhmsPerMeter;
    case(14)
        OhmsPerMeter=8.28/1000; %ohms per 1000 meters
        Perimeter=4*Side;
        WireUsed=Perimeter*nTurns;
        ResistanceFromCoil=WireUsed*OhmsPerMeter;
    case(12)
        OhmsPerMeter=5.21/1000; %ohms per 1000 meters
        Perimeter=4*Side;
        WireUsed=Perimeter*nTurns;
        ResistanceFromCoil=WireUsed*OhmsPerMeter;
    case(10)
        OhmsPerMeter=3.28/1000; %ohms per 1000 meters
        Perimeter=4*Side;
        WireUsed=Perimeter*nTurns;
        ResistanceFromCoil=WireUsed*OhmsPerMeter;
    otherwise
        fprintf('error\n' )
end

end

