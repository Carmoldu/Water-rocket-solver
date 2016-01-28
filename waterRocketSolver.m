function [h, Thrust, Acc, Mwater, Mair, Vexhaust, massFlow,P,Weight, t ] = waterRocketSolver( Vbottle, Vwater0, nozzleDiam, bottleDiam, P0, Patm, Cd, Mpl, Mstruc, gamma)
%This function takes data from the rocket design and returns its
%performance
%   INPUT
% Vbottle: bottle volume (L)
% Vwater: water volume (L)
% nozzleDiam: diameter of the nozzle (cm)
% bottleDiam: diameter of the bottle (cm)
% P0: initial internal pressure (bar)
% Patm: atmospheric pressure (bar)
% Cd: Drag coefficient of the rocket for frontal axial airstream
% Mpl: payload mass (g)
% Mstruc: structure mass (g)
% gamma: adiabatic constant
%
%   OUTPUT
% h: height profile matrix (m)
% Thrust: thrust profile matrix (N)
% Acc: acceletation profile matrix (m/s2)
% Mwater: mass of water left in the rocket (kg)
% Mair: mass of air left in the rocket (kg)
% Vexhaust: exhaust velocity of the propellants (m/s)
% massFLow: mass flow (kg/s)
% P: internal pressure profile (Pa)
% t: time matrix of the simulation (s)

%constants
g=9.80665; %m/s2
rhoWater=1000;  %Kg/m3
rhoAir=1.225;   %Kg/m3

LitersToM3=0.001; %m3/L
CmToM=0.01; %m/cm
GrToKg=0.001; %Kg/g
BarToPa=10e4; %Pa/bar

%Simulation constants
deltaT=0.1;  %s
h0=0;   %m
V0=0;   %m/s
t0=0;   %s

%Conversions and precomputations
Vbottle=Vbottle*LitersToM3 %m3
Vwater0=Vwater0*LitersToM3 %m3

nozzleDiam=nozzleDiam*CmToM    %m
bottleDiam=bottleDiam*CmToM;   %m

AreaNozzle=pi*(nozzleDiam/2)^2; %m^2
AreaBottle=pi*(bottleDiam/2)^2; %m^2

P0=BarToPa*P0        %Pa
Patm=BarToPa*Patm    %Pa

Mpl=GrToKg*Mpl       %Kg
Mstruc=GrToKg*Mstruc %Kg


%INITIAL STATE t=0: The state is computed for instant 0 as the variables
%will be considered constant during the first step and used to compute the
%consecuent states
V(1)=V0;
h(1)=h0;
t(1)=t0;
P(1)=P0;

Vexhaust(1)=sqrt((P0-Patm+rhoWater*g*Vwater0/AreaBottle)/(0.5*rhoWater*(1-(AreaNozzle/AreaBottle)^2)));     %m/s

massFlow(1)=rhoWater*AreaNozzle*Vexhaust(1);  %kg/s

Thrust(1)=Vexhaust(1)*massFlow(1);               %N

Mair(1)=Vbottle*rhoAir;         %kg  it must be highlighted that the mass of air present in the bottle is the same as previous to compresion, thus the mass for the full bottle volume at ambient temperature density
Mwater(1)=Vwater0*rhoWater;               %kg       As water is basically uncompressible, we can use a standard value for its density without introducing much error
TotalMass(1)=Mair(1)+Mwater(1)+Mpl+Mstruc;     %kg

Weight(1)=TotalMass(1)*g;                     %N

NetForce(1)=Thrust(1)-Weight(1);                 %N

Acc(1)=NetForce(1)/TotalMass(1);              %m/s2

i=1;
figure
plot(t(1),Mwater(1),'b',t(1),Mair(1),'r')
hold on
title('Water and air mass remaining vs time')

while h(i)>=0
    t(i+1)=t(i)+deltaT;
    V(i+1)=V(i)+Acc(i)*deltaT;                %m/s
    h(i+1)=h(i)+V(i)*deltaT+Acc(i)*deltaT^2 ; %m
    

    if Mwater(numel(Mwater))>0              %This will compute the water-propelled stage until "water-out"
        Mwater(i+1)=Mwater(i)-massFlow(i)*deltaT;   %Kg
        Vwater=Mwater(i+1)/rhoWater;                %m3

        P(i+1)=P0*((Vbottle-Vwater0)/(Vbottle-Vwater))^gamma;     %Pa  %Instead of calculating the mid point solution with an iterative algorithm (such as Newton-Rhapson) to solve P and Vexhaust, we will take the previous node solution

        hwater=Vwater/AreaBottle;       %m

        Vexhaust(i+1)=sqrt((P(i+1)-Patm+rhoWater*g*hwater)/(0.5*rhoWater-(AreaNozzle/AreaBottle)^2));   %m/s

        massFlow(i+1)=rhoWater*AreaNozzle*Vexhaust(i+1);        %Kg/s we suppose only water exits the bottle at this stage
        
        Mair(i+1)=Mair(i);         %kg  As only water exits the bottle, this keeps still
        
 
        plot(t(i),Mwater(i),'b',t(i),Mair(i),'r')
        hold on
        drawnow
        
%      elseif P(numel(P))>Patm               %This will compute the air-propelled stage until the pressure inside the bottle is equal to the atmospheric pressure
%          
%          P(i+1)=
%          
%          Vexhaust(i+1)=sqrt((P(i+1)-Patm)/(0.5*rhoAir-(AreaNozzle/AreaBottle)^2));   %m/s   with air, we don't consider the contribution of the column of air contained inside the column
%          massFlow(i+1)=rhoAir*AreaNozzle*Vexhaust(i+1);        %Kg/s (now only air exits the bottle)
%          
%          Mair(i+1)=Mair(i)-massFlow*deltaT;
%          Mwater(i)=0;
         
     else            %this will compute the free fall stage (no thrust)
         Vexhaust(i+1)=0;
         massFlow(i+1)=0;
         
         Mwater(i+1)=0;
         Mair(i+1)=0;
         
         P(i+1)=Patm;
         
     end
    
    
    Thrust(i+1)=Vexhaust(i+1)*massFlow(i+1);

    TotalMass(i+1)=Mair(numel(Mair))+Mwater(numel(Mwater))+Mpl+Mstruc;

    Weight(i+1)=TotalMass(i+1)*g;                     %N
    
    Drag(i+1)=0.5*Cd*rhoAir*AreaBottle*abs(V(i+1))*V(i+1);

    NetForce(i+1)=Thrust(i+1)-Weight(i+1)-Drag(i+1);                 %N

    Acc(i+1)=NetForce(i+1)/TotalMass(i+1);    
    
    i=i+1;
end

figure
plot(t,h,'b', t,V,'r')
title('Height & Velocity vs time')
legend('Height','Velocity')

figure
plot(t,Acc,'b',t,Vexhaust,'r')
title('Acceleration & exhaust velocity vs time')
legend('Acceleration','Exhaust velocity')

figure
plot(t,Thrust)
title('Thrust vs time')

figure
plot(t,Weight)
title('Weight vs time')

figure
plot(t,P)
title('Pressure vs time')



end

