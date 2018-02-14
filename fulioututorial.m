%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running Fu-Liou radiative transfer code in Matlab.
% Angie Pendergrass, 2/3/2014
% Things you need: 
%   This script: fulioututorial.m 
%   The data file with a profile in it: inputprofile.mat
%   A script that does more heavy lifting: rad_trans_mark.m 
%   Fu-Liou code: jacket.f, dynamic.f, Makefile 
%   A computer with Matlab and the pgf77 fortran compiler (e.g., challenger, olympus, pynchon) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This tutorial will hopefully get you up and running the Fu-Liou radiative transfer code. I have included a profile that should work.  First there is an example of a clear-sky calculation, and a second example includes clouds too. 

load inputprofile.mat za ts p w ta o3  % This profile should work.  We will go through the components below.  

za % this is the solar zenith angle in degrees used for the calculation.  Bryce Harrop has a program that calculates the time-averaged solar-weighted zenith angle which is useful sometimes.  A zenith angle of 0 is directly overhead.  The incoming solar radiation can be changed by combining the zenith angle with the solar constant (which you have to change in the fortran code itself). 

ts % this is the surface temperature in Kelvin. 

p  % this is a pressure profile in mb.  It is important that the order of this (and all other profile variables) is surface first, top last.  I started using the surface pressure for the first pressure level, along with surface air temperature for the first temperature value (tas in CMIP5 archives). 

ta % air temperature profile. The temperatures should correspond to the pressure levels in p.  

w % Water vapor mixing ratio in kg/kg. 

o3 % ozone profile. mixing ratio in kg/kg.  Its easy to be very far off on these.  You should do a sanity check on the column SW absorption to make sure its reasonable after running the code ( SWdn(end)-SWdn(1)- (SWup(end)-SWup(1))). It should be 50 Wm2 with the profile Ive given you and a solar constant of 375.  

co2=348;       % CO2 concentration in ppm.

compileflag=0; % This is a flag that tells the other matlab script to recompile the fortran code.  You should set it to 0 whenever you have a profile with a different number of levels than the previous one, or if you change anything in the fortran code.  

dummy=zeros(size(p)); % This is a dummy variable with all zeros that we will use to make clear-sky calculations.  We put it where the cloud variables would go. 

% This next command calls a script that does the heavy lifting.  Inputs are your profiles (the dummy variables are for cloud fields; set them to zero for clear-sky calculations).   The script writes your profile to the file CURRENT_PROFILE.TXT to be read by the FORTRAN fu-liou program.  If the compile flag is set to zero, it also writes para.file which tells it how many levels your profile has, and recompiles the code.  

%%%%% THIS IS WHERE IT ALL HAPPENS
[SWdn,SWup,LWdn,LWup,dT]=rad_trans_mark(za,ts,p,ta,w,o3,dummy,dummy,dummy,dummy,compileflag,co2);

% If nothing went wrong, we should have 5 output variables. These are also listed in jacket.f

SWdn % Downwelling SW radiation at pressure levels p in Wm^-2.  The first is surface, last is TOA.
SWup % Upwelling SW, Wm-2
LWdn % Downwelling LW, Wm-2 
LWup % Upwelling LW, Wm-2
dT % Net heating rate, K/day

SWdn(end) % Incoming SW at TOA
SWdn(end)-SWup(end) % Absorbed SW, total
SWdn(1)-SWup(1) % SW surface absorptio
LWdn(end) % Outgoing LW radiation (OLR)

% You can go into the code (mainly jacket.f, but also dynamic.f) and change more options.  The main one you will want to change is the solar zenith angle, ss. Others might be ee, the surface LW emissivity, as, the solar surface albedo, umch4, the methane concentration, and umn2o, the N20 concentration.  jacket.f also describes all the inputs to the program and gives their units, so you should take a look at it.  


% This loop just calculates heights using the thickness equation, because for some cloud calculations we have to work in meters instead of mb. 
L=2.5*10^6; 
R=287;g=9.81;cp=1004;Rv=461.5;
z=zeros(size(p));
for i=2:length(p)
   z(i)=-R.*ta(i-1)/g.*log(p(i)./p(1));
end

%%% you can uncomment these lines and they should show the profile of pressure versus height
%clf
%plot(p,z/1000,'-ko');
%set(gca,'xdir','rev')
%set(gca,'ytick',[0:1:15],'xtick',flipud(p));ylim([0 15]);grid on;xlabel('Pressure (mb)');
%ylabel('Height (km)')

dz=diff(z); % height of each layer

i1=2; % this is the bottom cloud layer
i2=3; % this is the top cloud layer
dzh=dz((i1-1):(i2-1));
dzt=sum(dzh); % this is the total cloud height
lwc=4.0./dzt;  % this is the liquid water content you need to get a liquid water path of 4 g/m2 over the whole cloud. 

lwcp=dummy;
lwcp(i1:i2)=lwc; %% this sets the liquid water content of the cloud layers in g/m3
rliq=dummy;
rliq(i1:i2)=14; %% effictive radius of 14 micrometers in all the cloud layers

% You can do similar things for ice clouds: specify the ice water content and effective radius of ice. Experiment on your own! 

% Here is the call to radiaiton code including the cloud. 
compileflag=1 % We dont need to recompile because we have the same number of levels and before. 
[SWdncld,SWupcld,LWdncld,LWupcld,dTcld]=rad_trans_mark(za,ts,p,ta,w,o3,lwcp,dummy,rliq,dummy,compileflag,co2);

% There you go.  Have fun! 
