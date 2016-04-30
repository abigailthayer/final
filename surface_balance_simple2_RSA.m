function [bw ,bs ,zela] = surface_balance_simple2(zmid,bwmax,zmid_bw,zsig_bw,Tmid,ztemp)
% mass balance function
% tinkered with by rsa 2/5/2011 
% and again 3/20/2011 by adding daily cycle
% and made into a function call 3/20/2011
% and simplified to just bw, bs, zela by rsa 2/21/2014 for when we arent
% worrying about thermal effects... and cleaned up by rsa 10/29/2015

% bs is set by a called ztemp and Tmid, Tmid being the MAT at the elevation
% ztemp
% bw is set by a called bwmax, and a zmid_bw, the elevation at which bw
% ramps to lower bw's over an elevation band prescribed by zsig_bw
% all other constants are set internally to this function, including the
% lapse rate, the pdd factors, the amplitudes of seasonal and daily T cycles

%%
% the necessary constants

L = 3.34e5; %latent heat of fusion of ice
c_snow = 2100; %heat capacity of ice...
rho_snow = 450;
rho_pc = 830; %pore closure density (see Pfeffer 1991 jgr)
pdd_star = 0.3; %governs how quickly pdd changes with snow depth
%bwmax = 1.0;
%zsig_bw = 500;
%zmid_bw = 2700;
%ztemp = 3700; %elevation at which we assign temps
% at D1 we have -3.7C now in modern times
%Tmid = -7.7; %mean annual temperature at ztemp
Tamp = 8; %amplitude of annual temp swing
Tamp_daily = 4; %daily amplitude
%lapse = -0.0065; % lapse rate in C/m a typical value
lapse = -0.0085;
%pdd_factor = 7e-3; %pdd factor in m melt/pdd
pdd_snow = 0.004; %m/day of melt per PDD
pdd_ice = 0.007;
    %note that hock 1999 uses 6.3 for snow and 4.4 for ice
    % which seems bass akwards...i expect faster melt for ice
    % and indeed in hock 2003 she reports in her big table 
    % values more like 4 and 6, or 6 and 8 for snow and ice, respectively
retention1 = 1.0; %fraction of melt-related heat retained in the snowpack
% where above the saturation limit
retention2 = 0.7; %fraction ... below the saturation limit
    
% dz = 10;
% zmin = 2500
% zmax = 4300;
% z = [zmin:dz:zmax];
imax = length(zmid);
   
dt = 1/10;
%dt = 1;
t = dt:dt:365;
P = 365;
Pday = 1;
jmax = length(t);


%%
% the calculations

% first the winter balance, a tanh function 
bw = (bwmax/2) + ((bwmax/2) * tanh((zmid-zmid_bw)/zsig_bw));
 
bs = zeros(size(bw));
% and second the summer balance
Tbar = Tmid + lapse*(zmid-ztemp);
for i = 1:imax; %loops through elevation
     T = Tbar(i) + Tamp*sin(2*pi*t/P) + Tamp_daily*sin(2*pi*t/Pday);
     T = max(0,T);
     depth = bw(i);
     sum_melt = 0;
     for j = 1:jmax; %loops thru time both day and year
%          if(depth==0)
%              pdd_factor = pdd_ice;
%          else
%              pdd_factor = pdd_snow;
%          end
% as of august 29th 2011 replacing this binary pdd call with an exponential
% model for increase in pdd as snow thins to mimic aging of snow
      pdd_factor = pdd_snow + (pdd_ice-pdd_snow).*exp(-depth./pdd_star);
      melt(j) = pdd_factor * T(j) * dt;
      sum_melt = sum_melt + melt(j);
      depth = depth - (melt(j));
      if(depth<=0)
          depth = 0;
      end
 
     end    
     bs(i) = -sum_melt; %negative because it is melt
end
 
b = bw+bs; %annual net balance
 
% now evolve the snowpack temperatures
T_snowfall = 0.6*(Tbar - Tamp); %assumed temp of snowfall
%here taken to be roughly mean daily temp in midwinter

cold_content = bw .* rho_snow * c_snow .* (-T_snowfall);
heat_melt = rho_snow * L .*(-bs);

% now deal with retention
retention = retention2 * ones(size(zmid));
unsaturated = find(-bs<((rho_pc - rho_snow)/rho_pc)*bw); 
% finds all sites where melt volume is insufficient to fill pores
% enough that they connect and allow runoff. after pfeffer 1991
retention(unsaturated) = retention1;

cold_remaining = cold_content - (retention .* heat_melt);
cold_remaining = max(0,cold_remaining);
 
new_temp = -cold_remaining ./(b .*rho_snow * c_snow);
Ts = new_temp; %surface temp b.c. for thermal code
 %just turn the cold energy back into a temperature of the 
 %net snowpack - i.e. b, not bw
 
 % now assess runoff
 % do not allow runoff if melt is insufficient to both remove cold content
 % AND fill pores. see pfeffer 1991 appendix
 runoff_condition = ((c_snow/L) * bw .* (-T_snowfall)) + (b .* ((rho_pc-rho_snow)/rho_snow));
 no_runoff = find(abs(bs) <= runoff_condition); 
 runoff = -bs;
 runoff(no_runoff) = 0;
 if(isempty(no_runoff)==1)
    runoff_line = max(zmid);
 else
    runoff_line = zmid(no_runoff(end));
 end
%ancillary calculations reported to the screen
%calculate mean db/dz over ablation zone:
ablation_zone = find(b<=0);
Ts(ablation_zone) = -1; % this is a cheat for now on ablation zone temps
% still needs to be fixed as of nov 17 2011
runoff_zone = find(-bs>0.7*bw); %from pfeffer 1991 jgr
%runoff_line = z(runoff_zone(1));
zela = zmid(ablation_zone(1));
dbdz = -b(1)/(zela-zmid(1));
