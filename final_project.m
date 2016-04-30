% simple 1d glacier model
% written june 12 2015 rsa and cole 
% modified by Sarah Crump and Simon Pendleton 11/16/15
% modified by AGT 12/15 for Cryosphere
% modified by AGT 4/16 for Modeling

% cleaned up code, loaded d18O data, taking out unnecessary code 
% topo data loaded and used for two cross sections

clear all
figure(1)
clf
%figure(2)
%clf
figure(3)
clf
figure(4)
clf
figure(5) 
clf


%% initialize

% material properties
rho_i = 917; %density of ice
g = 9.81; %gravity
%A = 2.1e-16; % Pa-3 yr-1
A = 0.01e-16; % do this to stiffen ice 11/18/15 flow law parameter
slide_ratio = 0.05;% ratio of sliding speed to internal defm speed for cold based glacier

%set up distance array
dx = 200;
%xmax = 20000; %meters for east/west cross section
xmax = 10000; %for north/south cross section
xedge = -xmax:dx:xmax;
x = -xmax+(dx/2):dx:xmax-(dx/2); % these are middles of cells
right = find(x>0); %right of center distance
left = find(x<=0); %left of center distance

%set up time array
dt = 0.02; %0.02 fo 2000 years
tmax = 2000; %years
%dt = 0.01;
%tmax = 150; %thousand years
t = 0:dt:tmax;


%% topo data

%load .mat file
load Renland_bed_topography.mat
load RECAP_lat.mat
load RECAP_lon.mat

y_elev = (topo(:,25))'; %cross section of bedrock for longitude, ' makes it a row
x_elev = topo(20,:); %cross section of bedrock for latitude
lat = RECAP_lat(:,25); 
lon = RECAP_lon(20,:);
y1 = linspace(-10000,10000,numel(lat));
x1 = linspace(-20000,20000,numel(lon));

%zb = interp1(x1,x_elev,x); %east/west cross section bedrock profile
zb = interp1(y1,y_elev,x); %north/south cross section bedrock profile

zbmin = min(zb); %min altitude of bedrock 
zbmax = max(zb); %max altitude of bedrock 
offcap = find(zb == zbmin); % finds all locations at base of pedestal
xcap = xmax;

%preallocate ice thickness
H = zeros(size(x)); %initial ice thickness
Hedge = zeros(size(xedge)); %initial thickness of ice at cell edges
z = zb+H; %initial ice height

% ice cap width as a function of distance is constant  
W0 = 3000; %ice cap width
W = ones(size(x)) * W0; %width as an array the size of x
Wedge = [W W(end)]; %define edges of W cells


%% temperature data 

Temp = zeros(size(t)); %equilibrium line altitude m

%load data file
load pleist_del18O_2000yrs.txt %data file has d180 every 1000 years for 2000 kyrs
age_data = pleist_del18O_2000yrs(:,1); %age column
d18O_data = pleist_del18O_2000yrs(:,2); %d18O column

%convert d18O to ela
age0 = age_data; %age in kyrs, multiply by 1000?
age=flipud(age0); %flip ages so it goes from past to present
rangeO = range(d18O_data); %range of isotope data 
scaled = -((d18O_data./rangeO) - min(d18O_data./rangeO)); %scales d18O data from 0-1
Temp0 = -15+(scaled.*(10)); %scales data to 10 C of temp change %-9.5, -15
Temp = interp1(age,Temp0,t); %interpolate temp for each age being calculated

%mass balance stuff 
tb = 5; 
bwmax= 0.35; %winter mass balance
zmid_bw= zbmax/2; %middle of ice5cap
zsig_bw= 700;
ztemp = zbmax+500; %temp at the top of the ice cap


%etc
imax = length(t);
nplots = 100;
tplot = tmax/nplots;
nframe = 0;
n=1;
tic


%% run

for i = 1:imax

    Tmid = Temp(i); %temp at time of iteration
    Hedge = H(1:end-1)+0.5*diff(H); %interpolates ice thickness to cell edges
    zedge = z(1:end-1) + 0.5*diff(z);
    S = (diff(z)/dx); % slope of ice surface calculated at cell edges
    Smid = diff(zedge)/dx;
   
    if rem(t(i),tb)==0; %every five iterations
        n=n+1;
        %call balance function to determine mass balance 
        bw = zeros(size(z)); 
        bs = zeros(size(z));
        [bw ,bs ,zela] = surface_balance_simple2_RSA(z,bwmax,zmid_bw,zsig_bw,Tmid,ztemp);

        b = bs + bw;
        ELA(n) = zela;
        tela(n) = t(i);
    end
    
    %calculate ice flux
    Udef = -sign(S).*(A/5).*((rho_i*g*abs(S)).^3).*(Hedge.^4); %mean deformation speed
    Q = Udef .* Hedge; %ice flux
    Qsl = (slide_ratio * Udef) .* Hedge; % ice discharge attributable to sliding...v crude here
    Q = Q + Qsl; %total ice flux across x
    Q =[0 Q 0]; %takes care of boundary conditions
    Qtotal = Q.*Wedge; %total ice flux in x and y directions
    
    taub = rho_i*g*abs(S).*Hedge; % magnitude of basal shear stress
    
    %calculate new ice thickness
    dHdt = b - (1./W).*(diff(Q.*Wedge)/dx); %continuity allowing width to vary
    H = H + (dHdt*dt); %updates ice thickness
    H = max(H,0); %keeps ice thickness from being negative
    H(offcap)=0; % prevents any glacier off the pedestal
    z = zb+H; %updates topography
    
    glacier = find(H>0); %where glacier exists
    if isempty(glacier)==1
        term(i) = 0;
        V_ice(i) = 0;
        Hmax(i) = 0;
    else
        term(i) = x(glacier(end)); % terminus position
        V_ice(i) = sum(H(glacier).*W(glacier))*dx; % glacier volume
        Hmax(i) = max(H);
    end

    % now for some plotting
    if rem(t(i),tplot)==0
        nframe = nframe + 1;

        figure(1) 
        subplot('position',[0.1 0.1 0.65 0.85]) %[left bottom width height]            
        plot(x/1000,zb,'k','linewidth',2)
        hold on
        plot(x/1000,z)
        plot(x/1000,zela*ones(size(x)),'g--','linewidth',2)
        axis([min(x)/1000 xmax/1000 zbmin-50 zbmax+1000])
        xlabel('Distance (km)','fontname','arial','fontsize',18)
        ylabel('Elevation (m)','fontname','arial','fontsize',18)
        set(gca,'fontsize',14,'fontname','arial')
        hold off

        subplot('position',[0.80 0.1 0.15 0.85])
        plot(b,z,'m','linewidth',2)
        hold on
        %plot(b0,zb,'k--','linewidth',2)
        plot(bs,z,'r','linewidth',1)
        plot(bw,z,'b','linewidth',1)
        plot(zeros(size(z)),z,'g--','linewidth',2)
        axis([min(bs) 1.2*max(bw) zbmin-50 zbmax+1000])
        xlabel('b (m/yr)','fontname','arial','fontsize',18)
        set(gca,'fontsize',14,'fontname','arial')
        % now stamp times
        time=num2str(t(i));
        timetext=strcat(time,' kyr');
        text(0.75*min(bs),zbmax+300,timetext,'fontsize',16)
        hold off
        %M(:,nframe) = getframe(gcf);

        %pause

        % now works only for positive x's...
        Qanal = cumsum(b(right).*W(right))*dx; %analytic solution for ss ice discharge in m3/yr
        Qanal = max(Qanal,0);

           
        %figure(2) 
        %plot(xedge/1000,Qtotal/1e6)
        %hold on
        %plot(x(right)/1000,Qanal/1e6,'g--')
        %axis([0 xmax/1000 0 (bcap/2)*(xmax/2)*mean(W)/1e6])
        %xlabel('Horizontal Distance (km)','fontname','arial','fontsize',18)
        %ylabel('Ice discharge (Mm^3/yr)','fontname','arial','fontsize',18)
        %set(gca,'fontsize',14,'fontname','arial')
        %hold off

        figure(4) 
        subplot(2,1,1)
        plot(x/1000,H)
        hold on
        xlabel('Distance (km)','fontname','arial','fontsize',18)
        ylabel('Ice thickness (m)','fontname','arial','fontsize',18)
        set(gca,'fontsize',14,'fontname','arial')
        subplot(2,1,2)
        plot(xedge(2:end-1)/1000,taub/1e5)
        hold on
        xlabel('Distance (km)','fontname','arial','fontsize',18)
        ylabel('Taub (bars)','fontname','arial','fontsize',18)
        set(gca,'fontsize',14,'fontname','arial')
        pause(0.05)
       % M(:,nframe) = getframe(gcf);

    end

end

toc


%% finalize

figure(3) 
%subplot(4,1,1)
% plot(tela,ELA,'r')
%     xlabel('Time (years)','fontname','arial','fontsize',18)
%     ylabel('ELA (m)','fontname','arial','fontsize',18)
%     set(gca,'fontsize',14,'fontname','arial')
subplot(3,1,1)
plot(t,Hmax,'r')
    xlabel('Time (kyr)','fontname','arial','fontsize',18)
    ylabel('H max (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
subplot(3,1,2)
plot(t,V_ice/1e9,'b')
    xlabel('Time (kyr)','fontname','arial','fontsize',18)
    ylabel('Volume (km^3)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
subplot(3,1,3)
plot(t,Temp,'b')
    xlabel('Time (kyr)','fontname','arial','fontsize',18)
    ylabel('Temp (C)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')


timescale = max(H)/max(b);

Hanalytic = max(H)*(1-(abs(x)/xcap)).^0.25;
outside = find(abs(x)>xcap);
Hanalytic(outside)=0;
figure(5)
plot(x/1000,H,'k')
hold on
plot(x/1000,Hanalytic,'g--','linewidth',2)
   xlabel('Distance (km)','fontname','arial','fontsize',18)
   ylabel('Ice thickness (m)','fontname','arial','fontsize',18)
   set(gca,'fontsize',14,'fontname','arial')
   
