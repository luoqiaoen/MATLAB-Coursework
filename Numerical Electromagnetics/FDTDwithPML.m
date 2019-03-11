clear all;
clc;
%% Choose if PML is used
% if 1, PML condition is used, if 0, no PML layers. 
PML = 1;

%% Arbitrary point location
xx = 100;
yy = 50;
%% Testing Parameters
% The user can give a frequency of his choice for sinusoidal waves in Hz 
frequency=1.5e+13;

% Grid numbers in x (xdim) and y (ydim) directions
xdim=300;
ydim=300;

%Total no of time steps
time_tot=500;

%Position of the source (center of the domain)
xsource=150;
ysource=150;

%Courant stability factor
S=1/(2^0.5);

% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;

% Spatial grid step length (spatial grid step= 1 micron and can be changed)
delta=1e-6;
% Temporal grid step obtained using Courant condition
deltat=S*delta/c;

% Initialization of field matrices
Ez=zeros(xdim,ydim);
Ezx=zeros(xdim,ydim);
Ezy=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);

% Initialization of permittivity and permeability matrices
epsilon=epsilon0*ones(xdim,ydim);
mu=mu0*ones(xdim,ydim);

% Initializing electric conductivity matrices in x and y directions
sigmax=zeros(xdim,ydim);
sigmay=zeros(xdim,ydim);
sigma_starx =zeros(xdim,ydim);
sigma_stary =zeros(xdim,ydim);

if PML == 1
    %Boundary width of PML in all directions
    bound_width=25;
    
    %Order of polynomial on which sigma is modeled
    gradingorder=6;
    
    %Required reflection co-efficient
    refl_coeff=1e-8;
    
    %Polynomial model for sigma
    sigmamax=(-log10(refl_coeff)*(gradingorder+1)*epsilon0*c)/(2*bound_width*delta);
    boundfact1=((epsilon(xdim/2,bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact2=((epsilon(xdim/2,ydim-bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact3=((epsilon(bound_width,ydim/2)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact4=((epsilon(xdim-bound_width,ydim/2)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    x=0:1:bound_width;
    for i=1:1:xdim
        sigmax(i,bound_width+1:-1:1)=boundfact1*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
        sigmax(i,ydim-bound_width:1:ydim)=boundfact2*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
    end
    for i=1:1:ydim
        sigmay(bound_width+1:-1:1,i)=boundfact3*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1))';
        sigmay(xdim-bound_width:1:xdim,i)=boundfact4*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1))';
    end
    
    %Magnetic conductivity matrix obtained by Perfectly Matched Layer condition
    %This is also split into x and y directions in Berenger's model
    sigma_starx=(sigmax.*mu)./epsilon;
    sigma_stary=(sigmay.*mu)./epsilon;
end

%Multiplication factor matrices for H matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed
C=((mu-0.5*deltat*sigma_starx)./(mu+0.5*deltat*sigma_starx)); 
D=(deltat/delta)./(mu+0.5*deltat*sigma_starx);
A=((mu-0.5*deltat*sigma_stary)./(mu+0.5*deltat*sigma_stary)); 
B=(deltat/delta)./(mu+0.5*deltat*sigma_stary);
                          
%Multiplication factor matrices for E matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                          
E=((epsilon-0.5*deltat*sigmax)./(epsilon+0.5*deltat*sigmax)); 
F=(deltat/delta)./(epsilon+0.5*deltat*sigmax);   
G=((epsilon-0.5*deltat*sigmay)./(epsilon+0.5*deltat*sigmay)); 
H=(deltat/delta)./(epsilon+0.5*deltat*sigmay);

%% record arbitrary point in space over time
arb = zeros(1,time_tot);


%% Update loop begins
for n=1:1:time_tot
    % Setting time dependent boundaries to update only relevant parts of the 
    % matrix where the wave has reached to avoid unnecessary updates.
    if n<xsource-2
        n1=xsource-n-1;
    else
        n1=1;
    end
    if n<xdim-1-xsource
        n2=xsource+n;
    else
        n2=xdim-1;
    end
    if n<ysource-2
        n11=ysource-n-1;
    else
        n11=1;
    end
    if n<ydim-1-ysource
        n21=ysource+n;
    else
        n21=ydim-1;
    end
    
    %matrix update instead of for-loop for Hy and Hx fields
    Hy(n1:n2,n11:n21)=A(n1:n2,n11:n21).*Hy(n1:n2,n11:n21)+B(n1:n2,n11:n21).*(Ezx(n1+1:n2+1,n11:n21)-Ezx(n1:n2,n11:n21)+Ezy(n1+1:n2+1,n11:n21)-Ezy(n1:n2,n11:n21));
    Hx(n1:n2,n11:n21)=C(n1:n2,n11:n21).*Hx(n1:n2,n11:n21)-D(n1:n2,n11:n21).*(Ezx(n1:n2,n11+1:n21+1)-Ezx(n1:n2,n11:n21)+Ezy(n1:n2,n11+1:n21+1)-Ezy(n1:n2,n11:n21));
    
    %matrix update instead of for-loop for Ez field
    Ezx(n1+1:n2+1,n11+1:n21+1)=E(n1+1:n2+1,n11+1:n21+1).*Ezx(n1+1:n2+1,n11+1:n21+1)+F(n1+1:n2+1,n11+1:n21+1).*(-Hx(n1+1:n2+1,n11+1:n21+1)+Hx(n1+1:n2+1,n11:n21));
    Ezy(n1+1:n2+1,n11+1:n21+1)=G(n1+1:n2+1,n11+1:n21+1).*Ezy(n1+1:n2+1,n11+1:n21+1)+H(n1+1:n2+1,n11+1:n21+1).*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1));
    

    tstart=1;
    N_lambda=c/(frequency*delta);
    Ezx(xsource,ysource)=0.5*sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
    Ezy(xsource,ysource)=0.5*sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));

   
    
    Ez=Ezx+Ezy;
    
    %Movie type colour scaled image plot of Ez
    imagesc(delta*1e+6*(1:1:xdim),(delta*1e+6*(1:1:ydim))',Ez',[-1,1]);colorbar;
    title(['E_z in a spatial domain with PML boundary and at time = ',num2str(round(n*deltat*1e+15)),' fs']); 
    xlabel('x (in \mu m)');
    ylabel('y (in \mu m)');
    axis image;
    drawnow;
    arb(n) = Ez(xx,yy);
end

figure(2)
plot(0:n*deltat*1e+15/(time_tot-1):n*deltat*1e+15,arb);
xlabel('time [fs]');
ylabel('E_z in selected fix point');
