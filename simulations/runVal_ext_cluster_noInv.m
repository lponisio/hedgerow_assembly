%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified 3/31/16: Invasions part removed to check a possible bug in the
% code given that the trajectories of the abundances change after the
% dynamics is stoped and run again even when no new species is introduced to the network 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VP, VA, MV, A]=runVal_ext_cluster_noInv(vectG,In,muAP)

global network_metadata J_pattern

tmax=3000; %tiempo maximo de simulacion

[m n]=size(In);
B=sparse(In);

% Parametros de la distribucion uniforme
varp=1e-1;
vara=1e-4;
%vara=0;
mC=0.2; vC=vara;
mE=0.8; vE=varp;
mb=0.4; vb=vara;
mU=0.06; vU=varp;
mw=1.2; vw=varp;
mB=0.2; vB=varp;
mG=2; vG=vara;
mg=0.4; vg=varp;
mphi=0.04; vphi=varp;
mtau=1; vtau=vara;
mepsilon=1; vepsilon=varp;
vmA=vara; vmP=varp;

if muAP==1
    mmA=0.01; mmP=0.002; % mueren animales s/AF
elseif muAP==2
    mmA=0.001; mmP=0.015; % mueren plantas s/AF (y con AF tb)
elseif muAP==3
    mmA=0.001; mmP=0.001; % todos viven s/AF
elseif muAP==4
    mmA=0.03; mmP=0.005; % mueren plantas & animales s/AF
end

% Valores de los pars del modelo obtenidos de una distr uniforme
% (10%meanP)-meanP+(10%meanP); (0.01%meanA)-meanA+(0.01%meanA)
c=uniform_rand(mC,vC,m,n).*B;
e=uniform_rand(mE,vE,m,n).*B;
b=uniform_rand(mb,vb,m,n).*B;

u=uniform_rand(mU,vU,m,1);
Beta=uniform_rand(mB,vB,m,1);

G=uniform_rand(mG,vG,n,1).*vectG';
g=uniform_rand(mg,vg,m,1);
mu_a=uniform_rand(mmA,vmA,n,1);
mu_p=uniform_rand(mmP,vmP,m,1);
w=uniform_rand(mw,vw,m,1);
%w=distLin(sum(In,2));
phi=uniform_rand(mphi,vphi,m,1);

tau=uniform_rand(mtau,vtau,n,1);

epsilon=uniform_rand(mepsilon,vepsilon,m,1);


%Magic thing. Do not change the following line
network_metadata = create_metadata(B, e, mu_p, mu_a, c, b, u, w, Beta, G, g, phi, tau, epsilon) ;

%Give initial state
mz=0.5; vz=0;%vz=1e-1;

yzero=uniform_rand(mz,vz,2*m+n,1);

initial_plants=yzero(1:m);

initial_nectar=yzero(m+1:2*m);

initial_animals=yzero(2*m+1:2*m+n);

initial_alphas=B;

%initial_alphas(:,indxA)=0;

%Normalization an packing. You should not need to touch these two lines
initial_alphas=initial_alphas*diag(sum(initial_alphas).^(-1));
initial_alphas=initial_alphas(network_metadata.nz_pos) ;

%Integrate the damn thing
initial_state=full([initial_plants;initial_nectar;initial_animals;initial_alphas]);
tspan = [0 tmax];

%polinizacion_rhs(0,initial_state)
options = odeset('JPattern', J_pattern,'NonNegative',1:2*m+n) ;
[t y]=ode15s(@polinizacion_rhs,tspan,initial_state, options) ;

%tf=length(t);
yf = y(end,:)';

%% Comienzan las INVASIONES despues de llegar al equilibrio
[plantsf nectarf animalsf alphasf] = unpack(yf, network_metadata);

% Registro pre-invasion

tmp = alphasf * diag(sparse(animalsf.*tau)) ;
%sVi = sum(tmp,2); % total de visitas que recibe cada planta (per-capita)
sigma = diag(sparse(plantsf.*epsilon)) * alphasf ; %sigma nxm sparse
sigma = sigma * diag(sparse(1./(sum(sigma)+realmin))) ;
pol_event= sum(sigma .* tmp, 2);% suma per-capita
%meansigma= mean(sigma,2);

tmp = alphasf * diag(sparse(tau)) ;
%sj = sum(  diag(sparse(plantsf)) * tmp  )';% cant de visitas de un animal (per-capita)

tmp = (diag(sparse(nectarf)) * tmp) .* b;
N_extract = sum(tmp)'; % sum of the resources that each individual extracts

EUp=2e-2;
EUa=1e-3;
%
vzp=plantsf<EUp;
vza=animalsf<EUa;
%% Resultados
A=alphasf;
MV=diag(plantsf)*alphasf*diag(animalsf.*tau);
VP=[vzp plantsf nectarf pol_event];% sparse npl*4
VA=[vza animalsf N_extract];% sparse npol*3



% Graficando

figure
subplot (3,1,1)
plot(t,y(:,1:m))
subplot (3,1,2)
plot(t,y(:,m+1:2*m))
subplot (3,1,3)
plot(t,y(:,2*m+1:2*m+n))

end
