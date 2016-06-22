%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified 8/06/15 to add invasion of plants in networks wihtout AF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VP VA A]=runVal_ext_cluster(indxP,indxA,vectG,In,muAP,fyzeroI,fbetaI,fepsilonI,ftauI)

global network_metadata indRemP indRemA J_pattern

tmax=3000; %tiempo maximo de simulacion

[m n]=size(In);
B=sparse(In);
VP=cell(1,2);
VA=cell(1,2);
MV=cell(1,2);
A=cell(1,2);

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
Beta(indxP)=fbetaI*Beta(indxP);

G=uniform_rand(mG,vG,n,1).*vectG';
g=uniform_rand(mg,vg,m,1);
mu_a=uniform_rand(mmA,vmA,n,1);
mu_p=uniform_rand(mmP,vmP,m,1);
w=uniform_rand(mw,vw,m,1);
%w=distLin(sum(In,2));
phi=uniform_rand(mphi,vphi,m,1);

tau=uniform_rand(mtau,vtau,n,1);
tau(indxA)=ftauI*tau(indxA);
epsilon=uniform_rand(mepsilon,vepsilon,m,1);
epsilon(indxP)=fepsilonI*epsilon(indxP);

%Magic thing. Do not change the following line
network_metadata = create_metadata(B, e, mu_p, mu_a, c, b, u, w, Beta, G, g, phi, tau, epsilon) ;

%Give initial state
mz=0.5; vz=0;%vz=1e-1;

yzero=uniform_rand(mz,vz,2*m+n,1);

initial_plants=yzero(1:m);
initial_plants(indxP)=0;

initial_nectar=yzero(m+1:2*m);
initial_nectar(indxP)=0;

initial_animals=yzero(2*m+1:2*m+n);
initial_animals(indxA)=0;
initial_alphas=B;
initial_alphas(indxP,:)=0;
%initial_alphas(:,indxA)=0;

%Normalization an packing. You should not need to touch these two lines
initial_alphas=initial_alphas*diag(sum(initial_alphas).^(-1));
initial_alphas=initial_alphas(network_metadata.nz_pos) ;

%Integrate the damn thing
initial_state=full([initial_plants;initial_nectar;initial_animals;initial_alphas]);
tspan = [0 tmax];

indRemP=indxP;
indRemA=indxA;

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
A{1}=alphasf;
MV{1}=diag(plantsf)*alphasf*diag(animalsf.*tau);
VP{1}=[vzp plantsf nectarf pol_event];% sparse npl*4
VA{1}=[vza animalsf N_extract];% sparse npol*3

assert(all(abs(plantsf(indxP))<1e-4));
assert(all(abs(nectarf(indxP))<1e-4));
assert(all(abs(animalsf(indxA))<1e-4));
assert(all(all(alphasf>-1e-4 & alphasf<1.0001)));
assert( all( abs(sum(alphasf) - 1.0)<1e-3 ) ) ;

%% Invasion

indxP2=find(vzp(length(indxP)+1:end));
indxA2=find(vza(length(indxA)+1:end));

indRemP=indxP2;
indRemA=indxA2;

initial_plants=plantsf;
initial_plants(indxP)=fyzeroI*mean(plantsf);
initial_plants(indxP2)=0;

initial_nectar=nectarf;
initial_nectar(indxP)=fyzeroI*mean(nectarf);
initial_nectar(indxP2)=0;

initial_animals=animalsf;
initial_animals(indxA)=fyzeroI*mean(animalsf);
initial_animals(indxA2)=0;

% Introduciendo a la planta en la matriz Alfa

if sum(indxP)>0
    
    if sum(vectG)==0
        
        initial_alphas=B; % Getting back the whole adjacency matrix (including the interactions of the introduced plant that where removed above)
        initial_alphas=initial_alphas*diag(sum(initial_alphas).^(-1));
        initial_alphas=initial_alphas(network_metadata.nz_pos) ;
    
    else
        
        a_focal=alphasf( :,logical(B(indxP,:)) );
        f_max=a_focal==ones(size(a_focal))*diag(max(a_focal));
        a_focal(f_max)=a_focal(f_max)-0.0001;
        a_focal(1,:)=0.0001;
        alphasf( :,logical(B(indxP,:)) )=a_focal;
        initial_alphas=alphasf(network_metadata.nz_pos) ;
        
    end
    
end


%Integrate the damn thing
initial_state=full([initial_plants;initial_nectar;initial_animals;initial_alphas]);
%tspan = linspace(0,3000,3001);


%polinizacion_rhs(0,initial_state)
[t2 y2]=ode15s(@polinizacion_rhs,tspan,initial_state,options) ;

%tf2=length(t2);
yf2 = y2(end,:)';
[p N a Alpha] = unpack(yf2, network_metadata ) ;
%p(indxP)
%a(indxA)

assert(all(all(Alpha>-1e-4 & Alpha<1.0001)));
assert( all( abs(sum(Alpha) - 1.0)<1e-3 ) ) ;

%% Calculo de variables para analizar los mecanismos explicativos

% Plants
tmp = Alpha * diag(sparse(a.*tau)) ;
sVi = sum(tmp,2); % total de visitas que recibe cada planta (per-capita)
sigma = diag(sparse(p.*epsilon)) * Alpha ; %sigma nxm sparse
sigma = sigma * diag(sparse(1./(sum(sigma)+realmin))) ;
%meansigma= mean(sigma,2);
pol_event= sum(sigma .* tmp, 2);% suma per-capita
%Gamma = g .* (1 - u'*p - w.*p + u.*p) ;% non sparse
%NpVi=N./(sVi+realmin);

% Animals
tmp = Alpha * diag(sparse(tau)) ;
%sj = sum(diag(sparse(p)) * tmp)';
tmp = (diag(sparse(N)) * tmp) .* b;
N_extract = sum(tmp)'; % sum of the resources that each individual extracts

% Graficando condicion post-extincion

Y=[y;y2];
T=[t;t2+tmax];

figure
subplot (3,1,1)
plot(T,Y(:,1:m))
subplot (3,1,2)
plot(T,Y(:,m+1:2*m))
subplot (3,1,3)
plot(T,Y(:,2*m+1:2*m+n))

%% Number of plant and animal extincts (below EU)
%
vzp=p<EUp;
vza=a<EUa;

%% Resultados

A{2}=Alpha;
MV{2}=diag(p)*Alpha*diag(a.*tau);
VP{2}=[vzp p N pol_event];% sparse npl*4
VA{2}=[vza a N_extract];

end
