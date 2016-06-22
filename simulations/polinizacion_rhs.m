%ODE's RHS. 
% New version. Tries to avoid 0/0 cases to simplify simulations
% where species are killed off or are to be added in the middle of the run

function dx = polinizacion_rhs(t,x)

global network_metadata indRemP indRemA

plant_qty  = network_metadata.plant_qty ;
animal_qty = network_metadata.animal_qty ;
nz_pos  = network_metadata.nz_pos ;
e       = network_metadata.e ;
mu_p    = network_metadata.mu_p ;
mu_a    = network_metadata.mu_a ;
c       = network_metadata.c ;
b       = network_metadata.b ;
u       = network_metadata.u ;
w       = network_metadata.w ;
Beta    = network_metadata.Beta ;
G       = network_metadata.G ;
g       = network_metadata.g ;
phi     = network_metadata.phi ;
tau     = network_metadata.tau ;
epsilon = network_metadata.epsilon ;
In      = network_metadata.In ;

[p N a Alpha] = unpack(x, network_metadata ) ;
p(indRemP)=0;% forzamos elementos a cero para evitar poblemas de resucitacion con el integrador stiff
N(indRemP)=0;
a(indRemA)=0;

% Model's specific computation begins here

sigma = diag(sparse(p.*epsilon)) * Alpha ; %sigma nxm sparse
sigma = sigma * diag(sparse(1./(sum(sigma)+realmin))) ;

Gamma = g .* (1 - u'*p - w.*p + u.*p) ;% non sparse

tmp = Alpha * diag(sparse(a.*tau)) ; %tmp nxm sparse

dp = ( (Gamma .* sum(e .* sigma .* tmp, 2)) - mu_p) .* p;
tmp = (diag(sparse(N)) * tmp) .* b ;
da = sum(c .* tmp, 1)' - mu_a .* a ;
dN = Beta .* p - phi.*N - sum(tmp, 2) ;

% Adaptation magic starts here

%Fernanda's fitness function
DH = diag(sparse(N)) * sparse(c.*b) ; %nxm sparse
DH(Alpha<0)=-DH(Alpha<0) ;

wavg = sum(Alpha.*DH) ; %Weights for average. nxm sparse

%This is the replicator equation
dAlpha = Alpha.*DH - Alpha*diag(sparse(wavg)) ;
dAlpha = dAlpha*diag(sparse(G)) ;


% Now pack the answer
dx = full([dp; dN; da; dAlpha(nz_pos)]) ;
