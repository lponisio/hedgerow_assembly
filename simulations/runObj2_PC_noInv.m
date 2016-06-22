function [VP, VA, MV, A]=runObj2_PC_noInv(r_i)

global J_pattern

% Inputs

muAP=1;
ifrG_native=4;
sem=0;

rand('seed',sem+r_i);

names=dir('*.txt');
nfiles=length(names);

for i=1:nfiles
    n=names(i).name;
    In=load(n);
    [rows cols]= size(In);
    
    J_pattern = J_zero_pattern(In) ;
    frG_native=ifrG_native*0.25;
    vectG=Gprim(frG_native,cols);
    
    [VP, VA, MV, A]=runVal_ext_cluster_noInv(vectG,In,muAP);
    VP=full(VP);
    VA=full(VA);
    MV=full(MV);
    A=full(A);
    
    n=n(1:end-4);
    
    dlmwrite(sprintf('VP_%s.cvs',n),VP)
    dlmwrite(sprintf('VA_%s.cvs',n),VA)
    dlmwrite(sprintf('matV_%s.cvs',n),MV)
    dlmwrite(sprintf('matA_%s.cvs',n),A)

end

end