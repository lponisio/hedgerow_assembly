function [VP, VA, A]=runObj2_PC(r_i)

global J_pattern

% Inputs

% ifrG_native: 0 & 4
% ifrG_Inv: 0 & 4
% PoA: 0 & 1
% nInv: 1
% k_level: 0 & 1
% opW: 0, 1 & 2
% febetaI: 1 & 4
% fepsilonI: 1 & 4
% ftauI: 1 & 2
% muAP: 1

ifrG_native=4;
ifrG_Inv=0;
PoA=0;
nInv=1;
k_level=1;
opW=0;
fbetaI=4;
fepsilonI=4;
ftauI=1;
muAP=2;
fyzeroI=0.1;
sem=0;

rand('seed',sem+r_i);

M=load(sprintf('m%04d.txt',r_i));% in each folder there are 800 matrices
[sA iA]=sort(sum(M));
[sP iP]=sort(sum(M,2));
M=M(iP,iA);


[In indx]=invasion(M,PoA,nInv,k_level,opW);
    
if PoA==0
    indxP=indx;indxA=[];
elseif PoA==1
    indxA=indx;indxP=[];
end

[rows cols]= size(In);
J_pattern = J_zero_pattern(In) ;

frG_native=ifrG_native*0.25;
frG_Inv=ifrG_Inv*0.25;
vectG_native=Gprim(frG_native,cols-length(indxA));
vectG_Inv=Gprim(frG_Inv,length(indxA));
vectG=[vectG_Inv vectG_native];
    
[VP, VA, A]=runVal_ext_cluster(indxP,indxA,vectG,In,muAP,fyzeroI,fbetaI,fepsilonI,ftauI);

%out_dir='res_800_N_M0/'; 
%save(sprintf('VP%d_vGinv%d_PoA%d_K%d_B%d_E%d_T%d_m%04d.mat',ifrG_native,ifrG_Inv,PoA,k_level,fbetaI,fepsilonI,ftauI,r_i),'VP')
%save(sprintf('VA%d_vGinv%d_PoA%d_K%d_B%d_E%d_T%d_m%04d.mat',ifrG_native,ifrG_Inv,PoA,k_level,fbetaI,fepsilonI,ftauI,r_i),'VA')
%save(sprintf('MV%d_vGinv%d_PoA%d_K%d_B%d_E%d_T%d_M%d%03d.mat',ifrG_native,ifrG_Inv,PoA,k_level,fbetaI,fepsilonI,ftauI,r_i),'MV')
%save(sprintf('Alfa%d_vGinv%d_PoA%d_K%d_B%d_E%d_T%d_m%04d.mat',ifrG_native,ifrG_Inv,PoA,k_level,fbetaI,fepsilonI,ftauI,r_i),'A')
    
end
