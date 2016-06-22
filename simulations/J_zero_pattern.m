%Last mod Oct 7 2011. 
% p_alpha fixed
% alpha_N Still not ready

function Jz=J_zero_pattern(In)

[m n]=size(In);
nz_pos=find(In); %indices in linear space
alpha_qty=length(nz_pos);

%preallocate submatrices that build the answer
alpha_alpha=sparse(m,n);
a_alpha=sparse(n,alpha_qty);
p_alpha=sparse(m,alpha_qty);

% Integer division trick
% nz_div(k)==i iff nz_pos(k) corresponds to an alpha of animal k
nz_idiv=idivide(uint32(nz_pos)-1,uint32(m))+1;

for j=1:n
    alpha_indx_aj=sparse(find(nz_idiv==j));
    a_alpha(j,alpha_indx_aj)=1;
    alpha_alpha(alpha_indx_aj,alpha_indx_aj)=1;
end

Tmatrix=(In*1) ;
Tmatrix(nz_pos)=1:length(nz_pos);

for i=1:m
   idx = find(In(i,:)==0) ;
   tmp=In ;
   tmp(:,idx)=0 ;
   p_alpha(i,Tmatrix(find(tmp))) = 1 ;
end

alpha_N=p_alpha';
N_alpha=p_alpha ;

Jz=[sparse(ones(m,m))   sparse(m,m)       In                  p_alpha;
    sparse(1:m,1:m,1)   sparse(1:m,1:m,1) In                  N_alpha;
    sparse(n,m)         In'               sparse(1:n,1:n,1)   a_alpha;
    sparse(alpha_qty,m) alpha_N           sparse(alpha_qty,n) alpha_alpha];

end
