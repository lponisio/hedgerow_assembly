function vectG=Gprim(fG,npol)

% 19 de Sep 2010
% Determina qué animales son los forrajeadores adaptativos
% Sorteo al azar de una fracción fG de todos los pols

pol=randperm(npol);
nAF=round(fG*npol);
polAF=pol(1:nAF);
vectG=zeros(1,npol);
vectG(polAF)=1;

end
