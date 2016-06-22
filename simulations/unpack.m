function [plants nectar animals alphas] = unpack(y, network_metadata)

n = network_metadata.plant_qty ;
m = network_metadata.animal_qty ;

plants=y(1:n,1) ;
nectar=y(n+1:2*n,1);
animals=y(2*n+1: 2*n+m,1);

alphas=sparse( n, m ) ;
alphas(network_metadata.nz_pos) = y(2*n+m+1:end,1) ;

end
