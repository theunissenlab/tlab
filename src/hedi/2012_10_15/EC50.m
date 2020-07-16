function n=EC50(spec);
% return the freq for 50pc of the energy 
% spec slice is obtained via timefreq
r=sum(spec,2);
z=cumsum(r)/sum(r); 
n=length(z(z<0.5));
