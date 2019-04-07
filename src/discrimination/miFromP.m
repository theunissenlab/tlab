function mi = miFromP(p, ns, mioff)
    % Returns mi for a confunsion matrix where the diagonals are equal to p
    % and the off-digonals are equal (1-p)/(ns-1)
    % ns is the size of the matrix
    % mioff is subracted from the output to use this function to solve for p
    
    mi = p*(2*log2(p)-log2(p/ns)) + (1-p)*(2*log2((1-p)/(ns-1))-log2((1-p)/((ns-1)*ns)));  
    mi = mi - mioff;
end