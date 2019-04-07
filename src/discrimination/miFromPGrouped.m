function mi = miFromPGrouped(p, ncat, mioff)
    % Returns mi for a confunsion matrix where the probability in blocked
    % diagonals are equal to p/n_i
    % and the off diagonals are (1-p)/(n - n_i)
    % ncat is the number in each group (or category)
    % mioff is subracted from the output to use this function to solve for p
    
    ng = length(ncat);    % Number of groups
    ns = sum(ncat);       % Number of stim
    
    mi = 0.0;     
    
    
    for ig=1:ng
        % The diagonal terms
        mi = mi + (ncat(ig)/ns)*p*( 2*log2(p/ncat(ig)) - log2(p/(ncat(ig)*ns)) );
        % The off-diagonal terms
        mi = mi + (ncat(ig)/ns)*(1-p)*( 2*log2((1-p)/(ns-ncat(ig))) - log2((1-p)/((ns-ncat(ig))*ns)) ); 
    end
    
    mi = mi - mioff;
end