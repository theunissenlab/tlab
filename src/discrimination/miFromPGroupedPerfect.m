function mi = miFromPGroupedPerfect(p, ncat, mioff)
    % Returns mi for a confunsion matrix where the diagonals are equal to p
    % and the off diagonals are 1-p/(number in each group - 1) and then
    % zero
    % ns is the size of the matrix
    % ncat is the number in each group (or category)
    % mioff is subracted from the output to use this function to solve for p
    
    ng = length(ncat);    % Number of groups
    ns = sum(ncat);       % Number of stim
    
    mi = p*(2*log2(p)-log2(p/ns));     % This comes from the diagonal
    
    % The off diagonal terms
    for ig=1:ng
        mi = mi + (ncat(ig)/ns)*(1-p)*( 2*log2((1-p)/(ncat(ig)-1)) - log2((1-p)/((ncat(ig)-1)*ns)) ); 
    end
    
    mi = mi - mioff;
end