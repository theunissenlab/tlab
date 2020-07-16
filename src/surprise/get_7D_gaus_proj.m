function out = get_7D_gaus_proj(square_from_small,square_from_large);
%If thee's a Gaussian centered at coordinates giving square_from_small,
%square_from_large, this function will give its values at the verticies of
%the hypercube.
diffs =  square_from_large - square_from_small;
out = zeros(2,2,2,2,2,2,2);
for ii1 = 1:2
    for ii2 = 1:2
        for ii3 = 1:2
            for ii4 =1:2
                for ii5 = 1:2
                    for ii6 = 1:2
                        for ii7 = 1:2

                            out(ii1,ii2,ii3,ii4,ii5,ii6,ii7) = exp(-sum(square_from_small + diffs.*([ii1; ii2; ii3; ii4; ii5; ii6; ii7]-1)));

                        end
                    end
                end
            end
        end
    end
end

nconst = sum(sum(sum(sum(sum(sum(sum(out)))))));
out = out / nconst;
