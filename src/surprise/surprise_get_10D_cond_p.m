function cond_p = surprise_get_10D_cond_p(joint_p)
sz = size(joint_p);
cond_p = zeros(sz,'single');
if ndims(joint_p) ~=10
    error(strcat(['Function single_get_10D_cond_p called with ' num2str(ndims(joint_p)) ' argument.']));
end

for ii2 = 1:sz(2)
    for ii3 = 1:sz(3)
        for ii4 = 1:sz(4)
            for ii5 = 1:sz(5)
                for ii6 = 1:sz(6)
                    for ii7 = 1:sz(7)
                        for ii8 = 1:sz(8)
                            for ii9 = 1:sz(9)
                                for ii10 = 1:sz(10)
                                    total_p = single(sum(squeeze(joint_p(:,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10)))); % = P(D)
                                    if total_p > 0
                                        % Apply Bayes' rule: P(S|D) =
                                        % P(S,D) / P(D)
                                        cond_p(:,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10) = single(joint_p(:,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10) / total_p);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
