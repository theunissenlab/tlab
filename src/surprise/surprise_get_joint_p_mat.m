function [joint_p,which_cubes,what_position] = surprise_get_joint_p_mat(one_stim_one_bin_PCA_Tmat,PC_Maxes,PC_Mins,PC_Nbins,big_familiarity)
%  This function finds the conditional probability matrix for predicting
%  one frequency band, using length(PC_Maxes)-1 of the total number of PCs,
%  with maxes for each PC PC_Maxes, mins for each PC PC_Mins, and the
%  number of bins for each PC PC_Nbins.  Here the PC vectors are all
%  actually the vectors of the target plus the PC outputs in increasing
%  eigenvalue strength.
%  NB PC_Maxes and PC_Mins MUST be at least as wide as the PC data;
%  otherwise this program will crash.  (It's possible to fix this problem
%  by making this routine slower, but I won't give up the speed.  Full
%  power, Mr. Scott!)

if any([length(PC_Mins) ~= length(PC_Maxes) length(PC_Nbins) ~= length(PC_Maxes)])
    error(strcat(['Error in function surprise_get_joint_p_mat: passed maxes, mins and nbins with different number of dimensions.' char(10) ...
        'Specifically, length(PC_Maxes) = ' num2str(length(PC_Maxes)) ' length(PC_Mins) = ' num2str(length(PC_Mins)) 'length(PC_Nbins) = ' num2str(length(PC_Nbins)) '.']));
end
%I'm replacing the line below (originally used in Gill et al. 2008) with
%the line after it because I don't know how other people's spectrograms are
%going to be scaled.
%epsilon = 1e-6;  %Make the maxes and the mins just a bit biger to handle the equality cases.
epsilon = 1e-5 * (max(PC_Maxes) - min(PC_Mins));

PC_Maxes = PC_Maxes + epsilon;
PC_Mins = PC_Mins - epsilon;

nD = size(one_stim_one_bin_PCA_Tmat,1);

to_take = [1 (nD - length(PC_Maxes) + 2):nD];

working = one_stim_one_bin_PCA_Tmat(to_take,:);
if size(PC_Mins,1) == 1
    PC_Mins = PC_Mins';
end

if size(PC_Maxes,1) == 1
    PC_Maxes = PC_Maxes';
end

if size(PC_Nbins,1) == 1
    PC_Nbins = PC_Nbins';
end
joint_p = zeros(PC_Nbins','single');

%  Scale and stretch working so that the integer parts of working
%  correspond to the bin number they belong in minus 1

working = working - PC_Mins*ones(1,size(working,2));
working = working ./ ((PC_Maxes-PC_Mins)*ones(1,size(working,2)));
working = working .* ((PC_Nbins-1)*ones(1,size(working,2)));

which_cubes = ceil(working);  %  What are the coordinates of the lower corner of the hypercube of the joint probability estimate to increment?
what_position = working+1-which_cubes;  %  Where in the hypercube is the center of probability?
square_from_small = what_position .^2;
square_from_large = (1-what_position).^2;
switch length(PC_Nbins)
    %  In the released package, only 10-D calculations (9 PCs in D and 1 of
    %  stim) are supported for now.  Feel free to edit out a bunch of loops
    %  and colons to get analogous functions with fewer dimensions; email
    %  me at patgill@berkeley.edu if you'd like some help.
    case 2
        for t = 1:size(working,2)
            cube = which_cubes(:,t);
            joint_p(cube(1) + (0:1), cube(2) + (0:1)) ...
                = joint_p(cube(1) + (0:1), cube(2) + (0:1)) ...
                + single(big_familiarity(t)*get_2D_gaus_proj(square_from_small(:,t),square_from_large(:,t)));
        end
    case 3
        for t = 1:size(working,2)
            cube = which_cubes(:,t);
            joint_p(cube(1) + (0:1), cube(2) + (0:1),cube(3) + (0:1)) ...
                = joint_p(cube(1) + (0:1), cube(2) + (0:1),cube(3) + (0:1)) ...
                + single(big_familiarity(t)*get_3D_gaus_proj(square_from_small(:,t),square_from_large(:,t)));
        end
    case 7
        for t = 1:size(working,2)
            cube = which_cubes(:,t);
            joint_p(cube(1) + (0:1), cube(2) + (0:1),cube(3) + (0:1),cube(4) + (0:1),cube(5) + (0:1),cube(6) + (0:1),cube(7) + (0:1)) ...
                = joint_p(cube(1) + (0:1), cube(2) + (0:1),cube(3) + (0:1),cube(4) + (0:1),cube(5) + (0:1),cube(6) + (0:1),cube(7) + (0:1)) ...
                + single(big_familiarity(t)*get_7D_gaus_proj(square_from_small(:,t),square_from_large(:,t)));
        end
    case 10
        for t = 1:size(working,2)
            cube = which_cubes(:,t);
            joint_p(cube(1) + (0:1), cube(2) + (0:1),cube(3) + (0:1),cube(4) + (0:1),cube(5) + (0:1),cube(6) + (0:1),cube(7) + (0:1),cube(8) + (0:1),cube(9) + (0:1),cube(10) + (0:1)) ...
                = joint_p(cube(1) + (0:1), cube(2) + (0:1),cube(3) + (0:1),cube(4) + (0:1),cube(5) + (0:1),cube(6) + (0:1),cube(7) + (0:1),cube(8) + (0:1),cube(9) + (0:1),cube(10) + (0:1))...
                + single(big_familiarity(t)*get_10D_gaus_proj(square_from_small(:,t),square_from_large(:,t)));
        end
    otherwise
        error(strcat(['Error in function surprise_get_joint_p_mat:  Arrays of dimension ' num2str(length(PC_Nbins)) ' are not yet supported.  Do some more coding.']));
end

