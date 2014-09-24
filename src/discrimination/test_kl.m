       
% Testing the formulat for the KL divergence

mean_other = 40;
std_other = 5:5:20;
std_self = 5:5:20;
n_stds = length(std_self);


% Effect of differences in mean

mean_self = 5:40;
ntest = length(mean_self);
kl_divergence = zeros(ntest, n_stds);

for in=1:n_stds
    for it=1:ntest
        kl_divergence(it,in) = log(std_other(in)./std_self(in)) + ...
            (mean_other - mean_self(it))^2./(2*std_other(in).^2) + ...
            (std_self(in).^2 - std_other(in).^2)./(2*std_other(in).^2);
        kl_divergence(it,in) = kl_divergence(it,in)./log(2);
    end
end

figure(1);
plot(mean_self, kl_divergence);
xlabel('Mean of Self (40 for Other)');
ylabel('KL Divergence (bits)');
legend('5', '10', '15', '20');
title('Std other = Std self (See Legend for value)');

clear all;
mean_self = 20;
mean_other = 40;
std_other = 20;
std_self = 5:40;

ntest = length(std_self);
kl_divergence = zeros(ntest, 1);

for it=1:ntest
    kl_divergence(it) = log(std_other./std_self(it)) + ...
        (mean_other - mean_self)^2./(2*std_other.^2) + ...
        (std_self(it).^2 - std_other.^2)./(2*std_other.^2);
    kl_divergence(it) = kl_divergence(it)./log(2);
end

figure(2);
plot(std_self, kl_divergence);
xlabel('Std of Self (20 for Other)');
ylabel('KL Divergence (bits)');
title('Mean other = 40 Mean self = 20');


clear all;
mean_self = 20;
mean_other = 40;
std_self = 20;
std_other = 5:40;

ntest = length(std_other);
kl_divergence = zeros(ntest, 1);

for it=1:ntest
    kl_divergence(it) = log(std_other(it)./std_self) + ...
        (mean_other - mean_self)^2./(2*std_other(it).^2) + ...
        (std_self.^2 - std_other(it).^2)./(2*std_other(it).^2);
    kl_divergence(it) = kl_divergence(it)./log(2);
end

figure(3);
plot(std_other, kl_divergence);
xlabel('Std of Other (20 for Self)');
ylabel('KL Divergence (bits)');
title('Mean other = 40 Mean self = 20');




