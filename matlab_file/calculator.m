iter = (0:30);
theta = atan(0.5.^iter);
theta = theta * 180 / pi;
An = 1 / prod(cos(atan(0.5.^iter)))
FWL = 11;
const = round((1 / An) * (2 ^ FWL))
q_const = dec2bin(const)
%%
Nsigma_x = 0.50059;
SNR = 10*log10(1) - 10*log10(2 * Nsigma_x^2)

%%
SNR = 28;
Nsigma_x = sqrt((10^((10*log10(1) - SNR)/10)) / 2)

%% 
x_avg_radii = 2^0.5 + 10^0.5 + 10^0.5 + 18^0.5
n_avg_radii = (sigma^2 + sigma^2)^0.5
H = [1.3333 0 0.2712 0.6801; 0 1.3333 -0.6801 0.2712; -0.6318 -0.2621 -0.2118 -0.9433; 0.2621 -0.6318 0.9433 -0.2118];
norm(H(1,:))
norm(H(3,:))
y = H * [1; 1; 1; 3]

%% test
a = (1:3)'
b = [1:3; 2:4; 3:5]
c = a - b

c = c .* c
e = sum(c)
[m in] = min(e)

%% scater test
xtest = [1 0; 1 0];
xtest(:)
ytest = [1 2; 3 4];
scatter(xtest, ytest);

%%
scatter(1,2)

%% Normal distribution
format long
x = [-4, -3, -2, -1, 0];
mu = 0;
% sigma = 0.7071;
sigma = 1;
p = normcdf(x,mu,sigma) * 2

%%
format short
x = [-3 -1 1 3];
FWL = 11;
x_normal_q = round((x / sqrt(10)) * (2^FWL))

