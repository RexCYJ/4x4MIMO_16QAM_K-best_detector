%% Decomposition Script
%   date: 2022-04-25

%%
format short

%% clear all data
clear

%% Read CSV file
% ------replace the path 
path = 'D:\Rex\NTHU\Implementation\cpp_simul\matrix.csv';   
data = csvread(path);

%% SNR 
Nsigma_x = 0.0707;
SNR = 10*log10(1) - 10*log10(2 * Nsigma_x^2);

%% analyze information
% CSV formation:
%       H: column1, column2, ..., columnN
%       R: column1, column2, ..., columnN
%    ┌ -----------------------------------------
%    | H1(1,1), H1(2,1), H1(3,1), ..., H1(N,N)
%    | R1(1,1), R1(2,1), R1(3,1), ..., R1(N,N)
%    | Y
%    | QhY
%    | X
%    | H2(1,1), H2(2,1), H2(3,1), ..., H2(N,N)
%    | R2(1,1), R2(2,1), R2(3,1), ..., R2(N,N)
%    | ...
%    └ -----------------------------------------

[Dsize, esize] = size(data);    % number of data, data length
dline = 5;
dsize = Dsize / dline;
H           = data((1:dline:Dsize), :);     % get data of H
R           = data((2:dline:Dsize), :);     % get data of R
N           = sqrt(length(H(1,:)));         % derive matrix size, assume M = N
Hsize       = [N N];                        % record matrix size
Y           = data((3:dline:Dsize), (1:N));
Y_h_cpp     = data((4:dline:Dsize), (1:N));
X           = data((5:dline:Dsize), (1:N));
QAM_16_sol_x = (-3:2:3);

%% RUN ALL MATRICES
QH_ml = zeros(size(H));     % matlab solved QH
R_ml = zeros(size(H));      % matlab solved R
Herror = zeros(Hsize);      % error between matlab solved H and original H
Yh = zeros(size(Y));        % Qh*Y, the Y transformed by QR decomposition
X_recover = zeros(size(Y)); % Try to use R^(-1) to recover X

% evaluated parameters
Hvar = 0;                   % 
errorThreshold = 0.0001;
RerrCNT = 0;                % Number of wrong R (error over errorThreshold)
maxRerr = 0;
maxRerrIndex = 0;
Yerr = 0;

for index = (1:dsize)
    [QH_ml(index, :), R_ml(index, :), Yh(index, :)] = Hdecompo(H(index, :), N, Y(index, :));
    rerror = reshape(R(index, :), Hsize) - reshape(R_ml(index, :), Hsize);
    rerrorVal = norm(rerror);
    RerrCNT = RerrCNT + (rerrorVal > errorThreshold);
    if maxRerr < rerrorVal
        maxRerr = rerrorVal;
        maxRerrIndex = index;
    end
    Herrori = (reshape(H(index, :), Hsize) - reshape(QH_ml(index, :), Hsize)' * ...
             reshape(R(index, :), Hsize)).^2;
    Herror = Herror + Herrori;
    Hvar = Hvar + sum(Herrori, 'all');
    Yerr = Y_h_cpp(index, :) - Yh(index, :);
    X_recover(index, :) = ((reshape(R(index, :), Hsize)) \ (Y_h_cpp(index, :))' .* sqrt(10))';
end

Yerror = norm(Y_h_cpp - Yh) / sqrt(dsize * N)
Rerror = norm(R - R_ml) / sqrt(dsize * esize)  % error of R between cpp code and Matlab code
RerrCNT
maxRerr
maxRerrIndex;
Herror = (Herror./ dsize).^(1/2)
Hvar = (Hvar / dsize / 16)^0.5

%% PLOT RECOVERED X

[Xstar, Ystar] = ndgrid(QAM_16_sol_x, QAM_16_sol_x);
X_plot = [X_recover(:, [1 2]); X_recover(:, [3 4])];
color = ([X(:, [1 2]) -3*ones(dsize, 1); X(:, [3 4]) -3*ones(dsize, 1)] + 3) / 6.3;
color = [color(:,1) color(:,3) color(:,2)];

figure
g1 = scatter(X_plot(:,1), X_plot(:, 2), 10, color, 'filled');
hold on;
g2 = scatter(Xstar(:), Ystar(:), 100, 'w', 'LineWidth', 3);
axis equal;
grid on;
title(sprintf('SNR %.2f', SNR))

%% ML detection
[X1, X2, X3, X4, X5, X6, X7, X8] ...
        = ndgrid(QAM_16_sol_x, QAM_16_sol_x, QAM_16_sol_x, QAM_16_sol_x, ...
                 QAM_16_sol_x, QAM_16_sol_x, QAM_16_sol_x, QAM_16_sol_x);
Xall = [X1(:), X2(:), X3(:), X4(:), X5(:), X6(:), X7(:), X8(:)]';

[Xsol, Sol_err] = MaxLike(H, Y, Xall, N);

X_err = abs(X - Xsol)/2;
large_err_point = X_err > 2;
X_err(large_err_point) = 1;
BER = sum(X_err, 'all') / (2 * numel(X_err))

%% test one matrix indexed by 'index'
%-----EDIT THIS-----%
    index = 524;      % <=====
%-------------------%
disp('Test');

[QH_test, R_test] = Hdecompo(H(index, :), N);
Rcpp = reshape(R(index, :), Hsize)
Rtest = reshape(R_test, Hsize)
Hrec = reshape(QH_test, Hsize)' * Rcpp
Horig = reshape(H(index, :), [4 4])
Herr = Horig - Hrec
Rerr = Rtest - Rcpp

%% EXAMPLE
H_test = [1.333 0 -0.6318 0.2621 0 1.333 -0.2621 -0.6318 0.2712 -0.6801 -0.2118 0.9433 0.6801 0.2712 -0.9433 -0.2118];
Htest = reshape(H_test, [4 4])
[QH_test, R_test] = Hdecompo(H_test, 4);

Rtest = reshape(R_test, [4 4])
QHtest = reshape(QH_test, [4 4])
Hrec = (QHtest' * Rtest)
Herr = Htest1 - Hrec

%% ML 

function [Xsol, err] = MaxLike(H, Y, Xi, N)
% H = NxM, Y = Nx1, X = Mx2^M
    Xsol = zeros(size(H, 1), N);
    err = zeros(size(H, 1), 1);
    Xall_norm = Xi ./ sqrt(10);
    for i = (1:size(H, 1))
        Yideal = reshape(H(i,:), [N N]) * Xall_norm;
        errvec = Y(i, :)' - Yideal;
        [err(i), Xindex] = min(sum(errvec .* errvec));
        Xsol(i, :) = Xi(:, Xindex);
    end
end

%% Decomposition function

function [QH, R, Yh] = Hdecompo(H, N, Y)
% Decompose the channel matrix H
    QH = eye(N);
    R = reshape(H, [N N]);
    Yh = Y;
    for j = (1:N-1)
        for i = (j+1:N)
            x = R(j, (j:N));
            y = R(i, (j:N));
            angle = 0;
            if x(1) < 0
                if y(1) > 0
                    angle = -pi/2;
                    temp = y(1);
                    y(1) = -x(1);
                    x(1) = temp;
                elseif y(1) < 0
                    angle = pi/2;
                    temp = y(1);
                    y(1) = x(1);
                    x(1) = -temp;
                end
            else 
                angle = 0;
            end
            angle = angle - atan(y(1)/x(1));
            x(1) = norm([x(1) y(1)]);
            y(1) = 0;
            RotateM = [cos(angle) -sin(angle); sin(angle) cos(angle)];
            yr = RotateM * [Yh(j); Yh(i)];
            Yh(j) = yr(1);
            Yh(i) = yr(2);
            NewXY = RotateM * [x(2:length(x)); y(2:length(y))];
            QH([j i], :) = RotateM * QH([j i], :);
            x(2:length(x)) = NewXY(1, :);
            y(2:length(y)) = NewXY(2, :);
            R(j, (j:N)) = x;
            R(i, (j:N)) = y;
        end
    end
    R = R(:);
    QH = QH(:);
end