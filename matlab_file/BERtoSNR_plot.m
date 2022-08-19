% %%
% save('BER.mat', 'SNRpoint', 'BER_ML', 'BER_2best', 'BER_2best_final', ...
%         'BER_2best_withSortedMatrix', 'BER_2best_withZigZag', 'BER_3best');
%%
clear

%% Read Data
matdata = matfile('BER.mat', 'Writable',true);

%% SAVEV
% matdata.BER_FPGA_16bit = zeros(12, 1);

%% Load
SNRpoint = matdata.SNRpoint;
BER_4best_q = matdata.BER_4best_q;
BER_4best = matdata.BER_4best;
BER_decom_FXP_FWL = matdata.BER_decom_FXP;
BER_4best_PED_FWL = matdata.BER_4best_PED_FWL;
BER_FPGA_15b = matdata.BER_FPGA_15bit;
BER_FPGA_16b = matdata.BER_FPGA_16bit;

% plot property
color = ['#0072BD'; '#D95319'; '#EDB120'; '#77AC30'; '#4DBEEE'; '#7E2F8E'; '#A2142F';...
         '#FF00FF'; '#00FF00'; '#FFFF00'; '#0000FF'; '#FF0000'];
linewidth = 1.4;    

%%
figure
semilogy(matdata.SNRpoint, matdata.BER_ML, '-ok', 'LineWidth', linewidth); hold on;
% semilogy(matdata.SNRpoint, matdata.BER_2best, '-ob');
% semilogy(matdata.SNRpoint, matdata.BER_2best_withZigZag, '-sb');
% semilogy(matdata.SNRpoint, matdata.BER_2best_withSortedMatrix, '-xb');
semilogy(matdata.SNRpoint, matdata.BER_2best_final, '-^b', 'LineWidth', linewidth);
semilogy(matdata.SNRpoint, matdata.BER_3best, '-^r', 'LineWidth', linewidth);
semilogy(matdata.SNRpoint, matdata.BER_4best, '-^', 'color', [0 .7 .2], 'LineWidth', linewidth);
% semilogy(matdata.SNRpoint, matdata.BER_optimized, '-sm', 'LineWidth', linewidth);

axis([0 32 1e-7 1]);
grid on;
legend('ML', 'K-best (K=2)', 'K-best (K=3)', 'K-best (K=4)', 'K-best (K=4, 16, 3, 3, ..., 3)');
%         'K-best (K=2) with zigzag sorting', ...   
%         'K-best (K=2) with sorted matrix', ...

xlabel('SNR (dB)');
ylabel('BER');
title('4x4 MIMO 16-QAM');

%% 
BERdecom = matdata.BER_vs_bit_decom;
figure
subplot(2, 1, 1);
semilogy(BERdecom(:,1), BERdecom(:,2), ...
            '-ob', 'LineWidth', linewidth);
axis([5.5 12.5 1e-4 1]);
grid on;
xlabel('Fraction word length');
ylabel('BER');
title('BER of quantilized decomposition with different FWLs');

BERwhole = matdata.BER_vs_bit_whole;
subplot(2, 1, 2);
semilogy(BERwhole(:,1), BERwhole(:,2), ...
            '-ob', 'LineWidth', linewidth);
axis([7.5 12.5 1e-4 1]);
grid on;
xlabel('Fraction word length');
ylabel('BER');
title('BER of quantilized detector with different FWLs');

%% quantized FWL for main numbers

figure
for i = 1:9
semilogy(SNRpoint, BER_4best_q((1:12),i), '-^',...
         'Color', color(i,:),'LineWidth', linewidth);   hold on;
end
semilogy(SNRpoint, BER_4best, '--ok', 'LineWidth', linewidth); 

axis([2 31 1e-4 1]);
grid on;
legendstring = sprintf('Fixed-point K-best (K=4) with FWL=') + string(BER_4best_q(13,(1:9))');

legend([legendstring; 'Floating-point K-best (K=4)']);
xlabel('SNR (dB)');
ylabel('BER');
title('BER of floating-point and fixed-point design');

%% FWL of column power

figure
semilogy(SNRpoint, BER_decom_FXP_FWL((1:12),1), '-s',...
         'Color', 'r', 'LineWidth', linewidth);   hold on;
for i = 2:7
semilogy(SNRpoint, BER_decom_FXP_FWL((1:12),i), '-^',...
         'Color', color(i,:), 'LineWidth', linewidth);
end
semilogy(SNRpoint, BER_4best, '--ok', 'LineWidth', linewidth);

axis([2 31 1e-4 1]);
grid on;
legendstring = sprintf('Column norm FWL=') + string(BER_decom_FXP_FWL(13,(2:7))');

legend(['Fixed-point K-best without column sorting'; legendstring; 'Floating-point K-best (K=4)']);
xlabel('SNR (dB)');
ylabel('BER');
title('BER of fixed-point K-best with different FWL of column norm');

%% FWL of PED

figure
for i = 7:-1:1
semilogy(SNRpoint, BER_4best_PED_FWL((1:12),i), '-^',...
         'Color', color(i,:), 'LineWidth', linewidth);   hold on;
end
semilogy(SNRpoint, BER_4best, '--ok', 'LineWidth', linewidth);

axis([2 31 1e-4 1]);
grid on;
legendstring = sprintf('PED FWL=') + string(BER_4best_PED_FWL(13,(7:-1:1))');

legend([legendstring; 'Floating-point K-best (K=4)']);
xlabel('SNR (dB)');
ylabel('BER');
title('BER of fixed-point K-best with different FWL of PED');

%% BER of FPGA
figure
semilogy(SNRpoint, BER_FPGA_15b((1:12),2), '-s', 'Color', [0 .5 1], 'LineWidth', linewidth); hold on;
semilogy(SNRpoint, BER_FPGA_16b((1:12),1), '-s', 'Color', 'b', 'LineWidth', linewidth); %hold on;
semilogy(SNRpoint, BER_decom_FXP_FWL((1:12),4), '-s', 'Color', [0 .5 0 .5], 'LineWidth', linewidth); 
semilogy(matdata.SNRpoint, matdata.BER_4best, '-^', 'color', [1 0 0 .5], 'LineWidth', linewidth);
% semilogy(SNRpoint, BER_4best_q((1:12),6), '-s', 'Color', 'g', 'LineWidth', linewidth); 
% semilogy(SNRpoint, BER_decom_FXP_FWL((1:12),5), '-^', 'Color', color(1,:), 'LineWidth', linewidth);
% semilogy(SNRpoint, BER_FPGA((1:12),1), '-s', 'Color', 'b', 'LineWidth', linewidth); 
semilogy(matdata.SNRpoint, matdata.BER_ML, '-ok', 'LineWidth', linewidth);
grid on;

legend('FPGA implementation last week', 'FPGA implementation this week'...
		, 'Fixed-point K-best (FWL=12)', 'Floating-point K-best', 'ML');
xlabel('SNR (dB)');
ylabel('BER');
title('BER of K-best MIMO Detector');
