%% RESET
clear

%% Declare Matrix
MAXDATA = 200000;
BER_var = zeros(12, MAXDATA);

%% open file

SNR = ["3", "5", "7", "10", "13", "15", "17", "20", "22", "25", "28", "30"];
path = "..\data\FPGA_result\" + SNR' + "dB\BER.csv"

for i = (1:12)
	data = csvread(path(i));
	BER_var(i, (1:length(data))) = data;
	if (length(data) < MAXDATA)
		BER_var(i, 10000 + 1: MAXDATA) = data(10000);
	end
end

matdata = matfile('BER.mat', 'Writable', false);
SNRpoint = matdata.SNRpoint;
BER_4best = matdata.BER_4best;
BER_decom_FXP_FWL = matdata.BER_decom_FXP;
linewidth = 1.4;

% Video Parameters
timelen = 5;
fps = 10;
totalframe = timelen * fps;
idx = round(logspace(0, log10(MAXDATA), totalframe));
% idx = round(linspace(1, 200000, totalframe));

%% 
idx = [10 linspace(100, 1000, 10) linspace(2000, 10000, 9), linspace(20000, 200000, 19)];
totalframe = length(idx);

%% GIF Parameters
GIFfilename = "BERvariation.gif";

GIFframedelay = 0.5;

%% Video parameters
vid = VideoWriter('BERvariation.avi');
vid.Quality = 100;
vid.FrameRate = fps;

%% Generate animation

open(vid);

figure
for i = 1:totalframe
	semilogy(matdata.SNRpoint, matdata.BER_4best, ...
				'-^', 'color', 'r', 'LineWidth', linewidth);
	hold on;
	semilogy(SNRpoint, BER_decom_FXP_FWL((1:12),4), ...
				'-s', 'Color', [0 .5 .1], 'LineWidth', linewidth); 
	semilogy(SNRpoint, BER_var(:,idx(i)), ...
				'-s', 'Color', 'b', 'LineWidth', linewidth);
			
	grid on;

	legend('Floating-point K-best', 'Fixed-point K-best (FWL=12)', ...
			num2str(idx(i), 'FPGA implementation n = %d'));
	xlabel('SNR (dB)');
	ylabel('BER');
	title('BER of K-best MIMO Detector');
	
% 	pause(0.02)
	hold off;
	frame = getframe(gcf);
	% output GIF
	im = frame2im(frame);
	[imind, colormap] = rgb2ind(im, 256);	
	if i == 1
		imwrite(imind, colormap, GIFfilename, ...
					'gif', 'DelayTime', GIFframedelay, 'LoopCount', inf);
	else
		imwrite(imind, colormap, GIFfilename, ...
					'gif', 'DelayTime', GIFframedelay,'WriteMode', 'append');
	end
	% output video
	writeVideo(vid, frame);
end

imwrite(imind, colormap, GIFfilename, ...
		'gif', 'DelayTime', 2,'WriteMode', 'append');

vid.close();
