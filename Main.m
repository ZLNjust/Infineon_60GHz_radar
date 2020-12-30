clc
clear all
close all




[FileName,PathName] = uigetfile('*.txt','Select the .txt file from radar');

[Radar_Parameter,Frame_Number,NumRXAntenna,Frame]=data_import([PathName FileName]);


%% Singal Processing
c = 3e8; % Speed of light (m/s)
fc = (Radar_Parameter.Lower_RF_Frequency_kHz+Radar_Parameter.Upper_RF_Frequency_kHz)/2*1000; % Center frequency (Hz)
CPI = 1/250; % Coherent processing interval (s). Seconds to be processed (CPI must be lower to N*fs)
CRR = 1/Radar_Parameter.Chirp_Time_sec; % Chirp repetition rate (Hz)
FRR=1/Radar_Parameter.Frame_Period_sec;% Frame repetition rate (Hz)
BW = (Radar_Parameter.Upper_RF_Frequency_kHz-Radar_Parameter.Lower_RF_Frequency_kHz)*1000; % Bandwidth (Hz)
gamma = BW*CRR; % Chirp rate (Hz/s)
range_res = c/(2*BW);
max_range = range_res*fix(Radar_Parameter.Sampling_Frequency_kHz*1e3/CRR)/2;
NL = 1024;


%% FFT
for RXAntenna = 1: NumRXAntenna  % RX 
    
    tau = (0:(Frame_Number-1))/FRR; % Slow time (s)
    raw_data_matrix =zeros(Frame_Number,Radar_Parameter.Samples_per_Chirp);
    for FrameN = 1:Frame_Number
        raw_data_matrix((FrameN),:)= Frame(:,1,RXAntenna,FrameN)';
    end
    
    % subtract DC
    avgDC=nanmean(raw_data_matrix,2);               
    for jj = 1:size(raw_data_matrix,1)
        raw_data_matrix(jj,:) = raw_data_matrix(jj,:) - avgDC(jj);
    end
    
    [a,b]= size(raw_data_matrix);
    eje_dis = (0:NL-1)/NL*c*Radar_Parameter.Sampling_Frequency_kHz*1000/2/gamma;
    win=rectwin(Radar_Parameter.Samples_per_Chirp);
    win_2=win(:,ones(Frame_Number,1)); %add window
    raw_data_matrix_2 = raw_data_matrix.*win_2';
    per = fft(raw_data_matrix_2,NL,2); % Range profiles
    raw_data_matrix_antenna(:,:,RXAntenna)=raw_data_matrix_2;
    % MTI处理
    per_MTI=diff(per,1,1);
    per_data_matrix_antenna(:,:,RXAntenna)=per_MTI;    
end


startRange = 0.1;
startPoint = fix(startRange/(max_range*2/NL));
stopRange =4;
stopPoint=fix(stopRange/(max_range*2/NL));


figure;imagesc(eje_dis(10:NL/2),tau,20*log10(abs(per_data_matrix_antenna(:,10:NL/2,1))/max(max(abs(per_data_matrix_antenna(:,10:NL/2,1))))),[-50 0]); colorbar; colormap(gca,jet);

colormap(gca,jet);
xlabel('Range (m)');
ylabel('Slow time (s)');
title("2D Time-Range-Map Channal 1",'FontSize', 18,'FontWeight', 'bold');

figure;
imagesc(eje_dis(10:NL/2),tau,20*log10(abs(per_data_matrix_antenna(:,10:NL/2,2))/max(max(abs(per_data_matrix_antenna(:,10:NL/2,2))))),[-50 0]); colorbar; colormap(gca,jet);

colormap(gca,jet);
xlabel('Range (m)');
ylabel('Slow time (s)');
title("2D Time-Range-Map Channal 2",'FontSize', 18,'FontWeight', 'bold');

figure;
imagesc(eje_dis(10:NL/2),tau,20*log10(abs(per_data_matrix_antenna(:,10:NL/2,3))/max(max(abs(per_data_matrix_antenna(:,10:NL/2,3))))),[-50 0]); colorbar; colormap(gca,jet);

colormap(gca,jet);
xlabel('Range (m)');
ylabel('Slow time (s)');
title("2D Time-Range-Map Channal 3",'FontSize', 18,'FontWeight', 'bold');



%%DBF

%%矫正参数
A_Calibration=[1;0.9468;1.0988];
P_Calibration=[1;0.0968 - 0.9953i;0.1478 - 0.9890i];
Wr=[1;1;1];

for NR= 1: NumRXAntenna
   per_data_matrix_antenna_2(NR,:)= reshape(per_data_matrix_antenna(:,:,NR)',1,[]);
end


%% 360频扫
theta=linspace(0,180,360);
for  j=1:length(theta)
    per_DBF=DBF(A_Calibration,P_Calibration,Wr,per_data_matrix_antenna)
    per_DBF_ave=mean(abs(per_DBF));
    per_DBF_ave_matrix(j,:)=per_DBF_ave;
end

%%VitalSign Detection
location=1.6;
[Resipration,Heartbeat]=vitalsign_extract(per,location,Radar_Parameter,NL);

