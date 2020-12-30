function [Resipration,Heartbeat]=vitalsign_extract(Phase,location,Radar_Parameter,NL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Vital Sign extaction                           %
%                                                 %
%  Li Zhang  & Prateek                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e8; % Speed of light (m/s)
CRR = 1/Radar_Parameter.Chirp_Time_sec; % Chirp repetition rate (Hz)
BW = (Radar_Parameter.Upper_RF_Frequency_kHz-Radar_Parameter.Lower_RF_Frequency_kHz)*1000; % Bandwidth (Hz)
range_res = c/(2*BW);
max_range = range_res*fix(Radar_Parameter.Sampling_Frequency_kHz*1e3/CRR)/2;

VitalRange = location;
VitalPoint = fix(VitalRange/(max_range*2/NL));


VitalSigns = Phase(:,VitalPoint);

Fs=1/Radar_Parameter.Frame_Period_sec;

%%DACM
I_signal=real(VitalSigns);
Q_signal=imag(VitalSigns);
diff_I=diff(I_signal);
diff_Q=diff(Q_signal);
diff_Signal=(I_signal(1:end-1).*diff_Q-diff_I.*Q_signal(1:end-1));
for i = 2:length(diff_I)
    Rawsignal(i-1)=sum(diff_Signal(1:i));
end

%%呼吸信号
Resipration = highpass(Rawsignal,0.15,Fs,'Steepness',0.85,'StopbandAttenuation',85);
%Resipration = smoothdata(Rest_highpass,'lowess','SmoothingFactor',0.45,'SamplePoints',tau);

%心跳信号
Heartbeat = bandpass(Rawsignal,[0.5 2],Fs,'Steepness',0.85,'StopbandAttenuation',60);

end