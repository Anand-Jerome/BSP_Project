% Create Arduino object
a = arduino('COM4', 'Uno'); % Replace 'COMX' with your actual port

% ecg=load('ecg_hfn.dat');
% fs=1000; % Sampling frequency
% L=length(ecg)
% t=[1:L]/fs;
% figure;
% plot(t,ecg);
% 
% ecg_sel=ecg(1:700);
% x=[0,50,75,105,130,210,215,237,268,292,395,456,498,556,700];
% y=[0,0,0.201075269,0.201075269,0,0,-0.046236559,1.30517921,-0.501433692,0,0,0.59894731,0.59894731,0,0];
% y=max(ecg_sel)*y
% 
% 
% 
% % Generating piecewise linear model
% linModel=[];
% for i=1:14
%     if isequal (y(i+1),y(i));
%         a=y(i)*ones(1,x(i+1)-x(i));
%     else
%         a=y(i):(y(i+1)-y(i))/(x(i+1)-x(i)):y(i+1);
%         a=a(1:end-1);
%     end
%     linModel=[linModel,a];
% end
% t1=(1:length(linModel))/fs;
% figure;
% plot(t1,linModel);
% 
% %PSD of desired signal
% nfft=max(256,2^nextpow2(length(linModel)));
% [Pxx,F]=periodogram(linModel,[],nfft,fs);
% figure;
% plot(F,10*log10(Pxx));
% 
% % PSD of desired signal
% ECG=ecg(2776:2948);
% ECG=ECG-mean(ECG);
% [Pxx1,F]=periodogram(ECG,[],nfft,fs);
% % Average PSD of noise
% ECG2=ecg(4205:4381);
% ECG2=ECG2-mean(ECG2);
% [Pxx2,F]=periodogram(ECG2,[],nfft,fs);
% ECG3=ecg(7051:7230);
% ECG3=ECG3-mean(ECG3);
% [Pxx3,F]=periodogram(ECG3,[],nfft,fs);
% Pxx_avg=(Pxx1+Pxx2+Pxx3)/3;
% figure;
% plot(F,10*log10(Pxx_avg));
% 
% % Transfer function of Wiener filter
% W=zeros(1,length(F));
% for i=1:length(F)
%     W(i)=1/(1+(Pxx_avg(i)/Pxx(i)));
% end
% % W in time domain
% Y=ifftshift(abs(ifft(W,200)));
% Initialize variables

numSamples = 1000; % Adjust the number of samples as needed
% Create a figure for real-time plotting
% h = figure;
k=1;
h = animatedline;
% ax = axes('Parent', h);
W = load("Wiener_Filter_Parameter.mat","Y");
Weight = W.Y;
buffer_size = 20;
ecgData = zeros(buffer_size, 1);

for i = 1:buffer_size:numSamples
    for j = 1:buffer_size    
        ecgValue = readVoltage(a, 'A0');
        ecgData(j) = ecgValue;
        pause(0.001);
    end
    output = conv(ecgData,Weight);
    % Plot real-time data
    % plot(ax, ecgData, 'b', 'LineWidth', 1.5);
    addpoints(h,[k:k+(length(Weight)+buffer_size-2)],output);
    k=k+(length(Weight)+buffer_size-2);
    drawnow
    title('Real-Time ECG Data');
    xlabel('Sample');
    ylabel('Voltage');
    grid on;
    % drawnow;

    pause(0.01); % Adjust the pause time as needed
end

% Close Arduino connection
clear a;
