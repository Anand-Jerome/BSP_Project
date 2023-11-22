% MATLAB PROGRAM ecg_hfn.m
clear all               % clears all active variables
close all

% Collect data

data_hfn = load('ecg_hfn.dat');
data_lfn = load('ecg_lfn.dat');
nf = load('pec1.dat');
data_nf= nf(:,2);

%sampling rate
fs = 1000; %sampling rate = 1000 Hz

% calculate length of data
slen_hfn = length(data_hfn);
slen_lfn = length(data_lfn);
slen_nf = length(data_nf);

% % Crop 1000 data points from both ends
% crop_length = 1000;
% data_nf = nf(crop_length:slen_nf-crop_length,2);
% slen_nf = slen_nf - 2*crop_length + 1;

% Plot all 3 signals in time-domain and their FFTs

% Define time vectors for each signal
timeVector_hfn = (0:(slen_hfn - 1)) * (1 / fs);
timeVector_lfn = (0:(slen_lfn - 1)) * (1 / fs);
timeVector_nf = (0:(slen_nf - 1)) * (1 / fs);

% Plot the signals in the time domain
figure;
subplot(3, 2, 1);
plot(timeVector_hfn, data_hfn);
title('High-Frequency Noise');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(3, 2, 3);
plot(timeVector_lfn, data_lfn);
title('Low-Frequency Noise');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(3, 2, 5);
plot(timeVector_nf, data_nf);
title('Noise-Free Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Compute FFTs
fft_hfn = fft(data_hfn);
fft_lfn = fft(data_lfn);
fft_nf = fft(data_nf);

% Compute corresponding frequency vectors
frequencies_hfn = linspace(0, fs / 2, slen_hfn );
frequencies_lfn = linspace(0, fs / 2, slen_lfn );
frequencies_nf = linspace(0, fs / 2, slen_nf );

% Plot the FFTs
subplot(3, 2, 2);
plot(frequencies_hfn, abs(fft_hfn));
title('FFT of High-Frequency Noise');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0, fs/2]); % Set x-axis limits

grid on;

subplot(3, 2, 4);
plot(frequencies_lfn, abs(fft_lfn));
title('FFT of Low-Frequency Noise');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0, fs/2]); % Set x-axis limits
grid on;

subplot(3, 2, 6);
plot(frequencies_nf, abs(fft_nf));
title('FFT of Noise-Free Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0, fs/2]); % Set x-axis limits

grid on;

% Adjust layout
sgtitle('Time Domain and FFT Plots');

% Time-domain and frequency-domain filtering of hfn data.

% Implement 8 point MAF filter (Time-domain)
maf_data_hfn = [data_hfn];

for i = 8:(slen_hfn)
    maf_data_hfn(i)=0.125*(data_hfn(i)+data_hfn(i-1)+data_hfn(i-2)+data_hfn(i-3)+data_hfn(i-4)+data_hfn(i-5)+data_hfn(i-6)+data_hfn(i-7));
end

% Implement butterworth filter (Frequency-domain)

% Define the filter order
filterOrder = 4;

% Define the cutoff frequencies
lowFreq = 0.5;  % Lower cutoff frequency
highFreq = 50;  % Upper cutoff frequency

% Design the Butterworth bandpass filter
[b, a] = butter(filterOrder, [lowFreq, highFreq] / (fs / 2), 'bandpass');



% Apply the filter to the noisy ECG signal
butter_data_hfn = filter(b, a, data_hfn);


% Plot the results for the Moving-Average Filter
figure;

subplot(2, 2, 1);
plot(timeVector_hfn, data_hfn);
title('Original High-Frequency Noise');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(2, 2, 3);
plot(timeVector_hfn, maf_data_hfn);
title('Moving-Average Filter Output');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Compute FFT for Moving-Average Filter Output
fft_maf_data_hfn = fft(maf_data_hfn);
frequencies_maf_hfn = linspace(0, fs / 2, slen_hfn / 2);

subplot(2, 2, 2);
plot(frequencies_maf_hfn, abs(fft_maf_data_hfn(1:slen_hfn/2)));
title('FFT of Moving-Average Filter Output');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

% Plot the results for the Butterworth Filter
subplot(2, 2, 4);
plot(timeVector_hfn, butter_data_hfn);
title('Butterworth Filter Output');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Compute FFT for Butterworth Filter Output
fft_butter_data_hfn = fft(butter_data_hfn);
frequencies_butter_hfn = linspace(0, fs / 2, slen_hfn / 2);

figure;
subplot(2, 1, 1);
plot(frequencies_butter_hfn, abs(fft_butter_data_hfn(1:slen_hfn/2)));
title('FFT of Butterworth Filter Output');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
% Plot a zoomed-in version to better visualize the filter characteristics
plot(frequencies_butter_hfn, abs(fft_butter_data_hfn(1:slen_hfn/2)));
title('Zoomed-In FFT of Butterworth Filter Output');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0, 10]); % Adjust the range as needed
grid on;

% Adjust layout
sgtitle('Time Domain and FFT Plots after Filtering HF data');

% Time-domain and frequency-domain filtering of lfn data.

% Implement 8 point MAF filter

maf_data_lfn = [data_lfn];

for i = 8:(slen_lfn)
    maf_data_lfn(i)=0.125*(data_lfn(i)+data_lfn(i-1)+data_lfn(i-2)+data_lfn(i-3)+data_lfn(i-4)+data_lfn(i-5)+data_lfn(i-6)+data_lfn(i-7));
end

% Implement a derivative filter

der_data = [maf_data_lfn];

for i = 3:(slen_lfn)
    der_data(i)=500*(maf_data_lfn(i)-maf_data_lfn(i-2));
end

% Implement butterworth filter (Frequency-domain)

% Define the filter order
filterOrder = 4;

% Define the cutoff frequencies
lowFreq = 0.5;  % Lower cutoff frequency
highFreq = 50;  % Upper cutoff frequency

% Design the Butterworth bandpass filter
[b, a] = butter(filterOrder, [lowFreq, highFreq] / (fs / 2), 'bandpass');



% Apply the filter to the noisy ECG signal
butter_data_lfn = filter(b, a, data_lfn);

% Plot the results for the Derivative Filter
figure;

subplot(2, 2, 1);
plot(timeVector_lfn, data_lfn);
title('Original Low-Frequency Noise');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(2, 2, 3);
plot(timeVector_lfn, der_data);
title('Derivative Filter Output');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Compute FFT for Derivative Filter Output
fft_der_data_lfn = fft(der_data);
frequencies_der_lfn = linspace(0, fs / 2, slen_lfn );

subplot(2, 2, 2);
plot(frequencies_der_lfn, abs(fft_der_data_lfn));
title('FFT of Derivative Filter Output');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0, fs/2]); 
grid on;

% Plot the results for the Butterworth Filter
subplot(2, 2, 4);
plot(timeVector_lfn, butter_data_lfn);
title('Butterworth Filter Output');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Compute FFT for Butterworth Filter Output
fft_butter_data_lfn = fft(butter_data_lfn);
frequencies_butter_lfn = linspace(0, fs / 2, slen_lfn);

figure;
subplot(2, 1, 1);
plot(frequencies_butter_lfn, abs(fft_butter_data_lfn));
title('FFT of Butterworth Filter Output');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0, fs/2]); 
grid on;

subplot(2, 1, 2);
% Plot a zoomed-in version to better visualize the filter characteristics
plot(frequencies_butter_lfn, abs(fft_butter_data_lfn));
title('Zoomed-In FFT of Butterworth Filter Output');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0, 10]); % Adjust the range as needed
grid on;

% Adjust layout
sgtitle('Time Domain and FFT Plots after Filtering LF data');

% Wiener filtering of HFN data
ecg= [data_hfn];
ecg_sel=ecg(1:700);
x=[0,50,75,105,130,210,215,237,268,292,395,456,498,556,700];
y=[0,0,0.201075269,0.201075269,0,0,-0.046236559,1.30517921,-0.501433692,0,0,0.59894731,0.59894731,0,0];
y=max(ecg_sel)*y

% Generating piecewise linear model
linModel=[];
for i=1:14
    if isequal (y(i+1),y(i));
        a=y(i)*ones(1,x(i+1)-x(i));
    else
        a=y(i):(y(i+1)-y(i))/(x(i+1)-x(i)):y(i+1);
        a=a(1:end-1);
    end
    linModel=[linModel,a];
end
t1=(1:length(linModel))/fs;
figure;
plot(t1,linModel);

%PSD of desired signal
nfft=max(256,2^nextpow2(length(linModel)));
[Pxx,F]=periodogram(linModel,[],nfft,fs);
figure;
plot(F,10*log10(Pxx));

% PSD of desired signal
ECG=ecg(2776:2948);
ECG=ECG-mean(ECG);
[Pxx1,F]=periodogram(ECG,[],nfft,fs);
% Average PSD of noise
ECG2=ecg(4205:4381);
ECG2=ECG2-mean(ECG2);
[Pxx2,F]=periodogram(ECG2,[],nfft,fs);
ECG3=ecg(7051:7230);
ECG3=ECG3-mean(ECG3);
[Pxx3,F]=periodogram(ECG3,[],nfft,fs);
Pxx_avg=(Pxx1+Pxx2+Pxx3)/3;
figure;
plot(F,10*log10(Pxx_avg));


% Transfer function of Wiener filter
W=zeros(1,length(F));
for i=1:length(F)
    W(i)=1/(1+(Pxx_avg(i)/Pxx(i)));
end
% W in time domain
Y=ifftshift(abs(ifft(W,200)));
output=conv(ecg,Y);
t=[1:length(output)]/fs
plot(t,output)
title('Wiener filter output of HFN data');

% Assuming 'output' is the vector containing the output of the Wiener filter

% Compute FFT of the Wiener filter output
fft_output = fft(output);
frequencies_output = linspace(0, fs / 2, length(output) / 2);

% Plot the magnitude spectrum
figure;
plot(frequencies_output, abs(fft_output(1:length(output)/2)));
title('FFT of Wiener Filter Output for HFN data');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% Wiener filtering of LFN data
ecg= [data_lfn];
ecg_sel=ecg(1:700);
x=[0,50,75,105,130,210,215,237,268,292,395,456,498,556,700];
y=[0,0,0.201075269,0.201075269,0,0,-0.046236559,1.30517921,-0.501433692,0,0,0.59894731,0.59894731,0,0];
y=max(ecg_sel)*y


%PSD of desired signal
nfft=max(256,2^nextpow2(length(linModel)));
[Pxx,F]=periodogram(linModel,[],nfft,fs);
figure;
plot(F,10*log10(Pxx));

% PSD of desired signal
ECG=ecg(2776:2948);
ECG=ECG-mean(ECG);
[Pxx1,F]=periodogram(ECG,[],nfft,fs);
% Average PSD of noise
ECG2=ecg(4205:4381);
ECG2=ECG2-mean(ECG2);
[Pxx2,F]=periodogram(ECG2,[],nfft,fs);
ECG3=ecg(7051:7230);
ECG3=ECG3-mean(ECG3);
[Pxx3,F]=periodogram(ECG3,[],nfft,fs);
Pxx_avg=(Pxx1+Pxx2+Pxx3)/3;
figure;
plot(F,10*log10(Pxx_avg));

% Transfer function of Wiener filter
W=zeros(1,length(F));
for i=1:length(F)
    W(i)=1/(1+(Pxx_avg(i)/Pxx(i)));
end
% W in time domain
Y=ifftshift(abs(ifft(W,200)));
output=conv(ecg,Y);
t=[1:length(output)]/fs
plot(t,output)
title('Wiener filter output of LFN data');

% Assuming 'output' is the vector containing the output of the Wiener filter

% Compute FFT of the Wiener filter output
fft_output = fft(output);
frequencies_output = linspace(0, fs / 2, length(output) / 2);

% Plot the magnitude spectrum
figure;
plot(frequencies_output, abs(fft_output(1:length(output)/2)));
title('FFT of Wiener Filter Output for LFN data');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;



% % Calculate the time step
% timeStep = 1 /fs;
% 
% % Determine the number of data points
% numDataPoints = length(filteredECG_butter); % Replace with the length of your ECG data
% 
% % Create the time vector
% timeVector = (0:(numDataPoints - 1)) * timeStep;
% 
% % Plot the filtered ECG signal
% figure;
% plot(timeVector, filteredECG_butter);
% title('Filtered ECG Signal');
% xlabel('Time (seconds)');
% ylabel('Amplitude');
% grid on;
% 
% % Calculate the time vector
% timeStep = 1 / fs;
% numDataPoints = length(data_nf);
% timeVector = (1:21485) * timeStep;
% 
% % Plot the relatively noise-free ECG signal
% figure;
% plot(timeVector, data_nf);
% title('Relatively Noise-Free ECG Signal');
% xlabel('Time (seconds)');
% ylabel('Amplitude');
% grid on;

% 
% % Define the sampling frequency (in Hz)
% fs = 1000; % Replace with the actual sampling frequency of your ECG data
% 
% % Calculate the time vector
% timeStep = 1 / fs;
% numDataPoints = length(data);
% timeVector = (0:(numDataPoints - 1)) * timeStep;
% 
% % Extract the ECG data from the second line
% startIdx = 1; % Adjust the start and end indices to match your data structure
% endIdx = numDataPoints; % These indices are for illustration
% ecgData = data(startIdx:endIdx);
% 
% % Plot the ECG data
% figure;
% plot(timeVector, ecgData);
% title('ECG Data from the .dat File');
% xlabel('Time (seconds)');
% ylabel('Amplitude');
% grid on;


% % Assuming originalECG is your original ECG signal, and filteredECG is your filtered signal
% 
% % Compute the FFT of the original ECG signal
% fftOriginal = fft(data_hfn);
% % Compute the corresponding frequency values for the FFT
% frequencies = linspace(0, fs, length(fftOriginal));
% 
% % Compute the FFT of the filtered ECG signal
% fftFiltered = fft(filteredECG);
% 
% % Plot the FFT of the original ECG signal
% figure;
% subplot(2,1,1);
% plot(frequencies, abs(fftOriginal));
% title('FFT of Original ECG Signal');
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% grid on;
% 
% % Plot the FFT of the filtered ECG signal
% subplot(2,1,2);
% plot(frequencies, abs(fftFiltered));
% title('FFT of Filtered ECG Signal');
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% grid on;
% 


% Implement hanning filter
% hanning_data = [data_hfn];
% 
% for i = 3:(slen_hfn)
%     hanning_data(i)=0.25*(data_hfn(i)+2*data_hfn(i-1)+data_hfn(i-2));
% end
% 
% Implement 8 point MAF filter
% maf_data = [data_hfn];
% 
% for i = 8:(slen_hfn)
%     maf_data(i)=0.125*(data_hfn(i)+data_hfn(i-1)+data_hfn(i-2)+data_hfn(i-3)+data_hfn(i-4)+data_hfn(i-5)+data_hfn(i-6)+data_hfn(i-7));
% end
% 
% Plot everything
% 
% plot(freq_nf, abs(FFT_NF));
% title('Magnitude Spectrum Noise free');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% 
% hold on
% plot(freq_lfn, abs(FFT_LF));
% title('Magnitude Spectrum lfn');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% 
% legend('Noise free','Low freq noise');
% 
% tiledlayout(1,3);
% 
% nexttile
% plot(freq_hfn, abs(FFT_1));
% title('Magnitude Spectrum hfn');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% 
% nexttile
% plot(freq_lfn, abs(FFT_LF));
% title('Magnitude Spectrum lfn');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% 
% nexttile
% hold off
% FFT_3 = fft(hanning_data);
% plot(freq_hfn, abs(FFT_3));
% title('Magnitude Spectrum Hanning filter');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% 
% nexttile
% FFT_4 = fft(maf_data);
% plot(freq_hfn, abs(FFT_4));
% title('Magnitude Spectrum 8 MAF filter');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% 
% legend('hfn','hanning filter','MAF')
% 
% 
% t_hfn=[1:slen_hfn]/fs;
% t_lfn=[1:slen_lfn]/fs;
% t_nf = [1:slen_nf]/fs;
% 
% size(t)
% tiledlayout(1,3);
% 
% nexttile
% plot(t_lfn, data_lfn)
% axis tight;
% title("Raw data");
% xlabel('Time in seconds');
% ylabel('ECG');
% 
% plot(t_nf, data_nf);
% axis tight;
% title("Noise free data");
% xlabel('Time in seconds');
% ylabel('ECG');
% 
% 
% hold on
% 
% nexttile
% plot(t, hanning_data);
% axis tight;
% title("Hanning filter");
% xlabel('Time in seconds');
% ylabel('ECG');
% 
% nexttile
% plot(t, maf_data);
% axis tight;
% title("MAF");
% xlabel('Time in seconds');
% ylabel('ECG');
% 
% legend('raw','han','maf');