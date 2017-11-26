%read clean speech
clear all;

close all;
[mic1,fs]=audioread('cleanspeech.wav');
%calculate the size of mic1
s1=size(mic1);
t = round(linspace(1,256,256));

%N=input('Enter the frame length:');
N=256;
% N = 164;
% # of frames
M=s1(1);
K=fix(M/N);
preemph=[1 0.63];

table1 = [];
z = 1;
% loop for K times computation
for i=76%overlapping
n =(1:N)+(N*(i-1)/2);
x1=mic1(n);
x=x1.*hamming(N);
x=filter(1,preemph,x);
% compute and plot fft and lpc
X=fft(x);
[a,e]=lpc(x,10)
Xlpc=freqz(1,a);
% Xlpc = Xlpc(1:256);
rr=roots(a); % determines the roots
norm_freq=angle(rr); %frequency in radians
freq_Hz=(norm_freq*fs)/(2*pi) % freq in Hz
freq_Hz1 = abs(freq_Hz);
freq_Hz1 = sort(freq_Hz1);
%pause
subplot(2,2,1);plot(20*log10(abs(X(1:N/2))))% FFT plot
subplot(2,2,2);plot(20*log10(abs(Xlpc)))% LPC Plot 
title(['Frame Number: ' num2str(i)]);
subplot(2,2,3);plot(x);% time domain plot
subplot(2,2,4);zplane(1,a);% pole zero plot
%pause
table1(1,i) = i;
table1(2,z) = freq_Hz1(1);
table1(3,z) = freq_Hz1(3);
table1(4,z) = freq_Hz1(5);
z = z+1;
end

i =1;

for i=1:z-1
    temp = 0;
    if table1(2,i) == 0
        table1(2,i) = table1(4,i);
    end
    if table1(2,i)>=table1(3,i)
        temp = table1(2,i);
        table1(2,i) = table1(3,i);
        table1(3,i) = temp;
    end
end


figure();
plot(table1(1,:),table1(2,:)); 
hold on 
plot(table1(1,:),table1(3,:)); 
title('Formant frequency versus Frame Number');
xlabel('Frame Number'); ylabel('Frequency');
legend('Formant 1', 'Formant 2');


% i = 1;
% m = 1;
% for i=1:z+2
%     
%     table1(4,m) = (table1(2,m+i-1) + table1(2,m+i))*0.5;
%     table1(5,m) = (table1(3,m+i-1) + table1(3,m+i))*0.5;
%     m = m+1;
% end
% 
