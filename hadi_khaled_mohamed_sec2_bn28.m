%% initializing
clc;
close all
clearvars
save_data=true;
%save_data=false;
%% reading matrix data
data=readmatrix('time_pure_noisy.txt');
%% extracting signal data
t_vec=data(:,1);
s0=data(:,2);
s_all=data(:,3);
signals=data(:,2:end);
dt=t_vec(2)-t_vec(1);
fs=1/dt;
T=t_vec(end);
df=1/T;
N=length(t_vec);
%% plotting in the time domain
figure();
plot(t_vec,[s0,s_all]); sgtitle('pure,and all noise figure in time domain')
legend("pure","pure + nb + wgn");
xlabel('t[sec]')
%% ploting in the frequency domain
f_1side=(0:N/2)*df';
y=fft(signals);
p2= abs(y/N);
p1=p2(1:N/2+1,:);
p1(2:end-1,:)=2*p1(2:end-1,:);
figure(); sgtitle('pure,and all noise figure in frequency domain')
loglog(f_1side,p1);
legend("pure","pure + nb + wgn");
ylim([10e-4,10]);
xlabel('Hz')
%% designing first butterworth filter
n_butter_l_1=8;
fc_butter_l_1=20*2*pi;
[b_butter_l_1,a_butter_l_1]=butter(n_butter_l_1,fc_butter_l_1,"low",'s');
butter_l_1=tf(b_butter_l_1,a_butter_l_1);
n_butter_h_1=8;
fc_butter_h_1=10*2*pi;
[b_butter_h_1,a_butter_h_1]=butter(n_butter_h_1,fc_butter_h_1,"high",'s');
butter_h_1=tf(b_butter_h_1,a_butter_h_1);
bpf_butter_1=butter_l_1*butter_h_1;
%% designing second butterworth filter
n_butter_l_2=7;
fc_butter_l_2=20*2*pi;
[b_butter_l_2,a_butter_l_2]=butter(n_butter_l_2,fc_butter_l_2,"low",'s');
butter_l_2=tf(b_butter_l_2,a_butter_l_2);
n_butter_h_2=7;
fc_butter_h_2=10*2*pi;
[b_butter_h_2,a_butter_h_2]=butter(n_butter_h_2,fc_butter_h_2,"high",'s');
butter_h_1=tf(b_butter_h_2,a_butter_h_2);
bpf_butter_2=butter_l_2*butter_h_1;
%% filtering signal using 1st butter worth
s_filt_butter_1=lsim(bpf_butter_1,s_all,t_vec);
%% filtering signal using 2nd butterworth
s_filt_butter_2=lsim(bpf_butter_2,s_all,t_vec);
%% designing 1st chebyshev 1 filter
n_cheb_l_1=5;
fc_cheb_l_1=14*2*pi;
rp=3;
[b_cheb_l_1,a_cheb_l_1]=cheby1(n_cheb_l_1,rp,fc_cheb_l_1,"low",'s');
cheb_tf_l_1=tf(b_cheb_l_1,a_cheb_l_1);
n_cheb_h_1=5;
fc_cheb_h_1=13*2*pi;
[b_cheb_h_1,a_cheb_h_1]=cheby1(n_cheb_h_1,rp,fc_cheb_h_1,"high",'s');
cheb_tf_h_1=tf(b_cheb_h_1,a_cheb_h_1);
cheb_tf_1=cheb_tf_h_1*cheb_tf_l_1;
%% designing 2nd chebyshev 1 filter
n_cheb_l_2=6;
fc_cheb_l_2=14*2*pi;
rp=3;
[b_cheb_l_2,a_cheb_l_2]=cheby1(n_cheb_l_2,rp,fc_cheb_l_2,"low",'s');
cheb_tf_l_1=tf(b_cheb_l_2,a_cheb_l_2);
n_cheb_h_2=6;
fc_cheb_h_2=13*2*pi;
[b_cheb_h_2,a_cheb_h_2]=cheby1(n_cheb_h_2,rp,fc_cheb_h_2,"high",'s');
cheb_tf_h_2=tf(b_cheb_h_2,a_cheb_h_2);
cheb_tf_2=cheb_tf_h_2*cheb_tf_l_1;
%wp_cheb_2=[fc_cheb_h_2 fc__cheb_l_2];
%[b,a]=cheby1(6,3,wp_cheb_2,"bandpass","s");
%cheb_tf_2=tf(b,a);
%% filtering signal using 1st chebyshev1 filter
s_filt_cheb_1=lsim(cheb_tf_1,s_all,t_vec);
%% filtering signal using 2nd chebyshev1 filter 
s_filt_cheb_2=lsim(cheb_tf_2,s_all,t_vec);
%% fourier plot of 1st butterworth filter
f_2side=(0:N/2)*df';
y_f=fft(s_filt_butter_1);
p2_f= abs(y_f/N);
p1_f=p2_f(1:N/2+1,:);
p1_f(2:end-1,:)=2*p1_f(2:end-1,:);
figure(); sgtitle("1st butterworth vs noisy in frequency domain")
loglog(f_2side,[p1(:,2),p1_f])
legend("all noise","filtered")
ylim([10e-4,10]);
xlabel('Hz')
%figure()
%loglog(f_2side,[p1(:,1),p1_f])
%legend("pure","filtered")
%xlabel('Hz')
%ylim([10e-4,10]);
%% fourier plot of 2nd butterworth filter
f_2side=(0:N/2)*df';
y_f=fft(s_filt_butter_2);
p2_f= abs(y_f/N);
p1_f=p2_f(1:N/2+1,:);
p1_f(2:end-1,:)=2*p1_f(2:end-1,:);
figure(); sgtitle("2nd butterworth vs noisy in frequency domain")
loglog(f_2side,[p1(:,2),p1_f])
legend("all noise","filtered")
ylim([10e-4,10]);
xlabel('Hz')
%% fourier plot of 1st chebyshev1 filtered signal
f_2side=(0:N/2)*df';
y_f=fft(s_filt_cheb_1);
p2_f= abs(y_f/N);
p1_f=p2_f(1:N/2+1,:);
p1_f(2:end-1,:)=2*p1_f(2:end-1,:);
figure(); sgtitle("1st chebysheb1 vs noisy in frequency domain")
loglog(f_2side,[p1(:,2),p1_f])
legend("all noise","filtered")
ylim([10e-4,10]);
xlabel('Hz')
%% fourier plot of 2nd chebyshev1 filtered signal
f_2side=(0:N/2)*df';
y_f=fft(s_filt_cheb_2);
p2_f= abs(y_f/N);
p1_f=p2_f(1:N/2+1,:);
p1_f(2:end-1,:)=2*p1_f(2:end-1,:);
figure(); sgtitle("2nd chebysheb1 vs noisy in frequency domain")
loglog(f_2side,[p1(:,2),p1_f])
legend("all noise","filtered")
ylim([10e-4,10]);
xlabel('Hz')
%% time plot of 1st butterworth filtered signal
figure(); sgtitle('1st butterworth filtered signal vs noisy vs pure signal in time domain ')
plot(t_vec,[s_all,s_filt_butter_1,s0]);
legend("noisy","filtered","pure")
xlabel('T[sec]')
%% time plot of 2nd butterworth filtered signal
figure(); sgtitle('2nd butterworth filtered signal vs noisy vs pure signal in time domain ')
plot(t_vec,[s_all,s_filt_butter_2,s0]);
legend("noisy","filtered","pure")
xlabel('T[sec]')

%% time plot of 1st chebyshev filtered signal
figure(); sgtitle('1st chebyshev 1 filtered signal vs noisy vs pure signal in time domain ')
plot(t_vec,[s_all,s_filt_cheb_1,s0]);
legend("noisy","filtered","pure")
xlabel('T[sec]')
%% time plot of 2nd chebyshev filtered signal
figure(); sgtitle('2nd chebyshev 1 filtered signal vs noisy vs pure signal in time domain ')
plot(t_vec,[s_all,s_filt_cheb_2,s0]);
legend("noisy","filtered","pure")
xlabel('T[sec]')
%% frequency response of the filters
figure(); sgtitle("frequency response of the 1st lowpass butterworth filter")
freqs(b_butter_l_1,a_butter_l_1)
figure(); sgtitle("frequency response of the 2nd lowpass butterworth filter")
freqs(b_butter_l_2,a_butter_l_2)
figure(); sgtitle("frequency response of the 2nd lowpass chebyshev1 filter")
freqs(b_cheb_l_2,a_cheb_l_2)
figure(); sgtitle("frequency response of the 1st lowpass chebyshev1 filter")
freqs(b_cheb_l_1,a_cheb_l_1)
%% saving all in a matrix
signals_and_fil=[signals,s_filt_butter_1,s_filt_butter_2,s_filt_cheb_1,s_filt_cheb_2];
if save_data
writematrix(signals_and_fil,'all signals.txt')
end