close all;
%%
y = app.MyRun.MATDSS.Cont.CA(2).P0Set(2:end);
Ts = 0.1
Fs = 1/Ts;
N = length(y);
y = [y, zeros(1, N - length(y))]

f=-Fs/2:Fs/(N-1):Fs/2;                      %Frequency Vector to plot

dt = 1/Fs;

t = 0:dt:(N*dt)-dt;

h = figure(1);
title('')
subplot(3,1,1)
plot(t, y);


Y = fftshift(fft(y));

subplot(3,1,2)
plot(f, abs(Y))


subplot(3,1,3)
plot(f,angle(Y)*180/pi)

%%

y = app.MyRun.MATDSS.Cont.CA(2).P0Set(2:end);
Ts = 0.1
Fs = 1/Ts;
N = 1024;
y = [y, zeros(1, N - length(y))]

f=-Fs/2:Fs/(N-1):Fs/2;                      %Frequency Vector to plot

dt = 1/Fs;

t = 0:dt:(N*dt)-dt;

h = figure(2);
subplot(3,1,1)
plot(t, y);


Y = fftshift(fft(y));

subplot(3,1,2)
plot(f, abs(Y))


subplot(3,1,3)
plot(f,angle(Y)*180/pi)

%%

lpFilt=designfilt('Lowpassfir', ...       % Response type
'StopbandFrequency',1, ...     % Frequency constraints
'PassbandFrequency',0.5, ...
'StopbandAttenuation',80, ...    % Magnitude constraints
'PassbandRipple',4, ...
'DesignMethod','kaiserwin', ...  % Design method
'ScalePassband',false, ...       % Design method options
'SampleRate',10);               % Sample rate

% fvtool(lpFilt)

%%

y = app.MyRun.MATDSS.Cont.CA(2).P0Set(2:end);

Filteredmessage = filter(lpFilt,y);


Ts = 0.1
Fs = 1/Ts;
N = length(y);
y = [y, zeros(1, N - length(y))]

f=-Fs/2:Fs/(N-1):Fs/2;                      %Frequency Vector to plot

dt = 1/Fs;

t = 0:dt:(N*dt)-dt;
close all
h = figure(3);
plot(t, y);
hold on
plot(t, Filteredmessage);

legend('y','filtered')

%%

Tc = 0.1:0.1:3;
fc = 1./Tc;
Ts = 0.1;
wc = 2*pi*fc;
for i = 1 : length(wc)
    Hs(i) = {tf([1], [1./wc(i), 1])};
end


% freqs([1],[1./wc(end), 1],)



%% Filtering in Frequency domain
close all

% y_original = app.MyRun.MATDSS.Cont.CA(2).P0Set(2:end);
% Ts = 0.1;
% % Here, we will do testing on different LPF configs. 

%First, get a copy of the original signal without any filtering
y_original = app.MyRun.MATDSS.Cont.CA(2).P0Set(2:end);
Ts = 0.1;
Fs = 1./Ts;
N = 1024;
f=-Fs/2:Fs/(N-1):Fs/2;                      %Frequency Vector to plot
dt = 1./Fs;
t = 0:dt:(N*dt)-dt;
y = [y_original, zeros(1, N - length(y_original))];

Tc = 3;
fc = 1/Tc;

Hs = tf([1],[1./(2*pi*fc), 1]);

Hz = c2d(Hs, Ts, 'tustin');

[h,w] = freqz(Hz.Numerator{1}, Hz.Denominator{1}, 2^10);

figure(1)
subplot(3,1,1)
hold on
title('Frequency Response')
plot(Fs*w/(2*pi),abs(h))
xlabel('f(Hz)')
subplot(3,1,2)
hold on
title('Frequency Response (dB)')
plot(Fs*w/(2*pi),db(h))
set(gca,'YScale','log')
ylim([-200 -1])
grid on
xlabel('f(Hz)')
subplot(3,1,3)
hold on
title('Phase Response')
plot(Fs*w/(2*pi), angle(h)*180/pi)
xlabel('f(Hz)')




%Filtering the signal
y_filtered = filter(Hz.Numerator{:}, Hz.Denominator{:}, y_original);

figure(2)
plot(0:Ts:(length(y_original)-1)*Ts, y_original);
hold on
plot(0:Ts:(length(y_original)-1)*Ts, y_filtered);


%% Adding MATDSS fitlered signal


%%



Tc = [0.2, 0.3, 0.5, 1, 2];
for i = 1:length(Tc)
    h = figure(1);
    
    % Filter z-TF
    Hs = tf([1],[Tc(i)./(2*pi), 1]);
    Hz = c2d(Hs,Ts,'tustin');
    [h,w] = freqz(Hz.Numerator{:}, Hz.Denominator{:},1024);
    
    % figure
    % plot(w/pi, abs(h))

    
    subplot(length(Tc),3,1+(i-1)*3)
    plot(t(1:length(y_original)), y_original);
    % filtfilt
    
    subplot(length(Tc),3,2+(i-1)*3)
    Y = fftshift(fft(y));
    plot(f, abs(Y))
    
    
    
    YFiltered = Y.*IR;
    subplot(length(Tc),3,1+(i-1)*3)
    yfiltered = real(fftshift(ifft(YFiltered)));
    hold on
    plot(t(1:length(y_original)), yfiltered(1:length(y_original)));
    legend('y', 'yfiltered')
    xlim([0,10])

    subplot(length(Tc),3,3+(i-1)*3)
    plot(f,angle(Y)*180/pi)

end