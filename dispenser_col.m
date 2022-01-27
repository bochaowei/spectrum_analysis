filename='F0002';
filename1=[filename,'CH1.CSV'];
filename2=[filename,'Ch2.CSV'];
%sat=csvread(filename2,0,4,[0,4,2499,4]);
sat=readmatrix(filename2);
sat=sat(:,5);
%signal=csvread(filename1,0,4,[0,4,2499,4]);
signal=readmatrix(filename1);
signal=signal(:,5);
%signals=sgolayfilt(signal,9,65);
%figure(3)
%plot([signal,sat])
x=1:length(sat);
plot(x,sat,x,signal)
%plot(x,signal(x),x,sat)
%%
start=181; %%set a start point to find the peak to avoid wrong peaks
stop=2250;
sats=sgolayfilt(sat,10,85);
%sats=sat;
%figure(1)
minpeaks=100;
findpeaks(sats(start:stop),'MinPeakDistance',minpeaks);
[pks,locs]=findpeaks(sats(start:stop),'MinPeakDistance',minpeaks);
locs=locs+start; %%%the locs found here is without the start columns, need to add back
loc_0=5;
loc_1=4;
scale=133.3/(locs(loc_0)-locs(loc_1));
range=900;
%x=max(600,locs(1)-range):min(locs(1)+range+1050,2499);
x=start:stop;
frequency=transpose((x-locs(loc_0))*scale);
figure(2)
mysignal=signal(x);
sats=sats(x);
sats=(sats-min(sats));
sats=sats/max(sats);
mysignal=(mysignal-min(mysignal));
mysignal=mysignal/max(mysignal);
plot(frequency,mysignal,frequency,sats)
xlabel('Frequency (MHz)')
ylabel('Signal (normalized)')
%% longitudinal speed
lv=frequency/1.282/8*17;
plot(lv,mysignal,frequency,sats)
xlabel('Speed(m/s)')
ylabel('Signal (normalized)')
%%
range=500:1700;
lv_s=lv(range);
height=mysignal(range);
plot(lv_s,height)
ave_v=sum_c(lv_s.*height)/sum_c(height); %%%the dv cancel out in numerator and demoninator
%% gaussian fit
 %%%%%fit the function
 fitfun = fittype(@(a1,b1,c1,d,x) a1*exp(-((x-b1)/c1).^2)+d,'independent',{'x'});
 
 options = fitoptions(fitfun);
 options.StartPoint = [1,0,80,0];
 %options.Lower = [max(mysignal)-0.2,-50,0,0];
 options.Upper = [1.2,9,1000,0.2];
 fitdata = fit(frequency,mysignal,fitfun,options);
 figure(2)
 plot(fitdata,frequency,mysignal)
 xlabel('Frequency(MHz)')
 ylabel('Voltage(V)')
 

b1=fitdata.b1
c=fitdata.c1;
fwhm=c*2*sqrt(log(2))
title('fwhm is '+string(fwhm)+' MHz')
saveas(figure(2),'figure')

%% lorentz fit
fitfun = fittype(@(a1,b1,c1,d,a2,x) a1*1./(1+(2.*(x-b1)/c1).^2)+a2*x+d,'independent',{'x'});
options = fitoptions(fitfun);
figure(2)
options.StartPoint = [1,0,100,0,0.1/800];
options.Lower = [0.98,-100,0,-0.1,-0.1];
options.Upper = [1.4,40,500,0.3,0.5];
fitdata = fit(frequency,mysignal,fitfun,options);
hold on
box on
plot(fitdata,frequency,mysignal)
plot(frequency,sats,'linewidth',2)
legend('data','fitted curve','sat spectrum')
hold off
xlabel('Frequency(MHz)')
c=fitdata.c1;
fwhm=c
title('fwhm is '+string(fwhm)+' MHz')
saveas(figure(2),'Lorentz fit')
%fitdata
%%
figure(10)
hold on
%plot(frequency,fitdata.a1*1./(1+(2.*(frequency-fitdata.b1)/fitdata.c1).^2),'linewidth',2)
peak_s=mysignal-fitdata.a2*frequency-fitdata.d;
plot(frequency,peak_s,'linewidth',1)
plot(frequency,sats,'linewidth',2)
%legend('Lorentz fit','Signal')
xlabel('Relative frequency (MHz)')
ylabel('Signal (arb. unit)')
xlim([-600,480])
legend('FWHM=220 MHz')
hold off
ax = gca;
ax.FontSize = 16;
ax.FontWeight='bold';
ax.LineWidth = 1;
box on
%%
peak_s=mysignal-fitdata.a2*frequency-fitdata.d;
plot(frequency,peak_s,frequency,sgolayfilt(peak_s,10,405))
%% data loading for lifetime measurement from picoscope
data=load('5.5A.mat');  %%%load the data from lifetime measurement
peak_s=data.A;
sat=data.B;
cut=4300:8201;
cut=700:4000; %% for 5.5A
peak_s=peak_s(cut);
sat=sat(cut);
c23=6665;%%% for 7.5A
c23=2269; %for 5.5A
cross=6382;%%% for 7.5A
cross=2581;
scale=133.3/(c23-cross);
frequency=(cut-c23)*scale;
sat=sat-min(sat);
sat=sat/max(sat);
peak_s=peak_s-min(peak_s);
peak_s=peak_s/max(peak_s);
plot(frequency,sat,frequency,peak_s)
%%
Pprobe=2; %mw
w1=0.586; %mm
w2=0.586;
Isat=3.05*10^(-2); %mw/mm^2%
Iprobe=Pprobe/(pi*w1*w2);
sp=Iprobe/Isat;

peak_s=sgolayfilt(peak_s,10,405);
%%
hold on 
gamma=6.066;
frequencystep=(frequency(length(frequency))-frequency(1))/(length(frequency)-1);
Rsc=1./(1+sp+4*(frequency./gamma).^2);
%[u,w]=deconvolution(peak_s,Rsc); %% for oscillopcope data
[u,w]=deconvolution(peak_s,Rsc'); %%%this for picoscope data
sum_c=0;
for i=1:length(frequency)
sum_c=sum_c+w(i)*abs(frequencystep)/1.282;
end
w=w/sum_c;
v=u.*frequencystep./1.282;
plot(v,w,'linewidth',2)
ww=w;
%center=1037; %% for oscilloscope data
center=1950; %% for dispenser lifetime picoscope data
center=2500;
for ii=1:length(v)
    if v(ii)<0 && (2*center-ii<length(v))
        ww(ii)=ww(2*center-ii);
    end
end
plot(v,ww,'linewidth',1.5)
ax = gca;
ax.FontSize = 16;
ax.FontWeight='bold';
ax.LineWidth = 1;
box on
ylabel('Ratio of atoms')
xlabel('Velocity (m/s)')
legend('FWHM=120 m/s','Symetric')
ylim([-1*10^-3,6.5*10^-3])
%yticks([])
x=[-500,-375,-250,-125,0,125,250,375,500];
xticks(x)
S = string(x);
S(mod(x, 2)~=0) = "";
xticklabels(S)

hold off

detune=v*1.282;
tau=38.117*10^6;
R_average=0;
for i=1:length(frequency)
%R_average=R_average+w(i)*tau/2*sp/(1+sp+4*(detune(i)./gamma).^2)*abs(frequencystep);
R_average=R_average+ww(i)*tau/2*sp/(1+sp+4*(detune(i)./gamma).^2)*abs(frequencystep);%% lifetime high current
end


%% get the angular distribution from transverse velocity distribution and
%%% a longitudinal Maxwell distribution, by hand fitting

d_theta=0.001;
theta=0.0001:d_theta:pi/2;
sigma=0.22;
angular=1/sigma/sqrt(2*pi)*exp(-0.5*(theta/sigma).^2);
plot(theta,angular)

%velocity=-900:1:900;
v_trans=-500:500;
%Fv=2/mu*(velocity/mu).^3.*exp(-(velocity/mu).^2);
p=[];
for i=v_trans
    sum_c=0;
    for j=theta
        sum_c=sum_c+inter(i,j)*d_theta;
    end
    p=[p,sum_c];
end
nor=0;
for i=1:length(v_trans)
    nor=nor+p(i);
end
p=p/nor*0.92;
plot(v_trans,p,v,w)

%% calculate the ratio of atoms enter the probe region, used x y z integration
sigma=0.22;
d=6*10^-3;
w1=0.586*10^-3; %mm
w2=0.586*10^-3;
L=6*10^-3;
fun = @(x,y) 1/sigma/sqrt(2*pi).*exp(-0.5*(atan(sqrt(x.^2+y.^2)/d)/sigma).^2).*d./(x.^2+y.^2+d^2).^(3/2);
ratio=integral2(fun,-L-0.68*10^-3,-0.68*10^-3,-w1,w1);
%%%normalization
d_theta=0.001;
theta=0.0001:d_theta:pi/2;
sum2=0;
for x=theta
sum2=sum2+1/sigma/sqrt(2*pi)*exp(-0.5*(x/sigma).^2).*sin(x)*d_theta;
end
sum2=sum2*2*pi;
ratio_f=ratio/sum2

%% calculate the atoms ratio within a certain velocity range
sum_c=0;
for i=length(frequency)-1
    if abs(v(i))<=5
        sum_c=sum_c+w(i)*(v(i+1)-v(i));
    end
end
sum_c
%%
function val=inter(vtran,theta)
mu=302.1;
sigma=0.22;
ftheta=1/sigma/sqrt(2*pi)*exp(-0.5*(theta/sigma).^2);
velocity=vtran/sin(theta);
Fv=2/mu*(velocity/mu).^3.*exp(-(velocity/mu).^2);
val=ftheta*Fv/(vtran+0.00001);
end
