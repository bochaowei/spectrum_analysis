% pd_current_plot.m
% Code for atom dispenser output calibration using large data sets
% parts of analysis based off code by Chao Li / Xiao Chai

% Alex Crawford
% 12/8/2021

% Function for isolating 2->3' transition and extracting PD data (Alex):
% peakPointData.m

% code sections:
% 1. Variable min denomination for disp lifetime data
%       Run this section first ... it pulls data from the external hard
%       drive (hard coded path needs adjustment per user) and processes the
%       spectra to isolate the F=2 -> F'=3 transition data. The specific
%       layout of the files/ file names are hard coded as well.
%       There are divisions to this first section:
%       A. Experimental Parameters
%       B. Set data density --> specify hour denominations to data
%       C. Initialize data arrays and time-stamp
%       D. Data analysis - main folder G=10^9
%       E. Data analysis - main folder G=10^8
%       F. Data analysis - subfolder G=10^8
%       G. Noise floor
% 2. Generate throughput vs time plot and estimated Rb emitted
%       Run this section second ... it generates a plot of the throughput
%       vs time based on the previous section's data. Do not clear data in
%       memory before running, graph will not render.
% 3. test plot for data jump
%       Irrelevant section.
% 4. test plot for all spectra in each separate folder
%       Irrelevant section.


%% Variable min denomination for disp lifetime data

clear all;
close all
clc;

% ----------------------------------------------------------
% Experimental parameters

tt = [204]; % oven body temperature
T = tt+273.15; % oven body temperature in kelvin
kb = 1.38064852e-23; % boltzmann constant
m = 1.4431606011e-25; % mass of the 87Rb atom
vbeam = sqrt(3.*kb./m.*T); % beam most probably speed
hbar=1.0545718e-34; % planck constantmax

omega0 = 2*pi*384.230484468562*1e12; % laser frequency
windowtrans = 0.97; % transmission of the vacuum window
lenstrans = 0.927; % transmission of the lens set
sldang = 2*pi*(1-cos(14/180*pi)); % LCP01 lens mount; the spanned half angle is measured to be 14 degrees
dpsldang = 8*pi/3; % total solid angle formed by a dipole emission
intFilt = 63.9/100; % interference filter transmission
eta = sldang/dpsldang*lenstrans*windowtrans*intFilt; % efficiency

pofenter = 0.11; % probability for atoms entering the measuring region

resp = 0.6; % Si photodiode responsitivity at 780 nm
% gain of amplifier (volts/amp) found in each sub-section below

Ab = 0.278; % isotope abundance
hyperf = 5/8; % considering the degeneracy of hyperfine levels

% probe_exit = 2.8; % cm distance between the probe location and the collimator exit


% laser intensity and saturation parameter
power = 2; %mW
w_x = 0.0586; % cm  
w_y = 0.0586; % cm
intensity = 2*power/(pi*w_x*w_y); % 77 mW/cm^2
Isat = 3.05; % mW/cm^2
s_p = intensity/Isat*0.97; % saturation parameter is about 23

ss = s_p;
gamma = 38.117*1e6; % !!! scattering rate should have the factor of 2pi in it!!!
lambda = 780.2412e-9; % wavelength of laser
% rsc_vt = @(vt_r) gamma./2.*ss./(1+ss+4.*(vt_r./lambda./(gamma./(2.*pi))).^2); % scattering rate function
% rsc_vt = gamma./2.*ss./(1+ss); % approx scattering rate function zero detuning


% DIRECTLY FROM BOCHAO CODE OUTPUT FOR THIS EXPERIMENT
R_average = 6.6410e+06;
R_avg_2 = 2.9*10^6; % for currents higher than 8A ... broader distribution


% -----------------------------------------------------------
% Set data density 

% 10 min denomination =6 ... use =1 to sample every hour
hour_denom = 60;


% ----------------------------------------------------------
% Initialize data arrays and time-stamp

curr_array = [];
index_array = [];
ind_main = [];
ind_sub = [];
noise_array = [];
time_stamp_arr = [];

promVal = 0.015;


% first data file for relative start day
dataFileName = sprintf('/Volumes/USBDisk/lifetime/20211206-0001 (50).mat');
D=dir(dataFileName);
TimeStamp = {D.date};
t=string(TimeStamp);
c=textscan(t,'%f-%s %f:%f:%f');
initFileTime = (24*(c{1}-6))+c{3}+(c{4}/60);


% ----------------------------------------------------------
% Data analysis - main folder G=10^9


for k0 = 50:(1005/hour_denom):48131
    gain = 10^9;
    dataFileName = sprintf('/Volumes/USBDisk/lifetime/20211206-0001 (%d).mat',round(k0));
    % for matlab 2017a and before use exist 
    if exist(dataFileName, 'file') == 2

        D=dir(dataFileName);
        TimeStamp = {D.date};
        t=string(TimeStamp);
        c=textscan(t,'%f-%s %f:%f:%f');
        fileTime = ((24*(c{1}-6))+c{3}+(c{4}/60)) - initFileTime;
        time_stamp_arr = [time_stamp_arr,fileTime];

        dataLoad = load(dataFileName,"A","B","Length","Tinterval","Tstart");
        ChA = dataLoad.A;
        ChB = dataLoad.B;
        dt = dataLoad.Tinterval;
        start = dataLoad.Tstart;
        len = dataLoad.Length;
        testTime = start+(0:len-1)*dt;

        flat = find(ChB == max(ChB),1);
        isolate = ChB(flat-50:flat+50);
        stdEr = std(isolate);
        noise_array = [noise_array, stdEr];

        cut = find(testTime > -0.74,1);

        % finding and extracting data peaks
        % function has plot for troubleshooting ... uncomment if wanted
        [dataPt, chooseTime, chooseA] = peakPointData(ChA(cut:end),ChB(cut:end),testTime(cut:end), promVal);

        V_to_measured_flow = (hbar*omega0*resp*gain*eta*Ab*hyperf.*(pi.*(w_x/100)./2./vbeam).*R_average).^-1;
        measured_flow = dataPt*V_to_measured_flow;
        throughput = measured_flow./pofenter; % atoms/sec

        curr_array = [curr_array,throughput];
        index_array = [index_array,k0];
        ind_main = [ind_main,k0];

    else 
        disp('error loading data file:');
        disp(dataFileName);
    end

end 
disp(length(curr_array));
disp(length(time_stamp_arr));


% ------------------------------------------------------------------
% Data analysis - main folder G = 10^8


for k1 = 48132:(1005/hour_denom):48333
    gain = 10^8;
    dataFileName = sprintf('/Volumes/USBDisk/lifetime/20211206-0001 (%d).mat',round(k1));
    % for matlab 2017a and before use exist 
    if exist(dataFileName, 'file') == 2

        D=dir(dataFileName);
        TimeStamp = {D.date};
        t=string(TimeStamp);
        c=textscan(t,'%f-%s %f:%f:%f');
        fileTime = ((24*(c{1}-6))+c{3}+(c{4}/60)) - initFileTime;
        time_stamp_arr = [time_stamp_arr,fileTime];

        dataLoad = load(dataFileName,"A","B","Length","Tinterval","Tstart");
        ChA = dataLoad.A;
        ChB = dataLoad.B;
        dt = dataLoad.Tinterval;
        start = dataLoad.Tstart;
        len = dataLoad.Length;
        testTime = start+(0:len-1)*dt;

        flat = find(ChB == max(ChB),1);
        isolate = ChB(flat-50:flat+50);
        stdEr = std(isolate);
        noise_array = [noise_array, stdEr];

        cut = find(testTime > -0.74,1);

        % finding and extracting data peaks
        % function has plot for troubleshooting ... uncomment if wanted
        [dataPt, chooseTime, chooseA] = peakPointData(ChA(cut:end),ChB(cut:end),testTime(cut:end), promVal);
        
        V_to_measured_flow = (hbar*omega0*resp*gain*eta*Ab*hyperf.*(pi.*(w_x/100)./2./vbeam).*R_avg_2).^-1;
        measured_flow = dataPt*V_to_measured_flow;
        throughput = measured_flow./pofenter; % atoms/sec

        curr_array = [curr_array,throughput];
        index_array = [index_array,k1];
        ind_main = [ind_main,k1];

    else 
        disp('error loading data file:');
        disp(dataFileName);
    end

end 
disp(length(curr_array));
disp(length(time_stamp_arr));


% ------------------------------------------------------------------
% Data analysis - subfolder G=10^8


for k2 = 2:(1200/hour_denom):31763
    gain = 10^8;
    dataFileName = sprintf('/Volumes/USBDisk/lifetime/7.5/20211206-0001 (%d).mat',round(k2));
    % for matlab 2017a and before use exist 
    if exist(dataFileName, 'file') == 2

        D=dir(dataFileName);
        TimeStamp = {D.date};
        t=string(TimeStamp);
        c=textscan(t,'%f-%s %f:%f:%f');
        fileTime = ((24*(c{1}-6))+c{3}+(c{4}/60)) - initFileTime;
        time_stamp_arr = [time_stamp_arr,fileTime];

        dataLoad = load(dataFileName,"A","B","Length","Tinterval","Tstart");
        ChA = dataLoad.A;
        ChB = dataLoad.B;
        dt = dataLoad.Tinterval;
        start = dataLoad.Tstart;
        len = dataLoad.Length;
        testTime = start+(0:len-1)*dt;

        flat = find(ChB == max(ChB),1);
        isolate = ChB(flat-50:flat+50);
        stdEr = std(isolate);
        noise_array = [noise_array, stdEr];

        % finding and extracting data peaks
        % function has plot for troubleshooting ... uncomment if wanted
        [dataPt, chooseTime, chooseA] = peakPointData(ChA(cut:end),ChB(cut:end),testTime(cut:end), promVal);

        V_to_measured_flow = (hbar*omega0*resp*gain*eta*Ab*hyperf.*(pi.*(w_x/100)./2./vbeam).*R_avg_2).^-1;
        measured_flow = dataPt*V_to_measured_flow;
        throughput = measured_flow./pofenter; % atoms/sec

        curr_array = [curr_array,throughput];
        index_array = [index_array,k2];
        ind_sub = [ind_sub,k2];

    else
        disp('error loading data file:');
        disp(dataFileName);
    end
end 
disp(length(curr_array));
disp(length(time_stamp_arr));


% -----------------------------------------------------------------
% Noise floor 

noise = max(noise_array); 
V_to_measured_flow_N = (hbar*omega0*resp*10^9*eta*Ab*hyperf.*(pi.*(w_x/100)./2./vbeam).*R_average).^-1;
measured_flow_N = noise*V_to_measured_flow_N;
throughput_N = measured_flow_N./pofenter; % atoms/sec

% -----------------------------------------------------------------

%% Generate throughput vs time plot and estimated Rb emitted

figure(1)
semilogy(time_stamp_arr,curr_array, 'b.','MarkerSize',10);
% yline(throughput_N,'r--','linewidth',3);
ax = gca;
ax.FontSize = 16;
ax.FontWeight='bold';
ax.LineWidth = 1;
box on
ylabel('Throughput (atoms/s)');
xlabel('Time (hour)');
% legend('Signal', 'Noise Floor');
xlim([0,80]);
ylim([10^10, 10^14]);


% integrate under curve and estimate emitted Rb in mg
integ = 0;
curr_nan = (curr_array(~isnan(curr_array)));
for ii = 2:length(curr_nan)
    integ = integ + ((curr_nan(ii-1)+curr_nan(ii))/2*(60*60/hour_denom));
end
disp(integ);

atom_mg = 10^3*integ/(6.022*10^23)*85.4678;
disp(atom_mg);



%% ----------------------------------------------------------------
%% test plot for data jump

% function peakPointData_sub used in this test section

% temp arrays
ptm_arr = [];
ptam_arr = [];
F2_arr = [];

% disp(ind_main(end-12:end));

figure(2)
% test plot to check how many data points at end of main folder need G=10^8
for kk = ind_main(end-12:end)
    dataFileName = sprintf('/Volumes/USBDisk/lifetime/20211206-0001 (%d).mat',round(kk));
    dataLoad = load(dataFileName,"A","B","Length","Tinterval","Tstart");
    ChA = dataLoad.A;
    ChB = dataLoad.B;
    dt = dataLoad.Tinterval;
    start = dataLoad.Tstart;
    len = dataLoad.Length;
    testTime = start+(0:len-1)*dt;

    flat = find(ChB == max(ChB),1);
    isolate = ChB(flat-50:flat+50);
    stdEr = std(isolate);
    noise_array = [noise_array, stdEr];

    cut = find(testTime > -0.74,1);

    % finding and extracting data peaks
    [dataPt, chooseTime, chooseA, ptm, ptam, F2] = peakPointData_sub(ChA(cut:end),ChB(cut:end),testTime(cut:end), promVal);
    ptm_arr = [ptm_arr, ptm];
    ptam_arr = [ptam_arr, transpose(ptam)];
    F2_arr = [F2_arr, F2];

    hold on
    plot(testTime,ChA,ptm(F2),ptam(F2),'ro');
    plot(testTime,ChB,'--');
    hold off

end


%% ----------------------------------------
%% test plot for all spectra in each separate folder

% refer to test notebook disp_output_calib.mlx for other test plots

%  Test plotting drift ... irrelevant as of 1/2/2022
% takes main and sub index arrays from above

figure(3)
% reducing to sample = denom/4
samp = hour_denom/4;

% test plot everything from main folder
for j = 1:samp:length(ind_main)
    testFile = sprintf('/Volumes/USBDisk/lifetime/20211206-0001 (%d).mat',round(ind_main(j)));
    testLoad = load(testFile,"A","B","Length","Tinterval","Tstart"); 

    testChA = testLoad.A;
    testChB = testLoad.B;
    dt = testLoad.Tinterval;
    start = testLoad.Tstart;
    len = testLoad.Length;
    testTime = start+(0:len-1)*dt;
    [TF1,P] = islocalmax(testChB, 'MinProminence',promVal);
    abs_max = max(ChB);
    prom_ch = testChB(TF1);
    prom_t = testTime(TF1);
    find_maxL = find(prom_ch == abs_max,1);
    F2_ind_p = find_maxL - 3;

    hold on
%     plot(testTime,testChB,prom_t,prom_ch,'r*',prom_t(F2_ind_p),prom_ch(F2_ind_p),'bo');
    plot(testTime,testChB,prom_t(F2_ind_p),prom_ch(F2_ind_p),'ro');
end
title('main folder test plot');
hold off


% same for sub-folder
figure(4)
for js = 1:samp:length(ind_sub)
    testFile = sprintf('/Volumes/USBDisk/lifetime/7.5/20211206-0001 (%d).mat',round(ind_sub(js)));
    testLoad = load(testFile,"A","B","Length","Tinterval","Tstart"); 

    testChA = testLoad.A;
    testChB = testLoad.B;
    dt = testLoad.Tinterval;
    start = testLoad.Tstart;
    len = testLoad.Length;
    testTime = start+(0:len-1)*dt;
    [TF1s,Ps] = islocalmax(testChB, 'MinProminence',promVal);    
    abs_max = max(ChB);
    prom_chs = testChB(TF1s);
    prom_ts = testTime(TF1s);
    find_maxRs = find(prom_chs == abs_max,1);
    F2_ind_ps = find_maxRs + 3;

    hold on
%     plot(testTime,testChB,prom_ts,prom_chs,'r*',prom_ts(F2_ind_ps),prom_chs(F2_ind_ps),'bo');
    plot(testTime,testChB,prom_ts(F2_ind_ps),prom_chs(F2_ind_ps),'ro');

end
title('sub folder test plot');
hold off


