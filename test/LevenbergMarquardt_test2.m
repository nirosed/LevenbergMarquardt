% 本测试脚本使用随机生成的带高斯白噪声的虚拟信号
% 在寻找到峰值大致位置后进行基于LM算法的高斯拟合
% 在最后给出误差
clear;clc;
tic;
%% SignalGenerating
t_start = 0;
t_end = 80000;
index = linspace(1,t_end,t_end)';

noise_mu = 0.35;
noise_sigma = 0.0154;

peak_position = [
    12345 + rand();
    23456 + rand();
    67890 + rand()
    ];
amplitude = [
    1 + 0.2 * rand() * sign( rand() - 0.5 );
    1 + 0.2 * rand() * sign( rand() - 0.5 );
    1 + 0.2 * rand() * sign( rand() - 0.5 )
    ];
width_pulse = 20+rand();

[ noise , signal ] = SignalGenerator( index,noise_mu,noise_sigma,peak_position,width_pulse,amplitude );
signal = signal.^2 + noise;
disp(['width_pulse: ',num2str(width_pulse),', dc: ',num2str(noise_mu)]);
disp(['peak_position: ',num2str(peak_position(1)),', ',num2str(peak_position(2)),', ',num2str(peak_position(3)),'.']);
disp(['amplitude: ',num2str(amplitude(1)),', ',num2str(amplitude(2)),', ',num2str(amplitude(3)),'.']);
%% peak search
threshold = 0.5;
n_extremum = 3;
peak_or_valley = 0; % peak（0） or valley（1）
extrenum_distance = 10;
[ locs , nFound ] = ExtremumSearch( index , signal , n_extremum, peak_or_valley , threshold , extrenum_distance );
%% plot peak
% figure(1),plot( index , signal , index(locs) , signal(locs) , '*r' );
%% setup
% 迭代初值
width = 20;
dc = 0.05;
% 误差设定
epsilon_1 = 1e-10;
epsilon_2 = 1e-20;
tau = 1e-6;
% 最大迭代次数
N = 50;
% 拟合单边取得点数
M = 5;
%% LM
output = '';
locs_error = NaN * zeros(3,1);
amplitude_locs = NaN * zeros(3,1);
width_error = NaN * zeros(3,1);
dc_error = NaN * zeros(3,1);
for ii = 1 : nFound
    x0 = [signal(locs(ii));locs(ii);width;dc];
    t = locs(ii) + linspace( -M , M , 2*M+1 )';
    y = signal( t(1) : t(2*M+1) );
    
%     tt = linspace(t(1),t(length(t)),100);
%     figure(2),plot(t,y,'*',tt,x0(1) * exp( -( ( tt - x0(2) ) / x0(3) ).^2 / 2 ) + x0(4),'r');
%     title(['peak #',num2str(ii),', iterative:',num2str(0)]);
%     pause(0.2);

    x = LevenbergMarquardt( t , y , x0 , epsilon_1 , epsilon_2 , tau , N );
    output = [ 'A: ' , num2str( x(1) ) , ', mu(position): ' , num2str( x(2) ) ,...
        ', sigma(width): ' , num2str( x(3) ) , ', c(dc): ' , num2str( x(4) ) ];
    disp(output);
    locs_error(ii,1) = abs( ( x(2) - peak_position(ii) ) / peak_position(ii) );
    amplitude_locs(ii,1) = abs( ( x(1) - amplitude(ii) ) / amplitude(ii) );
    width_error(ii,1) = abs( ( x(3) - width_pulse ) / width_pulse );
    dc_error(ii,1) = abs( ( x(4) - noise_mu ) / noise_mu );
end
%% error about positions
disp(['locs_error: ',num2str(locs_error(1)*1000),'‰, ',num2str(locs_error(2)*1000),'‰, ',num2str(locs_error(3)*1000),'‰.']);
% disp(['amplitude_locs: ',num2str(amplitude_locs(1)*100),'%, ',num2str(amplitude_locs(2)*100),'%, ',num2str(amplitude_locs(3)*100),'%.']);
% disp(['width_error: ',num2str(width_error(1)*100),'%, ',num2str(width_error(2)*100),'%, ',num2str(width_error(3)*100),'%.']);
% disp(['dc_error: ',num2str(dc_error(1)*100),'%, ',num2str(dc_error(2)*100),'%, ',num2str(dc_error(3)*100),'%.']);
toc;
