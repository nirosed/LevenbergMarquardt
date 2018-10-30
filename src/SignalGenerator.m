function [ noise , signal ] = SignalGenerator( t_vector,mu,sigma,peak_position,width_pulse,amplitude )
%% signal_with_noise_generator
% normal distribution noise is generated in Box-Muller method

%% noise generation and display

N = length(t_vector);
noise = NaN*zeros(N,1);
for i = 1 : N
    noise(i,1) = NormalDistribution( 1 );
end
noise = noise * sigma + mu;
% disp(['mu = ',num2str(mean(noise))]);
% disp(['sigma^2 = ',num2str(var(noise))]);
%% signal generation
if ( min(peak_position) <= t_vector(1) ) || ( max(peak_position) >= t_vector(end) )
    error('peak_position must between t_start and t_end!');
end
signal = zeros(N,1);
for ii = 1 : length(peak_position)
    signal = signal + amplitude(ii) * exp( -( t_vector - peak_position(ii) ).^2 / 2 / width_pulse^2 );
end

end
