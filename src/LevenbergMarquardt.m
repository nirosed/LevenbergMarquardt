function [ x ] = LevenbergMarquardt( t , y , x0 , epsilon_1 , epsilon_2 , tau , N )
%% Levenberg-Marquardt算法，进行高斯曲线拟合
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入变量为 t , y , x0 , epsilon_1 , epsilon_2 , tau , N
% t：自变量
% y：因变量
% x0：参数初始猜测值
% epsilon_1：迭代终止误差，用户定义
% epsilon_2：迭代终止误差，用户定义
% tau：倍率，用户定义。初值接近极值时为1e-6，否则可设定为1或1e-3
% N：迭代最高次数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输出变量为 x
% x：迭代终止后参数的结果
if length(x0) ~= 4
    error('error x length!');
end
num = min(t) - 1;
t = t - num;
x0(2) = x0(2) - num;
%% Objective function vector
fun = @(x) y - ( x(1) * exp( -( ( t - x(2) ) / x(3) ).^2 / 2 ) + x(4) ); % gauss type
% fun = @(x) log( abs( y - x(4) ) ) - log( abs( x(1) ) ) + ( ( t - x(2) ) / x(3) ).^2 / 2; % log type
%% Jacobi matrix for gauss function
Jacobi = @(x) [
    -exp( -( ( t - x(2) ) / x(3) ).^2 / 2 ) ,...
    -x(1) * exp( -( ( t - x(2) ) / x(3) ).^2 / 2 ) .* ( t - x(2) ) / x(3)^2 ,...
    -x(1) * exp( -( ( t - x(2) ) / x(3) ).^2 / 2 ) .* ( t - x(2) ).^2 / x(3)^3 ,...
    -ones( length(t) , 1 )
    ]; % gauss type
% Jacobi = @(x) -1 * [
%     ones( length(t) , 1 ) / x(1) ,...
%     ( x(2) - t ) / x(3)^2,...
%     -( t - x(2) ).^2 / x(3)^3,...
%     1 ./ ( x(4) - y )
%     ]; % log type
%%
count = 0;
v = 2;
x = x0;
Jmatrix = Jacobi( x );
A = Jmatrix' * Jmatrix;
g = Jmatrix' * fun( x );
found = ( norm( g , Inf ) < epsilon_1 );
mu = tau * max( diag( A ) );

while ~found && ( count < N )
    count = count + 1;
%     disp(count)
%     ( A + mu * I )可能是病态矩阵，解h_lm可能需要高精度算法
    h_lm = ( A + mu * diag( ones( 4 , 1 ) ) ) \ ( -g );
    if ( norm( h_lm ) <= ( epsilon_2 * ( norm( x ) + epsilon_2 ) ) )
        found = 1;
    else
        x_new = x + h_lm;
        gain_radio = ( fun( x )' * fun( x ) - fun( x_new )' * fun( x_new ) ) / 2 /...
            ( 0.5 * h_lm' * ( mu * h_lm - g ) );
        if gain_radio > 0
            x = x_new;
            Jmatrix = Jacobi( x );
            A = Jmatrix' * Jmatrix;
            g = Jmatrix' * fun( x );
            found = ( norm( g , Inf ) <= epsilon_1 );
            mu = mu * max( [ 1/3 , 1 - ( 2 * gain_radio - 1 )^3 ] );
            v = 2;
        else
            mu = mu * v;
            v = 2 * v;
        end
    end

% 如下注释部分用于迭代过程中显示拟合线向实验数据的接近过程，如需显示取消注释即可
%     tt = linspace(t(1),t(length(t)),100);
%     figure(2),plot(t,y,'*',tt,x(1) * exp( -( ( tt - x(2) ) / x(3) ).^2 / 2 ) + x(4),'r');
%     title(['iterative:',num2str(count)]);
%     pause(0.5);
end
x(2) = x(2) + num;
end
