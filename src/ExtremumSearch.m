function [ locs , nFound ] = ExtremumSearch( x , y , n_extremum, peak_or_valley , threshold , extrenum_distance )
%% 搜索至多前n_extremum个极值，返回极值对应下标
% 所寻仅为大于阈值的极值
% -------------------------------------------------------------------------
% x                 - 数据横坐标
% y                 - 数据纵坐标
% threshold         - 阈值
% n_extremum        - 寻找的极值数目，序列的前n_extremum个极值
% extrenum_distance - 两极大/小值之间的最小距离，该距离内所有除最大/小值之外的其他极值均被认为是噪声引起的误差
% peak_or_valley    - peak（0） or valley（1）
% locs              - 返回极值的真实下标
% -------------------------------------------------------------------------
if extrenum_distance <=0
    error('extrenum_distance must be positive');
end

locs = NaN * zeros( n_extremum , 1 );
N = length( y );
%% 极值寻找
j = 0;
for i = 1 : N
    %% 只寻找大于阈值的极大值，或小于阈值的极小值
    if peak_or_valley == 0 % 寻找大于阈值的极大值
        if y(i) > threshold
            if ( i > 1 ) && ( i < N )
                if ( y(i) > y(i - 1) ) && ( y(i) > y(i + 1) )
                    j = j + 1;
                else
                    continue;
                end
            elseif i == 1
                if y(i) > y(i + 1)
                    j = j + 1;
                else
                    continue;
                end
            elseif i == N
                if y(i) > y(i - 1)
                    j = j + 1;
                else
                    continue;
                end
            end
        else
            continue;
        end
    elseif peak_or_valley == 1 % 寻找小于阈值的极小值
        if y(i) < threshold
            if ( i > 1 ) && ( i < N )
                if ( y(i) < y(i - 1) ) && ( y(i) < y(i + 1) )
                    j = j + 1;
                else
                    continue;
                end
            elseif i == 1
                if y(i) < y(i + 1)
                    j = j + 1;
                else
                    continue;
                end
            elseif i == N
                if y(i) < y(i - 1)
                    j = j + 1;
                else
                    continue;
                end
            end
        else
            continue;
        end
    end
    %% 依据设定的极值最小间距，排除噪声造成的毛刺，粗略获得极值位置
    if j == 1
        locs(j , 1) = i;
    else
        if ( x( i ) - x( locs(j - 1 , 1) ) ) < extrenum_distance
            if peak_or_valley == 0 && y( i ) > y( locs(j - 1 , 1) )
                j = j - 1;
                locs(j , 1) = i;
            elseif peak_or_valley == 1 && y( i ) < y( locs(j - 1 , 1) )
                j = j - 1;
                locs(j , 1) = i;
            else
                j = j - 1;
            end
        else
            if (j - 1) == n_extremum
                break;
            end
            locs(j , 1) = i;
        end
    end
end
%% 寻找到极值少于预定极值数时
if j < n_extremum
    locs = locs( 1 : j , 1 );
    nFound = j;
else
    nFound = n_extremum;
end

