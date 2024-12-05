function [smoothed_data] = smoothran(x,y,seconds)
% SMOOTHRAN smoothes data on a moving average
%
% [smoothed_data] = smoothran(x,y,seconds)
%   x = time_min, a time series in minutes
%   y = data
%   seconds = the seconds around the point averaged
ran = seconds/60; % Minutes
ysm = zeros(length(x),1);
for i = 1:length(x)
    if (x(i)-x(1)) <= ran/2
        ysm(i)=nan;
    elseif x(i) >= (x(end)-ran/2)
        ysm(i)=nan;
    else
        z = find(x>(x(i)-ran/2) & x<(x(i)+ran/2)); % Indicies of x +/- n
        ysm(i) = mean(y(z));
    end
end
smoothed_data = ysm;
end
