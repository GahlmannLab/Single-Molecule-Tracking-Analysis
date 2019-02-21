function y = nanMedian(x)
%NANMEDIAN Median value, ignoring NaNs.
%   y = nanMedian(x) returns the sample median of X, treating NaNs as
%   missing values.  


y = median(x(~isnan(x)));

end


