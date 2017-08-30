function [ Y ] = soft_thresh( X, thresh )
%SOFT_THRESH Summary of this function goes here
%   Detailed explanation goes here

Y = max(X-thresh,0);
Y = Y + min(X+thresh,0);


end

