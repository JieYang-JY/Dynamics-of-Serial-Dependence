function [sm_output] = getGaussResults(raw_input, sigma)

half_gaus = sigma * 3;
dist = (-(half_gaus):( half_gaus))';
gauss = 1 / (sigma*sqrt(2*pi))*exp(-((dist.^2)/(2*sigma^2)));

input_buffered = [raw_input(half_gaus:-1:1); raw_input; raw_input(end:-1:end-half_gaus+1)];
sm_output = convn(input_buffered, gauss, 'valid');


