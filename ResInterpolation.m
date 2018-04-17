
% this function is for interpolate nir infrared iamges,
function img_interpolated=ResInterpolation(green, channel, mask, h, v, eps)

% channel residual interpolation
channel_hat = guidedfilter(green, channel, mask,h,v,eps);
res_channel = (channel-channel_hat).*mask;
H = [1/4,1/2,1/4;1/4,1/2,1/4;1/4,1/2,1/4];
res_channel = imfilter(res_channel,H,'replicate');
img_interpolated = res_channel + channel_hat;

end