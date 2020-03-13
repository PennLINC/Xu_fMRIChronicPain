function sdata = my_MemSmooth64bit(data,fwhm,V,sdata)
% FORMAT sdata = MemSmooth64bit(data,fwhm,Vo)
% Vo    - mapped image describing the volume loaded in data
% fwhm  - FWHM of Guassian filter width in mm
% data  - unsmoothed 3D matrix
%
%_______________________________________________________________________
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_smoothto8bit.m 112 2005-05-04 18:20:52Z john $

vx   = sqrt(sum(V.mat(1:3,1:3).^2));
s    = (fwhm./vx./sqrt(8*log(2)) + eps).^2;
r    = cell(1,3);
for i=1:3,
	r{i}.s = ceil(3.5*sqrt(s(i)));
	x      = -r{i}.s:r{i}.s;
	r{i}.k = exp(-0.5 * (x.*x)/s(i))/sqrt(2*pi*s(i));
	r{i}.k = r{i}.k/sum(r{i}.k);
end;

buff = zeros([V.dim(1:2) r{3}.s*2+1]);

for i=1:V.dim(3)+r{3}.s,
	if i<=V.dim(3),
		buff(:,:,rem(i-1,r{3}.s*2+1)+1) = ...
			conv2(conv2(data(:,:,i),r{1}.k,'same'),r{2}.k','same');
	else,
		buff(:,:,rem(i-1,r{3}.s*2+1)+1) = 0;
	end;

	if i>r{3}.s,
		kern    = zeros(size(r{3}.k'));
		kern(rem((i:(i+r{3}.s*2))',r{3}.s*2+1)+1) = r{3}.k';
		sdata(:,:,i-r{3}.s)     = reshape(reshape(buff,[prod(V.dim(1:2)) r{3}.s*2+1])*kern,V.dim(1:2));
	end;
end;
