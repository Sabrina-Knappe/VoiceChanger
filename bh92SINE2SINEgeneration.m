function [bh92SINE2SINE,bh92SINE2SINEsize] = bh92SINE2SINEgeneration;
% function [bh92SINE2SINE,bh92SINE2SINEsize] =bh92SINE2SINEgeneration;
%
% ==> generation of the Blackman-Harris window
% output data:
% bh92SINE2SINEsize: size of the window
% bh92SINE2SINE: (sampled) window
bh92SINE2SINEsize = 4096;
bh92SINE2SINE = zeros(bh92SINE2SINEsize,1);
bh92N = 512;
bh92const = [.35875, .48829, .14128, .01168];
bh92Theta = -4*2*pi/bh92N;
bh92ThetaIncr = 8*2*pi/bh92N/bh92SINE2SINEsize;
for i=1:bh92SINE2SINEsize
 for m=0:3
 bh92SINE2SINE(i)=bh92SINE2SINE(i)-bh92const(m+1)/2*(sine2sine(bh92Theta-m*2*pi/bh92N,bh92N)+sine2sine(bh92Theta+m*2*pi/bh92N,bh92N));
 end;
 bh92Theta = bh92Theta + bh92ThetaIncr;
end;
bh92SINE2SINE = bh92SINE2SINE/bh92SINE2SINE(bh92SINE2SINEsize/2+1);