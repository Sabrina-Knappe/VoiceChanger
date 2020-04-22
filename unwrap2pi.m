function argunwrap = unwrap2pi (arg)
% function argunwrap = unwrap2pi (arg)
%
%==> unwrapping of t h e p h a s e , i n [ - p i , p i ]
% a r g : p h a s e t o unwrap
arg = arg - floor(arg/2/pi)*2*pi;
argunwrap = arg - (arg>=pi)*2*pi;