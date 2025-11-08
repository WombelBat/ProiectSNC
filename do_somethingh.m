function [y] = do_somethingh(y)
y = 20*log10(abs(fft(y))/length(y) )
end