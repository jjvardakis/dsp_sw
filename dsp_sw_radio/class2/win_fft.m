function [fout, faxis] = win_fft(iq, fs)

% def winfft(vec, fs=1., full_scale=1., beta=12., n_samp = None):
%     """
%     Returns Kaiser windowed fft with proper scaling of tones
% 
%     Dan Boschen 12/2/2018    
% 
%     Parameters
%     ----------
%     
%     vec : 1d list-like object
%         input time domain signal
%         
%     fs: float, optional
%         Sampling rate in Hz, default = 1.0
%         
%     full_scale: float, optional 
%         normalization for full scale, default = 1.0 (0 dB)
%         
%     beta: float, optional
%         Kaiser beta, default = 12.0 
%         
%     n_samp: integer, optional
%         Number of samples, if longer will zero pad vec
% 
%     Returns
%     -------
%     
%     faxis: ndarray, float
%         Frequency axis
%     
%     fout: ndarray, complex float
%         scaled frequency output values
% 
%     """
%     

full_scale = 1;
beta       = 12;
n_samp     = numel(iq);  

% if nargin <= 4
% 
% end
%    
% # length of input
input_len = numel(iq);
% 
%     if n_samp == None:
%         n_samp = input_len
%     elif n_samp < input_len:
%         vec = vec[:n_samp]
%         
% # only window data, not zero padded signal:
win_len = input_len;

% kaiser window
W = kaiser(win_len,beta).';
winscale = sum(W);

% windowed data
win_data = W .* iq;
%  

% # scaled fft
fout = 1./(winscale*full_scale) * fft(win_data, n_samp);
%     
%     # fftshift
%     fout = fft.fftshift(fout)
%     
% # create frequency axis
faxis = (0:n_samp-1) * fs / n_samp;
%     return faxis, fout

end