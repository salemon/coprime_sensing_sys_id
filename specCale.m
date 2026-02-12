function varargout=specCale(idata,Fs,NFFT)
% function dataOut=Xu_specCal(idata,Fs)
% this file calculates the spectrum of the input data
%   Fs: sampling frequency
% 
%   Copyright (c) 2008-, 
%   Xu Chen

b_autoNFFT = 0;
if nargin < 3
    b_autoNFFT = 1;
    if nargin < 2
        Fs = 26400;
    end
end
if size(idata,4)>1
    for ii = 1:size(idata,4)
        L       = length(idata(:,:,:,ii));
        if b_autoNFFT
            NFFT    = 2^nextpow2(L);
        end
        Y       = fft(idata(:,:,:,ii),NFFT)/L;
        f       = Fs/2*linspace(0,1,NFFT/2);
        
        amp = 2*abs(Y(1:NFFT/2));
        if 0
            %%
            figure, plot(amp)
        end
        if rem(NFFT, 2) % odd nfft excludes Nyquist point
            amp(1)      = amp(1)/2;
        else
            amp(end)    = amp(end)/2;
        end
        dataOut.f(:,ii)   = f'; %Hz
        % dataOut.y = 20*log10(amp); %dB
        dataOut.amp(:,ii) = amp;
        dataOut.pha(:,ii) = angle(Y(1:NFFT/2))/pi*180;
        ilegend{ii} = ['3\sigma = ',num2str(3*std(idata(:,:,:,ii)))];
    end
else
    if size(idata,2) > size(idata,1)
        idata = idata';
    end
    for ii = 1:size(idata,2)
        L       = length(idata(:,ii));
        if b_autoNFFT
            NFFT    = 2^nextpow2(L);
        end
        Y       = fft(idata(:,ii),NFFT)/L;
        f       = Fs/2*linspace(0,1,NFFT/2);
        
        amp = 2*abs(Y(1:NFFT/2));
        if 0
            %%
            figure, plot(amp)
        end
        if rem(NFFT, 2) % odd nfft excludes Nyquist point
            amp(1)      = amp(1)/2;
        else
            amp(end)    = amp(end)/2;
        end
        dataOut.f(:,ii)   = f'; %Hz
        % dataOut.y = 20*log10(amp); %dB
        dataOut.amp(:,ii) = amp;
        dataOut.pha(:,ii) = angle(Y(1:NFFT/2))/pi*180;
        ilegend{ii} = ['3\sigma = ',num2str(3*std(idata(:,ii)))];
    end
end
if nargout == 1
    varargout{1} = dataOut;
elseif nargout == 0
    %     figure,
    %     plot(f,amp)
    
    figure,
    plot(dataOut.f,dataOut.amp)
    xlim([0,max(dataOut.f(:))])
    legend(ilegend)
    %
    %     figure,
    %     semilogx(f',mag2db(amp))
else
    error('Output error');
end
%2011-07-18
% L       = length(idata);
% NFFT    = 2^nextpow2(L);
% Y       = fft(idata,NFFT)/L;
% f       = Fs/2*linspace(0,1,NFFT/2);
%
% amp = 2*abs(Y(1:NFFT/2));
% if rem(NFFT, 2) % odd nfft excludes Nyquist point
%     amp(1)      = amp(1)/2;
% else
%     amp(end)    = amp(end)/2;
% end
% dataOut.f   = f; %Hz
% % dataOut.y = 20*log10(amp); %dB
% dataOut.amp = amp;
