% Based on  https://kr.mathworks.com/matlabcentral/fileexchange/23179-melfilter
% HTK style ? 
%% Author Information
%   Pierce Brady
%   Smart Systems Integration Group - SSIG
%	Cork Institute of Technology, Ireland.
%  



close all;
clear

fs = 16000;
nfft = 512;

nhfft = 257;
N = 40;
F = linspace(0,fs/2,nhfft);

MelFrequencyVector = 2595*log10(1+F/700);

figure;
plot(MelFrequencyVector);title("mel frequency vector");

MaxF = max(MelFrequencyVector);           
MinF = min(MelFrequencyVector);          
MelBinWidth = (MaxF-MinF)/(N+1);

Filter = zeros([N nhfft]); % filter bank alloc
size(Filter)

figure()
hold on;
%% Construct filter bank
for i = 1:N
    iFilter = find(MelFrequencyVector>=((i-1)*MelBinWidth+MinF) & ...
                    MelFrequencyVector<=((i+1)*MelBinWidth+MinF));
    disp([num2str(numel(iFilter)) ' : ' num2str(min(iFilter)) ' ~ ' num2str(max(iFilter))])
    w = triangle(numel(iFilter));
    plot(min(iFilter):max(iFilter),w);
    Filter(i,iFilter) = w(numel(iFilter)); % Triangle window
end
hold off;
title("Triangluar Window");


function w = triangle(N)
    if rem(N,2) % odd
        w_temp = 2*(1:(N+1)/2)/(N+1);
        w = [w_temp w_temp((N-1)/2:-1:1)]';
    else % even
        w_temp = (2*(1:(N+1)/2)-1)/N;
        w = [w_temp w_temp(N/2:-1:1)]';
    end

end