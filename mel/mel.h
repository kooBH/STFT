

// librosa-Slaney style mel-filter bank
// https://librosa.org/doc/main/_modules/librosa/filters.html#mel  
//

/*
% librosa - Slaney-style
% https://librosa.org/doc/main/_modules/librosa/filters.html#mel
% confirmed exact output with librosa default mel-filterbank

%{
Note that there are different ways to compute MFCC
@article{article,
author = {Ganchev, Todor and Fakotakis, Nikos and George, Kokkinakis},
year = {2005},
month = {01},
pages = {},
title = {Comparative evaluation of various MFCC implementations on the speaker verification task},
volume = {1},
journal = {Proceedings of the SPECOM}
}
%}


close all;
clear

fs = 16000;
nfft = 512;

nhfft = 257;
n_mels = 40;
F = linspace(0,fs/2,nhfft);

%% fft_frequencies()
fftfreqs  = linspace(0, fs / 2, int32(1 + nfft / 2));

%% mel_frequencies()
min_mel  = hz_to_mel(0);
max_mel  = hz_to_mel(8000);
mels  = linspace(min_mel,max_mel,n_mels+2);
mel_f  = mel_to_hz(mels);

fdiff = diff(mel_f);

%  ramps = np.subtract.outer(mel_f, fftfreqs)
ramps  = mel_f' - fftfreqs;

weights = zeros(n_mels,nhfft);

for i = 1:n_mels
    % lower and upper slopes for all bins
    lower = -ramps(i,:) / fdiff(i);
    upper = ramps(i + 2,:) / fdiff(i + 1);

    % .. then intersect them with each other and zero
    weights(i,:) = max(0, min(lower, upper));
end

% Slaney-style mel is scaled to be approx constant energy per channel
enorm = 2.0 ./ (mel_f(3 : n_mels + 2) - mel_f(1:n_mels));
weights = times(weights,enorm');

disp(max(max(weights)));

figure()
hold on;
%% Construct filter bank
for i = 1:n_mels
    plot(weights(i,:));
end
hold off;
title("mel weights");


function w = triangle(N)
    if rem(N,2) % odd
        w_temp = 2*(1:(N+1)/2)/(N+1);
        w = [w_temp w_temp((N-1)/2:-1:1)]';
    else % even
        w_temp = (2*(1:(N+1)/2)-1)/N;
        w = [w_temp w_temp(N/2:-1:1)]';
    end

end

function mels = hz_to_mel(f)

% Fill in the linear part
f_min = 0.0;
f_sp = 200.0 / 3;

mels = (f - f_min) / f_sp;

% Fill in the log-scale part
min_log_hz = 1000.0;  % beginning of log region (Hz)
min_log_mel = (min_log_hz - f_min) / f_sp;  % same (Mels)
logstep = log(6.4) / 27.0;  % step size for log region

if f >= min_log_hz
    mels = min_log_mel + log(f / min_log_hz) / logstep;
end

end

function freqs = mel_to_hz(mels)

    % Fill in the linear scale
    f_min = 0.0;
    f_sp = 200.0 / 3;
    freqs = f_min + f_sp * mels;

    % And now the nonlinear scale
    min_log_hz = 1000.0;  % beginning of log region (Hz)
    min_log_mel = (min_log_hz - f_min) / f_sp;  % same (Mels)
    logstep = log(6.4) / 27.0;  % step size for log region

    if ~isscalar(mels)
        % If we have vector data, vectorize
        log_t = mels >= min_log_mel;
        freqs(log_t) = min_log_hz * exp(logstep * (mels(log_t) - min_log_mel));
    elseif mels >= min_log_mel
        % If we have scalar data, check directly
        freqs = min_log_hz * np.exp(logstep * (mels - min_log_mel));
    end
end
*/

class mel{

private:
  int n_mels;
  int sr;
  int fft;
  int nhfft;

  double** filter_bank;
public : 
  inline mel();
  inline ~mel();

  inline int filter(double**stft,double**out);

  inline static double** filterBank(int sr,int n_mels, int nfft);
};

mel::mel() {

}

mel::~mel() {

}

int mel::filter(double** stft, double** out) {

}

double** mel::filterBank(int sr_, int n_mels_, int nfft_) {

}