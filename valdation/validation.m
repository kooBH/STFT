close all;

[in,fs] = audioread('../build/input.wav');
[out,fs] = audioread('../build/output.wav');

t_in = in(1:150000-384,1);
t_out = out(385:150000,1);

figure();plot(t_in-t_out);title(" input - output ");