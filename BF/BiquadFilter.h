#ifndef _H_BIQUAD_H_
#define _H_BIQUAD_H_

#include <cmath>
#include <stdio.h>

/*
%% https://www.earlevel.com/main/2021/09/02/biquad-calculator-v3/

		V = 10^(abs(peakGain)/20);
		K = tan(pi*Fc/Fs);
		sqrt2 = sqrt(2);

	if type == 1

		%% case "one-pole lp":
		a1 = exp(-2.0 * pi * (Fc / Fs));
				b0 = 1.0 - a1;
				a1 = -a1;
		b1 = 0;
				b2 = 0;
				a2 = 0;

		elseif type == 2

		%% case "one-pole hp":
		a1 = -exp(-2.0 * pi * (0.5 - Fc / Fs));
				b0 = 1.0 + a1;
				a1 = -a1;
		b1 = 0;
				b2 = 0;
				a2 = 0;

		elseif type == 3

		%% case "lowpass 1p1z":
		norm = 1 / (1 / K + 1);
		b0 = norm;
				b1 = norm;
		a1 = (1 - 1 / K) * norm;
		b2 = 0;
				a2 = 0;

	elseif type == 4

		%% case "highpass 1p1z":
		norm = 1 / (K + 1);
		b0 = norm;
		b1 = -norm;
		a1 = (K - 1) * norm;
		b2 = 0;
				a2 = 0;

	elseif type == 5

		%% case "lowpass":
		norm = 1 / (1 + K / Q + K * K);
		b0 = K * K * norm;
		b1 = 2 * b0;
		b2 = b0;
		a1 = 2 * (K * K - 1) * norm;
		a2 = (1 - K / Q + K * K) * norm;

		elseif type == 6

		%% case "highpass":
		norm = 1 / (1 + K / Q + K * K);
		b0 = 1 * norm;
		b1 = -2 * b0;
		b2 = b0;
		a1 = 2 * (K * K - 1) * norm;
		a2 = (1 - K / Q + K * K) * norm;

	elseif type == 7

		%% case "bandpass":
		norm = 1 / (1 + K / Q + K * K);
		b0 = K / Q * norm;
		b1 = 0;
		b2 = -b0;
		a1 = 2 * (K * K - 1) * norm;
		a2 = (1 - K / Q + K * K) * norm;

		elseif type == 8

		%% case "notch":
		norm = 1 / (1 + K / Q + K * K);
		b0 = (1 + K * K) * norm;
		b1 = 2 * (K * K - 1) * norm;
		b2 = b0;
		a1 = b1;
		a2 = (1 - K / Q + K * K) * norm;

		elseif type == 9

		%% case "peak":
		if peakGain >= 0
			norm = 1 / (1 + 1/Q * K + K * K);
			b0 = (1 + V/Q * K + K * K) * norm;
			b1 = 2 * (K * K - 1) * norm;
			b2 = (1 - V/Q * K + K * K) * norm;
			a1 = b1;
			a2 = (1 - 1/Q * K + K * K) * norm;
		else
			norm = 1 / (1 + V/Q * K + K * K);
			b0 = (1 + 1/Q * K + K * K) * norm;
			b1 = 2 * (K * K - 1) * norm;
			b2 = (1 - 1/Q * K + K * K) * norm;
			a1 = b1;
			a2 = (1 - V/Q * K + K * K) * norm;
				end

	elseif type == 10

		%% case "lowShelf":
		if peakGain >= 0
			norm = 1 / (1 + sqrt2 * K + K * K);
			b0 = (1 + sqrt(2*V) * K + V * K * K) * norm;
			b1 = 2 * (V * K * K - 1) * norm;
			b2 = (1 - sqrt(2*V) * K + V * K * K) * norm;
			a1 = 2 * (K * K - 1) * norm;
			a2 = (1 - sqrt2 * K + K * K) * norm;

		else
			norm = 1 / (1 + sqrt(2*V) * K + V * K * K);
			b0 = (1 + sqrt2 * K + K * K) * norm;
			b1 = 2 * (K * K - 1) * norm;
			b2 = (1 - sqrt2 * K + K * K) * norm;
			a1 = 2 * (V * K * K - 1) * norm;
			a2 = (1 - sqrt(2*V) * K + V * K * K) * norm;
				end

	elseif type == 11

		%% case "highShelf":
				if peakGain >= 0
						norm = 1 / (1 + sqrt2 * K + K * K);
						b0 = (V + sqrt(2*V) * K + K * K) * norm;
						b1 = 2 * (K * K - V) * norm;
						b2 = (V - sqrt(2*V) * K + K * K) * norm;
						a1 = 2 * (K * K - 1) * norm;
						a2 = (1 - sqrt2 * K + K * K) * norm;

				else
						norm = 1 / (V + Math.sqrt(2*V) * K + K * K);
						b0 = (1 + sqrt2 * K + K * K) * norm;
						b1 = 2 * (K * K - 1) * norm;
						b2 = (1 - sqrt2 * K + K * K) * norm;
						a1 = 2 * (K * K - V) * norm;
						a2 = (V - sqrt(2*V) * K + K * K) * norm;
				end

	elseif type == 12

		%% case "lowShelf 1st":
		if peakGain >= 0
			norm = 1 / (K + 1);
			b0 = (K * V + 1) * norm;
			b1 = (K * V - 1) * norm;
			b2 = 0;
			a1 = (K - 1) * norm;
			a2 = 0;
		else
			norm = 1 / (K * V + 1);
			b0 = (K + 1) * norm;
			b1 = (K - 1) * norm;
			b2 = 0;
			a1 = (K * V - 1) * norm;
			a2 = 0;
				end

	elseif type == 13

		%% case "highShelf 1st":
				if peakGain >= 0
			norm = 1 / (K + 1);
			b0 = (K + V) * norm;
			b1 = (K - V) * norm;
			b2 = 0;
			a1 = (K - 1) * norm;
			a2 = 0;
				else
			norm = 1 / (K + V);
			b0 = (K + 1) * norm;
			b1 = (K - 1) * norm;
			b2 = 0;
			a1 = (K - V) * norm;
			a2 = 0;
				end

	elseif type == 14

		%% case "allpass":
		norm = 1 / (1 + K / Q + K * K);
		b0 = (1 - K / Q + K * K) * norm;
		b1 = 2 * (K * K - 1) * norm;
		b2 = 1;
		a1 = b1;
		a2 = b0;

		elseif type == 15

		%% case "allpass 1st":
		b0 = (1 - K) / (1 + K);
		b1 = -1;
		b2 = 0;
		a1 = -b0;
		a2 = 0;

		end

		a = [1.0; a1; a2];  % pole
		b = [b0; b1; b2];   % zero


*/
class BiquadFilter {

private :
	const double pi = 3.141591;

	double a[3] = { 1.0, 0.0, 0.0 }; // pole
	double b[3] = { 0.0, 0.0, 0.0 }; // zero
	const int p = 2;
	const int q = 2;
	const int pq = 2;
	double u[2] = { 0.0, 0.0 }; // internal state

	double* tmp = nullptr;
	size_t size_tmp=0;

public :
	BiquadFilter(int type, double Fc, double Fs, double Q, double peakGain);
	~BiquadFilter();
	enum TYPE{
		one_pole_lp,
		one_pole_hp,
		lowpass_1p1z,
		highpass_1p1z,
		lowpass,
		highpass,
		bandpass,
		notch,
		peak,
		lowShelf,
		highShelf,
		lowShelf_1st,
		highSgelf_1st,
		allpass,
		allpass_1st
	};

	void Filter(double* x, int n_sample);
	void Reset();
};

BiquadFilter::BiquadFilter(int type, double Fc, double Fs, double Q, double peakGain){
	double V, K, sqrt2;
	double norm, b0, b1, b2, a1, a2;

	V = std::pow(10,std::abs(peakGain)/20);
	K = std::tan(pi * Fc / Fs);
	sqrt2 = std::sqrt(2);

	switch (type) {
	case BiquadFilter::TYPE::highpass:
		norm = 1.0 / (1.0 + K/Q + K*K);
		b0 = 1 * norm;
		b1 = -2 * b0;
		b2 = b0;
		a1 = 2.0 * (K*K - 1) * norm;
		a2 = (1 - K/Q + K*K) * norm;

		break;
	case BiquadFilter::TYPE::peak:
		if (peakGain >= 0.0) {
			norm = 1.0 / (1.0 + 1.0/Q*K + K*K);
			b0 = (1 + V / Q * K + K * K) * norm;
			b1 = 2.0 * (K * K - 1) * norm;
			b2 = (1 - V / Q * K + K * K) * norm;
			a1 = b1;
			a2 = (1 - 1 / Q * K + K * K) * norm;
		}else{
			norm = 1.0 / (1.0 + V/Q*K + K*K);
			b0 = (1 + 1 / Q * K + K * K) * norm;
			b1 = 2.0 * (K * K - 1) * norm;
			b2 = (1 - 1 / Q * K + K * K) * norm;
			a1 = b1;
			a2 = (1 - V / Q * K + K * K) * norm;
		}
		break;
	default : 
		a1 = 0;
		a2 = 0;
		b0 = 0;
		b1 = 0;
		b2 = 0;
		printf("BiquadFilter type %d is not implemented\n",type);
		break;
	}

	a[1] = a1;
	a[2] = a2;
	b[0] = b0;
	b[1] = b1;
	b[2] = b2;
	Reset();
}

BiquadFilter::~BiquadFilter() {
		if (tmp)delete[] tmp;
}

/*
Implement Type-2 only

% Ref. : Boaz Porat, "A Course In Digital Signal Processing", p.440
%
% Synopsis: direct(type,b,a,x)
% Direct realization of rational transfer functions.
%
% Input parameters:
% typ: 1 for direct realization, 2 for transposed
% b, a: numerator and denominator polynomials
% x: input sequence.
%
% Output:
% y: output sequence.

p = length(a)-1;
q = length(b)-1;
pq = max(p,q);
a = a(2:p+1);
u = zeros(1,pq); % u: the internal state
y = zeros(length(x),1);

if typ == 1

		for i = 1 : length(x)
				unew = x(i) - sum(u(1:p).*a);
				u = [unew, u];
				y(i) = sum(u(1:q+1).*b);
				u = u(1:pq);
		end

elseif typ == 2

		for i = 1 : length(x)
				y(i) = b(1)*x(i)+u(1);
				u = [u(2:pq), 0];
				u(1:q) = u(1:q) + b(2:q+1) * x(i);
				u(1:p) = u(1:p) - a*y(i);
		end
*/



void BiquadFilter::Filter(double* x, int n_sample) {
	if (size_tmp < n_sample) {
		if (tmp)delete[] tmp;
		tmp = new double[n_sample];
		size_tmp = n_sample;
	}
	memcpy(tmp, x, sizeof(double) * n_sample);


	// type-2 
	for (int i = 0; i < n_sample; i++) {
		// y(i) = b(1)*x(i)+u(1);
		x[i] = b[0] * tmp[i] + u[0];

		// u = [u(2:pq), 0];
		u[0] = u[1];
		u[1] = 0;

		// u(1:q) = u(1:q) + b(2:q+1) * x(i);
		u[0] = u[0] + b[1] * tmp[i];
		u[1] = u[1] + b[2] * tmp[i];

		// u(1:p) = u(1:p) - a*y(i);
		u[0] = u[0] - a[1] * x[i];
		u[1] = u[1] - a[2] * x[i];
	}
}

void BiquadFilter::Reset() {
	u[0] = 0;
	u[1] = 0;
}

#endif