#pragma once
#include <iostream>
#include <complex>
#define PI 3.14159265359
using namespace std;
class FT
{
	public:
		FT();
		void DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
		void DFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

		void InverseDiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
		void InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

		void FastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
		void FFT(int N, std::complex<double>* x);

		void InverseFastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
		void InverseFFT(int N, std::complex<double>* x);

		void LowpassFilter(int** InputImage, int** OutputImage, double** Real, double** Img, int h, int w);
		void HighpassFilter(int** InputImage, int** OutputImage, double** Real, double** Img, int h, int w);

	private:

};



