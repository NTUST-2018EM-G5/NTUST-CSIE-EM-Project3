#include "FT.h"

FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			DFT(FreqReal, FreqImag, InputImage, M, N, j, i);
		}
	}
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v * y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i < M; ++i)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			InverseDFT(InverseReal, InverseImag, FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
			//存下反傅立葉實數與虛數部分
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; ++i)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v * y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FFT(int N, complex<double>* x)
{
	InverseFFT(N, x);
	for (int i = 0; i < N; ++i)
	{
		x[i] /= N;
	}
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M, N;
	M = h;
	N = w;
	complex<double>** result = new complex<double>*[M];

	for (int i = 0; i < N; ++i) result[i] = new complex<double>[N];

	for (int i = 0; i < M; ++i) for (int j = 0; j < N; ++j) result[i][j] = complex<double>(InputImage[j][i], 0);

	for (int i = 0; i < M; ++i)
	{
		complex<double> *input = new complex<double>[M];
		for (int j = 0; j < M; ++j) input[j] = result[j][i];
		FFT(M, input);
		for (int j = 0; j < M; ++j) result[j][i] = input[j];
	}
	for (int i = 0; i < N; ++i)
	{
		complex<double> *input = new complex<double>[N];
		for (int j = 0; j < N; ++j) input[j] = result[i][j];
		FFT(N, input);
		for (int j = 0; j < N; ++j) result[i][j] = input[j];
	}
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			FreqImag[i][j] = result[i][j].imag();
			FreqReal[i][j] = result[i][j].real();
			double temp1 = FreqImag[i][j] * FreqImag[i][j];
			double temp2 = FreqReal[i][j] * FreqReal[i][j];
			OutputImage[i][j] = sqrt(temp1 + temp2) * N;
		}
	}
}

void FT::InverseFFT(int N, complex<double>* x)
{
	int i = 1, j = 0;
	while (i < N)
	{
		int k = N >> 1;
		while (!((j ^= k)&k))
		{
			k = k >> 1;
		}
		if (i > j)
		{
			swap(x[i], x[j]);
		}
		++i;
	}
	for (int k = 2; k <= N; k <<= 1)
	{
		double omega = -2.0 * PI / (1.0*k);
		double cosOmega = cos(omega);
		double sinOmega = sin(omega);
		complex<double> dtheta(cosOmega, sinOmega);
		for (int j = 0; j < N; j += k)
		{
			complex<double> theta(1, 0);
			for (int i = j; i < j + (k / 2); ++i)
			{
				complex<double> a = x[i];
				complex<double> b = x[i + (k / 2)] * theta;
				x[i] = a + b;
				x[i + k / 2] = a - b;
				theta *= dtheta;
			}
		}
	}
}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M, N;
	M = h;
	N = w;
	double max_s = -100;
	complex<double>** result = new complex<double>*[M];

	for (int i = 0; i < N; ++i)result[i] = new complex<double>[N];

	for (int i = 0; i < M; ++i) for (int j = 0; j < N; ++j) result[i][j] = complex<double>(FreqReal[j][i], FreqImag[j][i]);

	for (int i = 0; i < M; ++i)
	{
		complex<double> *input = new complex<double>[M];
		for (int j = 0; j < M; ++j) input[j] = complex<double>(result[i][j].real()*N, result[i][j].imag()*N);
		InverseFFT(M, input);
		for (int j = 0; j < M; ++j) result[i][j] = input[j];
	}
	for (int i = 0; i < N; ++i)
	{
		complex<double> *input = new complex<double>[N];
		for (int j = 0; j < N; ++j) input[j] = complex<double>(result[j][i].real()*N, result[j][i].imag()*N);
		InverseFFT(N, input);
		for (int j = 0; j < N; ++j) result[j][i] = input[j];
	}
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			FreqImag[i][j] = result[i][j].imag();
			FreqReal[i][j] = result[i][j].real();
			double temp1 = result[i][j].real()*result[i][j].real();
			double temp2 = result[i][j].imag()*result[i][j].imag();
			OutputImage[i][j] = sqrt(temp1 + temp2);
			max_s < OutputImage[i][j] ? max_s = OutputImage[i][j] : 0;
		}
	}
	max_s = 1.0f / (max_s / 255);
	for (int i = 0; i < M; ++i) for (int j = 0; j < N; ++j) OutputImage[i][j] *= max_s;

	for (int i = M - 1; i >= M / 2; --i) for (int j = N - 1; j >= 0; --j) swap(OutputImage[i][j], OutputImage[M - i][j]);

	for (int i = M - 1; i >= 0; --i) for (int j = N - 1; j >= N / 2; --j) swap(OutputImage[i][j], OutputImage[i][N - j]);
}

void FT::HighpassFilter(int** InputImage, int** OutputImage, double** Real, double** Img, int h, int w)
{
	int threshold, x, y;
	threshold = h * 0.125;
	x = (h / 2);
	y = (w / 2);
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			int temp1 = (i - x)*(i - x);
			int temp2 = (j - y)*(j - y);
			int d = sqrt(temp1 + temp2);
			bool highpass = d > threshold;

			OutputImage[i][j] = highpass ? InputImage[i][j] : 0;
			if (!highpass) //NOT HIGH PASS
			{
				Real[i][j] = 0;
				Img[i][j] = 0;
			}
		}
	}
}

void FT::LowpassFilter(int** InputImage, int** OutputImage, double** Real, double** Img, int h, int w)
{
	int threshold, x, y;
	threshold = h * 0.125;
	x = (h / 2);
	y = (w / 2);
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			int temp1 = (i - x)*(i - x);
			int temp2 = (j - y)*(j - y);
			int d = sqrt(temp1 + temp2);
			bool lowpass = d < threshold;

			OutputImage[i][j] = lowpass ? InputImage[i][j] : 0;
			if (!lowpass) //NOT LOW PASS
			{
				Real[i][j] = 0;
				Img[i][j] = 0;
			}
		}
	}
}