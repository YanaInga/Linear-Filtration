#include <iostream>
#include <omp.h>
#include <string>
#include <Windows.h>
#include <math.h>
#include "BMPFileRW.h"
using namespace std;

double ArmV(double* time) //Среднее арифметическое значение
{
	double amvTime = 0;
	for (int i = 0; i < 20; i++)
	{
		amvTime += time[i];
	}
	return (double)amvTime / 20;
}
// Функция рассчета среднеарифметического значения в доверительном интервале
double AvgTrustedInterval(double& amvTime, double* time)
{
	double sd = 0, newAVg = 0;
	int newiter = 0;
	for (int i = 0; i < 20; i++)
	{
		sd += (time[i] - amvTime) * (time[i] - amvTime);
	}
	sd /= 19.0;
	sd = sqrt(sd);
	for (int i = 0; i < 20; i++)
	{
		if (amvTime - sd <= time[i] && time[i] <= amvTime + sd)
		{
			newAVg += time[i];
			newiter++;
		}
	}
	if (newiter == 0) newiter = 1;
	return newAVg / newiter;
}
void print(int count_of_threads, int size, double amvTime, string func)
{
	cout << "Ksize: " << size << func << ". Потоков - " << count_of_threads << " Время: " << amvTime << endl;
}
void findcoeff(double**& MATR_coef, int RH, int RW, int ksize)
{
	double sum = 0;
	for (int y = -RH; y <= RH; y++)
		for (int x = -RW; x <= RW; x++) {
			MATR_coef[y + RH][x + RW] = (1 / (2 * 3.14 * ksize * ksize)) * exp(-1 * (x * x + y * y) / (ksize * ksize * 2));
			sum += MATR_coef[y + RH][x + RW];
		}
	for (int y = -RH; y <= RH; y++)
		for (int x = -RW; x <= RW; x++) {
			MATR_coef[y + RH][x + RW] /= sum;
		}
}
void findcoeffPar(double**& MATR_coef, int RH, int RW, int ksize)
{
	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for (int y = -RH; y <= RH; y++)
		for (int x = -RW; x <= RW; x++) {
			MATR_coef[y + RH][x + RW] = (1 / (2 * 3.14 * ksize * ksize)) * exp(-1 * (x * x + y * y) / (ksize * ksize * 2));
			sum += MATR_coef[y + RH][x + RW];
		}
#pragma omp parallel for
	for (int y = -RH; y <= RH; y++)
		for (int x = -RW; x <= RW; x++) {
			MATR_coef[y + RH][x + RW] /= sum;
		}
}
double linefilterposl(int ksize, int height, int width, RGBQUAD** rgb1, RGBQUAD** rgb2)
{
	double time_start = omp_get_wtime();
	int size = ksize * ksize;
	int rh = (ksize - 1) / 2, rw = (ksize - 1) / 2;
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			int LinF_Valuer = 0;
			int LinF_Valueg = 0;
			int LinF_Valueb = 0;
			for (int dy = -rh; dy <= rh; dy++)
			{
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++)
				{
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;
					LinF_Valuer += rgb1[ky][kx].rgbRed;
					LinF_Valueg += rgb1[ky][kx].rgbGreen;
					LinF_Valueb += rgb1[ky][kx].rgbBlue;
				}
			}
			LinF_Valuer = LinF_Valuer / ((rh * 2 + 1) * (rw * 2 + 1));
			LinF_Valueg = LinF_Valueg / ((rh * 2 + 1) * (rw * 2 + 1));
			LinF_Valueb = LinF_Valueb / ((rh * 2 + 1) * (rw * 2 + 1));
			if (LinF_Valuer < 0) LinF_Valuer = 0;
			if(LinF_Valuer > 255) LinF_Valuer = 255;
			if (LinF_Valueg < 0) LinF_Valueg = 0;
			if (LinF_Valueg > 255) LinF_Valueg = 255;
			if (LinF_Valueb < 0) LinF_Valueb = 0;
			if (LinF_Valueb > 255) LinF_Valueb = 255;
			rgb2[y][x].rgbRed = LinF_Valuer;
			rgb2[y][x].rgbGreen = LinF_Valueg;
			rgb2[y][x].rgbBlue = LinF_Valueb;
		}
	}	
	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}
double linefilter(int ksize, int height, int width, RGBQUAD** rgb1, RGBQUAD**& rgb2)
{
	double time_start = omp_get_wtime();
	int size = ksize * ksize;
	int rh = (ksize - 1) / 2, rw = (ksize - 1) / 2;
#pragma omp parallel for schedule(dynamic, height/12)
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			int LinF_Valuer = 0;
			int LinF_Valueg = 0;
			int LinF_Valueb = 0;
			for (int dy = -rh; dy <= rh; dy++)
			{
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++)
				{
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;
					LinF_Valuer += rgb1[ky][kx].rgbRed;
					LinF_Valueg += rgb1[ky][kx].rgbGreen;
					LinF_Valueb += rgb1[ky][kx].rgbBlue;
				}
			}
			LinF_Valuer = LinF_Valuer / ((rh * 2 + 1) * (rw * 2 + 1));
			LinF_Valueg = LinF_Valueg / ((rh * 2 + 1) * (rw * 2 + 1));
			LinF_Valueb = LinF_Valueb / ((rh * 2 + 1) * (rw * 2 + 1));
			if (LinF_Valuer < 0) LinF_Valuer = 0;
			if (LinF_Valuer > 255) LinF_Valuer = 255;
			if (LinF_Valueg < 0) LinF_Valueg = 0;
			if (LinF_Valueg > 255) LinF_Valueg = 255;
			if (LinF_Valueb < 0) LinF_Valueb = 0;
			if (LinF_Valueb > 255) LinF_Valueb = 255;
			rgb2[y][x].rgbRed = LinF_Valuer;
			rgb2[y][x].rgbGreen = LinF_Valueg;
			rgb2[y][x].rgbBlue = LinF_Valueb;
		}
	}
	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}
double gaussfilter(int ksize, int height, int width, RGBQUAD** rgb1, RGBQUAD**& rgb2)
{
	double time_start = omp_get_wtime();
	int size = ksize * ksize;
	int rh = (ksize - 1) / 2, rw = (ksize - 1) / 2;
	double** Matr_coef = new double* [ksize];
#pragma omp parallel for
	for (int i = 0; i < ksize; i++)
		Matr_coef[i] = new double[ksize];
	findcoeffPar(Matr_coef, rh, rw, ksize);
#pragma omp parallel for schedule(dynamic, height/12)
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			double LinF_Valuer = 0.0;
			double LinF_Valueg = 0.0;
			double LinF_Valueb = 0.0;
			for (int dy = -rh; dy <= rh; dy++)
			{
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++)
				{
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;
					LinF_Valuer += rgb1[ky][kx].rgbRed * Matr_coef[dy + rh][dx + rw];
					LinF_Valueg += rgb1[ky][kx].rgbGreen * Matr_coef[dy + rh][dx + rw];
					LinF_Valueb += rgb1[ky][kx].rgbBlue * Matr_coef[dy + rh][dx + rw];
				}
			}
			if (LinF_Valuer < 0) LinF_Valuer = 0;
			if (LinF_Valuer > 255) LinF_Valuer = 255;
			if (LinF_Valueg < 0) LinF_Valueg = 0;
			if (LinF_Valueg > 255) LinF_Valueg = 255;
			if (LinF_Valueb < 0) LinF_Valueb = 0;
			if (LinF_Valueb > 255) LinF_Valueb = 255;
			rgb2[y][x].rgbRed = (int)LinF_Valuer;
			rgb2[y][x].rgbGreen = (int)LinF_Valueg;
			rgb2[y][x].rgbBlue = (int)LinF_Valueb;
		}
	}
	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}

double gaussposl(int ksize, int height, int width, RGBQUAD** rgb1, RGBQUAD** rgb2)
{
	double time_start = omp_get_wtime();
	int size = ksize * ksize;
	int rh = (ksize - 1) / 2, rw = (ksize - 1) / 2;
	double** Matr_coef = new double*[ksize];
	for (int i = 0; i < ksize; i++)
		Matr_coef[i] = new double[ksize];
	findcoeff(Matr_coef, rh, rw, ksize);
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			double LinF_Valuer = 0.0;
			double LinF_Valueg = 0.0;
			double LinF_Valueb = 0.0;
			for (int dy = -rh; dy <= rh; dy++)
			{
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++)
				{
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;
					LinF_Valuer += rgb1[ky][kx].rgbRed * Matr_coef[dy+rh][dx+rw];
					LinF_Valueg += rgb1[ky][kx].rgbGreen * Matr_coef[dy + rh][dx + rw];
					LinF_Valueb += rgb1[ky][kx].rgbBlue * Matr_coef[dy + rh][dx + rw];
				}
			}
			if (LinF_Valuer < 0) LinF_Valuer = 0;
			if (LinF_Valuer > 255) LinF_Valuer = 255;
			if (LinF_Valueg < 0) LinF_Valueg = 0;
			if (LinF_Valueg > 255) LinF_Valueg = 255;
			if (LinF_Valueb < 0) LinF_Valueb = 0;
			if (LinF_Valueb > 255) LinF_Valueb = 255;
			rgb2[y][x].rgbRed = (int)LinF_Valuer;
			rgb2[y][x].rgbGreen = (int)LinF_Valueg;
			rgb2[y][x].rgbBlue = (int)LinF_Valueb;
		}
	}
	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}
void formpictures(int ksize, RGBQUAD**& rgb, int width, int height, int countofthreads, BITMAPFILEHEADER header, BITMAPINFOHEADER bmiHeader)
{
	double* time = new double[20];
	RGBQUAD** rgbres = new RGBQUAD * [height];
	for (int i = 0; i < height; i++)
		rgbres[i] = new RGBQUAD[width];
	if (countofthreads == 1) {
		for (int i = 0; i < 20; i++)
		{
			time[i] = linefilterposl(ksize, height, width, rgb, rgbres) * 1000;
		}
		double amvTime = ArmV(time);//Среднее арифметическое значение
		amvTime = AvgTrustedInterval(amvTime, time);//Среднее арифметическое значение в доверительном интервале
		string name = "c:\\test\\linePosled" + to_string(ksize) + "_" + to_string(width) + "x" + to_string(height) + ".bmp";
		BMPWrite(rgbres, header, bmiHeader, name.c_str());
		print(countofthreads, ksize, amvTime, ". Линейный фильтр");
		for (int i = 0; i < 20; i++)
		{
			time[i] = gaussposl(ksize, height, width, rgb, rgbres) * 1000;
		}
		amvTime = ArmV(time);//Среднее арифметическое значение
		amvTime = AvgTrustedInterval(amvTime, time);//Среднее арифметическое значение в доверительном интервале
		name = "c:\\test\\gaussPosled" + to_string(ksize) + "_" + to_string(width) + "x" + to_string(height) + ".bmp";
		BMPWrite(rgbres, header, bmiHeader, name.c_str());
		print(countofthreads, ksize, amvTime, ". Фильтр Гаусса");

	}
	else
	{
		omp_set_num_threads(countofthreads);
		for (int i = 0; i < 20; i++)
		{
			time[i] = linefilter(ksize, height, width, rgb, rgbres) * 1000;
		}
		double amvTime = ArmV(time);//Среднее арифметическое значение
		amvTime = AvgTrustedInterval(amvTime, time);//Среднее арифметическое значение в доверительном интервале
		string name = "c:\\test\\lineParallel" + to_string(ksize) + "_" + to_string(countofthreads) +
			"_" + to_string(width) + "x" + to_string(height) + ".bmp";
		BMPWrite(rgbres, header, bmiHeader, name.c_str());
		print(countofthreads, ksize, amvTime, ". Линейный фильтр");
		for (int i = 0; i < 20; i++)
		{
			time[i] = gaussfilter(ksize, height, width, rgb, rgbres) * 1000;
		}
		amvTime = ArmV(time);//Среднее арифметическое значение
		amvTime = AvgTrustedInterval(amvTime, time);//Среднее арифметическое значение в доверительном интервале
		name = "c:\\test\\gaussparallel" + to_string(ksize) + "_" + to_string(countofthreads)
			+ "_" + to_string(width) + "x" + to_string(height) + ".bmp";
		BMPWrite(rgbres, header, bmiHeader, name.c_str());
		print(countofthreads, ksize, amvTime, ". Фильтр Гаусса");
	}
}
int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	double* time = new double[20];
	RGBQUAD** rgb1;
	BITMAPFILEHEADER header1;
	BITMAPINFOHEADER bmiHeader1;
	BMPRead(rgb1, header1, bmiHeader1, "c:\\test\\640x480.bmp");
	RGBQUAD** rgb2;
	BITMAPFILEHEADER header2;
	BITMAPINFOHEADER bmiHeader2;
	BMPRead(rgb2, header2, bmiHeader2, "c:\\test\\1280x800.bmp");
	RGBQUAD** rgb3;
	BITMAPFILEHEADER header3;
	BITMAPINFOHEADER bmiHeader3;
	BMPRead(rgb3, header3, bmiHeader3, "c:\\test\\1680x1050.bmp");
	RGBQUAD** rgb4;
	BITMAPFILEHEADER header4;
	BITMAPINFOHEADER bmiHeader4;
	BMPRead(rgb4, header4, bmiHeader4, "c:\\test\\1920x1200.bmp");
	RGBQUAD** rgbs[] = { rgb1, rgb2, rgb3, rgb4 };
	BITMAPFILEHEADER headers[] = { header1, header2,  header3,  header4 };
	BITMAPINFOHEADER bmiHeaders[] = { bmiHeader1, bmiHeader2, bmiHeader3, bmiHeader4 };
	string name[4]{ "640x480", "1280x800", "1680x1050", "1920x1200" };
	int ksize[3]{ 9,  15, 21 };
	for (int i = 0; i < 4; i++) {
		cout << name[i] << endl;
		for (int j = 0; j < 3; j++)
		{
			formpictures(ksize[j], rgbs[i], bmiHeaders[i].biWidth, bmiHeaders[i].biHeight, 1, headers[i], bmiHeaders[i]);
			formpictures(ksize[j], rgbs[i], bmiHeaders[i].biWidth, bmiHeaders[i].biHeight, 2, headers[i], bmiHeaders[i]);
			formpictures(ksize[j], rgbs[i], bmiHeaders[i].biWidth, bmiHeaders[i].biHeight, 3, headers[i], bmiHeaders[i]);
			formpictures(ksize[j], rgbs[i], bmiHeaders[i].biWidth, bmiHeaders[i].biHeight, 4, headers[i], bmiHeaders[i]);
		}
	}
	system("pause");
	return 0;
}