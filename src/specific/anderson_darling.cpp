#include <stdio.h>
#include <algorithm>
#include <cmath>

const double TWO_SQRT = sqrt(2);

double CDF_lognormal(double x, double sigma, double mean)
{
	return (1 + erf(  (std::log(x) - mean)/(sigma * TWO_SQRT)  )) / 2;
}

int anderson_darling_normal_test()
{
	double arr_data[] = { 1.5112,0.2129,0.2824,9.8631,3.3056,0.1673,0.3722,0.1027,0.1076,0.1517,0.4820,0.4490,0.9259,0.6246,0.3693,11.7244,0.1076,0.2382,0.0900,0.1144,0.2832,0.4416 };
	int N = sizeof(arr_data) / sizeof(arr_data[0]);
	double AD = 0.0;

	double mean = 0.0;
	double sigma = 0.0;

	/*
		1. SORT Data 
	 */

	for (int i = 0; i < N; ++i)
	{
		mean += arr_data[i];
	}
	mean = mean / N; 

	for (int i = 0; i < N; ++i)
	{
		sigma += (arr_data[i] - mean)*(arr_data[i] - mean);
	}
	sigma = sigma / (N - 1);

	return 0;
}

int main(int argc, char const *argv[])
{
	double arr_data[] = { 1.5112,0.2129,0.2824,9.8631,3.3056,0.1673,0.3722,0.1027,0.1076,0.1517,0.4820,0.4490,0.9259,0.6246,0.3693,11.7244,0.1076,0.2382,0.0900,0.1144,0.2832,0.4416  };
	int N = sizeof(arr_data) / sizeof(arr_data[0]);
	double AD = 0.0;

	double average = 0.0;
	double sigma = 0.0;

	/*============================================================
		1. Sort Data
		2. Calculating AD value 
	 *============================================================*/
	//1. Sort
	std::sort(&arr_data[0], &arr_data[N]);
	for (int i = 0; i < N - 1; ++i)
	{
		if (arr_data[i] > arr_data[i + 1])
		{
			printf("Error\r\n");
		}
	}
	printf("Sort-Verified\r\n");

	//2.AD Value 
	average = 0.0;
	sigma = 0.0;

	for (int i = 0; i < N; ++i)
	{
		average += arr_data[i];
	}
	average = average / N;

	for (int i = 0; i < N; ++i)
	{
		sigma += (arr_data[i] - average)*(arr_data[i] - average) / (N-1);
	}
	sigma = sqrt(sigma);


	AD = -N;
	for (int i = 1; i <= N; ++i)
	{
		AD = AD - 1 / (double)N * (
				(2*i-1) * 
				( std::log(CDF_lognormal(arr_data[i-1], sigma, average)) + 
					std::log(1-CDF_lognormal(arr_data[N-i], sigma, average))
				)
			);
	}

	// AD = 0.9491456;
	double p = std::exp(1.2937 - 5.709*AD)+0.0186*(AD*AD);

    printf("Anderson Darling. AD=%f\r\n. p = %f\r\n", AD, p);
    return 0;
}