#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

int read_standard_distribution(const char* szDistributionFile) {

	// Read file Data 
	double aX[5000];
	double aY[5000];
	double* aBuffer[2]{ &aX[0], &aY[0] };
	
	int n = 0;
	int i;
	int j;

	// Line processing
	int nl = 0; // Line length
	char chl;  // One char in line
	char szbuff[100];
	int nbuff = 0;
	bool bNumber;

	std::ifstream inputFile(szDistributionFile); // Replace "example.txt" with your file path
	double* pData;
	int nData;

	if (inputFile.is_open()) {
		std::string line;
		char chsep = ';';
		nbuff = 0;
		
		while (n < 2 && std::getline(inputFile, line)) {
			/*==================================================
			 * Process one line
			 * PL01. Split Number (or text) by Semicolon (;)
			 * PL02. State-Machine 
			 *		Process char-by-char
			 * 		+ szAccumulated 
			 *		+ Separated character 
			 *		+ End of line 
			 *==================================================*/

			// Split by Semi-Colon
			i = 0;
			nl = line.length();
			bNumber = true;
            nbuff = 0;
			pData = aBuffer[n++]; // Increase n
			nData = 0;

			for (i = 0; i < nl; ++i) {
				chl = line.at(i);

				if (chl == chsep) {
					// Process a Number
					szbuff[nbuff] = NULL; // End string

					if (bNumber && nbuff>0) {
						// Process a Number 
					    pData[nData++] = atof(szbuff);
					}
					else {
						// Ignore if szAccumulated is not a Number
                        bool bNotNumber = true;
					}

					// Reset State 
					bNumber = true;
					nbuff = 0;
				} 
				else {
					szbuff[nbuff++] = chl;
					bNumber &= (chl == '.' || chl == '-' || chl == 'e' || chl == 'E' || (chl >= '0' && chl <= '9'));  // Current Buffer is Number?
				}
			}

			// End of line -- still un processed 
			if (nbuff > 0 && bNumber) {
                szbuff[nbuff] = NULL; // End string
				pData[nData++] = atof(szbuff);
			}
		}
		inputFile.close();

		printf("Finished. n=%d; nData=%d\r\n", n, nData);

		double dTotalProb = 0.0;
		for (int i = 1; i < nData; ++i) {
			dTotalProb += (aX[i] - aX[i - 1])* aY[i];
		}

        // Found outliner
		printf("TotalProbability=%0.5f\r\n", dTotalProb);

		if (1.05 >= dTotalProb && dTotalProb > 0.98) {
			printf("ACCEPTTEDDDDD\r\n");
		}
		else {
			printf("ERRORRRRRRRRRRRRRRRRR\r\n");            
		}
	}
	else {
		std::cerr << "Error: Unable to open file!" << std::endl;
		return 1;
	}

	return 0;
}

double quantile(double *X, double *Y, double p) {
    double d1;
}

int main()
{
    // 
    double X[10000];
    double Y[10000];
    int n = 0;

    

    return 0;
}