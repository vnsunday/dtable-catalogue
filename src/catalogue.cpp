#include <conio.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

using namespace std;
#define PI 3.14159265358979323846

#define FILE_BINOMIAL_8_0_5 "C:\\Temp\\Binomial_8_0.5.csv"

// fs : std::ofstream
#define APPEND_ARRAY_TOFILE(fs, arr, n, separator) { for (int j = 0; j < n; j++) { \
    if (j > 0) {\
        fs << separator;\
    }\
        fs << arr[j];\
    }\
        fs << std::endl;  }\

#define STRING_TO_INT_ARRAY(fs, arr, n, separator) {}
#define STRING_TO_DOUBLE_ARRAY(fs, str, arr, n, separator) {

int detect_outliers_zscore(double* pA, int nA, int& nNormal, int& nOutlier, vector<int>& voutlier) {
    // Made a distribution
    // z-score 
    // Three sigma
    double dSum = 0;
    double dMean; 
    double dVariance = 0.0;
    double dSigma;
    double dThreeSigmaL;
    double dThreeSigmaH;


    dSum = std::accumulate(pA, pA + nA, 0);
    dMean = dSum / nA;

    for (int i = 0; i < nA; ++i) {
        double dVal = (pA[i] - dMean);
        dVariance += dVal * dVal;
    }

    dVariance = dVariance / nA;
    dSigma = sqrt(dVariance);

    dThreeSigmaL = dMean - 3 * dSigma;
    dThreeSigmaH = dMean + 3 * dSigma;

    nNormal = 0;
    nOutlier = 0;
    voutlier.clear();

    for (int i = 0; i < nA; ++i) {
        if (pA[i] >= dThreeSigmaL && pA[i] <= dThreeSigmaH) {
            nNormal++;
        }
        else {
            nOutlier++;
            voutlier.push_back(i);
        }
    }

    return 0;
}
            
int string_to_int_arrary(std::string str, char chseparator, int&n, int* pI) {

    int nl = str.length();
    char chl;
    char szbuff[1000];
    int nbuff = 0;
    bool bNumber = true;
    
    n = 0;
    for (int i = 0; i < nl; ++i) {
        chl = str.at(i);

        if (chl == chseparator) {
            // Process a Number
            szbuff[nbuff++] = NULL; // End string

            if (bNumber && nbuff > 0) {
                // Process a Number 
                pI[n++] = atoi(szbuff);
            }
            else {
                // Ignore if szAccumulated is not a Number
            }

            // Reset State 
            bNumber = true;
            nbuff = 0;
        }
        else {
            szbuff[nbuff++] = chl;
            bNumber &= (chl == '.' || (chl >= '0' && chl <= '9'));  // Current Buffer is Number?
        }
    }

    // End of line -- still un processed 
    if (nbuff > 0 && bNumber) {
        pI[n++] = atoi(szbuff);
    }
}

int string_to_int_array(std::string str, char chseparator, int&n, int* pI) {
    int nl = str.length();
    char chl;
    char szbuff[1000];
    int nbuff = 0;
    bool bNumber = true;

    n = 0;

    for (int i = 0; i < nl; ++i) {
        chl = str.at(i);

        if (chl == chseparator) {
            // Process a Number
            szbuff[nbuff++] = NULL; // End string

            if (bNumber && nbuff > 0) {
                // Process a Number 
                pI[n++] = atoi(szbuff);
            }
            else {
                // Ignore if szAccumulated is not a Number
            }

            // Reset State 
            bNumber = true;
            nbuff = 0;
        }
        else {
            szbuff[nbuff++] = chl;
            bNumber &= ((chl >= '0' && chl <= '9'));  // Current Buffer is Number?
        }
    }

    // End of line -- still un processed 
    if (nbuff > 0 && bNumber) {
        pI[n++] = atoi(szbuff);
    }
}

int string_to_double_arrary(std::string str, char chseparator, int&n, double* pD) {
    int nl = str.length();
    char chl;
    char szbuff[1000];
    int nbuff = 0;
    bool bNumber = true;

    n = 0;

    for (int i = 0; i < nl; ++i) {
        chl = str.at(i);

        if (chl == chseparator) {
            // Process a Number
            szbuff[nbuff++] = NULL; // End string

            if (bNumber && nbuff > 0) {
                // Process a Number 
                pD[n++] = atof(szbuff);
            }
            else {
                // Ignore if szAccumulated is not a Number
            }

            // Reset State 
            bNumber = true;
            nbuff = 0;
        }
        else {
            szbuff[nbuff++] = chl;
            bNumber &= (chl == '.' || (chl >= '0' && chl <= '9'));  // Current Buffer is Number?
        }
    }

    // End of line -- still un processed 
    if (nbuff > 0 && bNumber) {
        pD[n++] = atof(szbuff);
    }
}

#ifdef DYNSOCCLIB

// Supported Text View
int int_to_string(int* pA, int n, std::vector<std::string>& v, int colwidth, int align=DYNSOCC_TEXT_ALIGN_RIGHT) {
    char szBuffer[100];
    char szNumber[100];
    int nNumstrLen;
    int nBuffer;

    char chspace = ' ';
    

    if (align == DYNSOCC_TEXT_ALIGN_RIGHT) {
        for (int i = 0; i < n; ++i) {
            sprintf(szNumber, "%d", pA[i]);
            nNumstrLen = strlen(szNumber);

            // Empty string 
            memset(szBuffer, chspace, colwidth);
            szBuffer[colwidth] = NULL;

            if (nNumstrLen <= colwidth) {
                memcpy(szBuffer + colwidth - nNumstrLen, szNumber, nNumstrLen);
            }
            else {
                memcpy(szBuffer, szNumber, colwidth);
            }

            std::string strCell = szBuffer;
            v.push_back(strCell);
        }
    }
    else if (align == DYNSOCC_TEXT_ALIGN_LEFT) {
        for (int i = 0; i < n; ++i) {
            sprintf(szNumber, "%d", pA[i]);
            nNumstrLen = strlen(szNumber);

            // Empty string 
            memset(szBuffer, chspace, colwidth);
            szBuffer[colwidth] = NULL;

            if (nNumstrLen <= colwidth) {
                memcpy(szBuffer, szNumber, nNumstrLen);
            }
            else {
                memcpy(szBuffer, szNumber, colwidth);
            }

            std::string strCell = szBuffer;
            v.push_back(strCell);
        }
    }
    else {
        throw "UnImplemented";
    }

    return 0;
}

int double_to_string(double* pD, int n, std::vector<std::string>& v, int colwidth, int align = DYNSOCC_TEXT_ALIGN_RIGHT) {
    char szBuffer[100];
    char szNumber[100];
    int nNumstrLen;
    int nBuffer;

    char chspace = ' ';


    if (align == DYNSOCC_TEXT_ALIGN_RIGHT) {
        for (int i = 0; i < n; ++i) {
            sprintf(szNumber, "%0.3f", pD[i]);
            nNumstrLen = strlen(szNumber);

            // Empty string 
            memset(szBuffer, chspace, colwidth);
            szBuffer[colwidth] = NULL;

            if (nNumstrLen <= colwidth) {
                memcpy(szBuffer + colwidth - nNumstrLen, szNumber, nNumstrLen);
            }
            else {
                memcpy(szBuffer, szNumber, colwidth);
            }

            std::string strCell = szBuffer;
            v.push_back(strCell);
        }
    }
    else if (align == DYNSOCC_TEXT_ALIGN_LEFT) {
        for (int i = 0; i < n; ++i) {
            sprintf(szNumber, "%0.3f", pD[i]);
            nNumstrLen = strlen(szNumber);

            // Empty string 
            memset(szBuffer, chspace, colwidth);
            szBuffer[colwidth] = NULL;

            if (nNumstrLen <= colwidth) {
                memcpy(szBuffer, szNumber, nNumstrLen);
            }
            else {
                memcpy(szBuffer, szNumber, colwidth);
            }

            std::string strCell = szBuffer;
            v.push_back(strCell);
        }
    }
    else {
        throw "UnImplemented";
    }

    return 0;
}

#endif

// End of Supported Text View
int test1()
{
    const int n = 400; // Sample size
    const int k = 5;   // Category size
    double pEqual = 1/(double)k;
    double O[k] { 95, 65, 60, 80, 100 }; // Observations
    double p[k] { pEqual, pEqual, pEqual, pEqual, pEqual };  // Null Hypothesis
    double E[k]; 
    double significance_level = 0.01; // 1%
    double Xp2;
    int df = k-1;

    for (int i=0; i<k; ++i) {
        E[i] = p[i] * n;

        if (E[i] < 5) {
            printf("Data not enough\r\n");
            return -1;
        }
    }

    Xp2 = 0;
    for (int i=0; i<k; i++) {
        Xp2 += (O[i] - E[i])*(O[i] - E[i]) / E[i];
    }

    printf("X^2=%0.3f\r\n", Xp2);
    return 0;
}

int one_sample_t_test() {

    double x[100];
    double M; // mean
    double M0; // Specific value 
    double test_statistic; 
    
    return 0;
}

int critical_value_normal_distribution() {
	double d1 = 0;
	double crv = 0.05; // Critical-value

	/*================================================================================*
	    Normal-Distribution
			N(0,1)=
			Normal(m,s^2) = (1 /sqrt(2pi)) exp (-x^2/2)
	 *================================================================================*/
	double d_bin = 2000;
	double XL = -1.0;
	double XR = 1.0;
	double aX[5000];
	double aY[5000];
	int i = 0;

	double xStep = (XR - XL) / d_bin;
	double x = XL;
	while (x <= XR) {
		aX[i] = x;
		aY[i] = exp(-x * x / 2.0) / sqrt(2 * PI);
		x += xStep;
		i++;
	}

	std::ofstream outputFile("C:\\Temp\\StandardNormal.csv", std::ios::app);
	// Check if the file was successfully opened
	if (outputFile.is_open()) {
		for (int j = 0; j < i; j++) {
			if (j > 0) {
				outputFile << ";";
			}
			outputFile << aX[j];
		}
		outputFile << std::endl;
		for (int j = 0; j < i; j++) {
			if (j > 0) {
				outputFile << ";";
			}
			outputFile << aY[j];
		}
		outputFile.close();
	}
	else {
		// Handle the case where the file could not be opened
		std::cerr << "Error: Unable to open file for appending.\n";
	}

	return 0;
}

int build_normal_distribution_table(const char* szFile, double xL=-3, double xR=3, int nNumBin=1000) {
	// double XL = -10.0;
	// double XR = 10.0;
	double aX[5000];
	double aY[5000];
	int n = 0;

	double xStep = (xR - xL) / nNumBin;

	double x = xL;
	while (x <= xR) {
		aX[n] = x;
		aY[n] = exp(-x * x / 2.0) / sqrt(2 * PI);
		x += xStep;
		n++;
	}

	std::ofstream outputFile(szFile, std::ios::app);

	// Check if the file was successfully opened
	if (outputFile.is_open()) {
		for (int j = 0; j < n; j++) {
			if (j > 0) {
				outputFile << ";";
			}
			outputFile << aX[j];
		}
		outputFile << std::endl;
		for (int j = 0; j < n; j++) {
			if (j > 0) {
				outputFile << ";";
			}
			outputFile << aY[j];
		}
		outputFile.close();
	}
	else {
		// Handle the case where the file could not be opened
		std::cerr << "Error: Unable to open file for appending.\n";
	}
	
	return 0;
}

int assess_standard_distribution(const char* szDistributionFile) {

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

        int nNormal;
        int nOutlier;
        vector<int> vOutlier;
        detect_outliers_zscore(aY, nData, nNormal, nOutlier, vOutlier);

        printf("Detect Outlier (aY): nNormal=%d; nOutlier=%d; Outlier Values= ", nNormal, nOutlier);

        for (int i = 0; i < vOutlier.size(); ++i) {
            printf("%0.5f;", aY[vOutlier[i]]);
        }
        printf("\r\n");

        vOutlier.clear();
        detect_outliers_zscore(aX, nData, nNormal, nOutlier, vOutlier);
        printf("Detect Outlier (aX): nNormal=%d; nOutlier=%d; Outlier Values= ", nNormal, nOutlier);

        for (int i = 0; i < vOutlier.size(); ++i) {
            printf("[%d]=%0.5f;", vOutlier[i], aX[vOutlier[i]]);
        }

        printf("\r\n");
	}
	else {
		std::cerr << "Error: Unable to open file!" << std::endl;
		return 1;
	}

	return 0;
}

int print_binomial_table(int n, double p, const char* szFile) {
    int aX[100];
    double aY[100];

    if (n > 100) {
        printf("Invalid\r\n");
        return -1;
    }

    // Binomial Table  (8, 0.5)
    // double p = 0.5;
    int choiceNK = 1;
    double pow_k = 1;
    double pow_n_k = 1;
    int nX = 0;

    for (int i = 0; i < n; i++) {
        pow_n_k = pow_n_k * (1 - p);
    }
    for (int k = 0; k <= n; k++) {
        // 
        // P(X=k) = (n k) p^k (1-p)^(n-k)
        // 

        if (k > 0) {
            // Calculating C(n,k)
            choiceNK = choiceNK * (n - k + 1) / (k);
            pow_k = pow_k * p;
            pow_n_k = pow_n_k / (1 - p);
        }

        double pX = choiceNK * pow_k * pow_n_k;
        
        aX[nX] = k;
        aY[nX++] = pX;
    }

    /*================================================
     *  Print to File.
     *================================================*/
    std::ofstream outputFile(szFile, std::ios::app);
    if (outputFile.is_open()) {
        APPEND_ARRAY_TOFILE(outputFile, aX, n+1, ";");
        APPEND_ARRAY_TOFILE(outputFile, aY, n+1, ";");
        outputFile.close();
    }
    else {
        // Handle the case where the file could not be opened
        std::cerr << "Error: Unable to open file for appending.\n";
    }

    return 0;
}

int read_binomial_table(const char* szFile, int& n, int* paX, double* paY) {
    std::ifstream fs(szFile); // Replace "example.txt" with your file path
    double* pData;
    int nData;
    string strLine;
    bool bX;

    if (fs.is_open()) {
        int nLine = 0;

        bX = false;
        if (std::getline(fs, strLine)) {
            string_to_int_arrary(strLine, ';', nData, paX);
            bX = true;
        }

        if (bX && std::getline(fs, strLine)) {
            string_to_double_arrary(strLine, ';', nData, paY);
        }

        fs.close();
    }
    else {
        // Error
        printf("Can't open file %s\r\n", szFile);
    }
    return 0;
}

int example_chi_square_goodness_fit() {

    /*==================================================
        Observed and Expected and count
        1. Table of expected
        2. Degree of freedom
        3. Likelihood ratio statistic G
        4. Likelihood ratio statistic G
     *==================================================*/


    int aObsX[6] = { 0, 1, 2, 3, 4, 5 }; // Observed Outcomes
    int aObsY[6] = { 3, 10, 15, 13, 7, 3 }; // Observed Count
    int nObs = 6;
    double aExp[6]; // Expected 
    double aP[6]; // Possibility of outcomes 

    double aExpected[9];
    int n = std::accumulate(aObsY, aObsY + 6, 0); // Count of samples

    int pBinomX[9];
    double pBinomY[9];
    int nBinom;

    read_binomial_table(FILE_BINOMIAL_8_0_5, nBinom, &pBinomX[0], &pBinomY[0]);
    
    for (int i = 0; i < 5; i++) {
        aP[i] = pBinomY[i];
    }
    aP[5] = pBinomY[5] + pBinomY[6] + pBinomY[7] + pBinomY[8]; // X>=5
 
    for (int i = 0; i < 6; i++) {
        aExp[i] = aP[i] * n;
    }

    /*
    std::vector<int> vcolsize(7, 8); // 8 column. WIdth 8
    std::vector<int> vformat(7, DYNSOCC_TEXT_ALIGN_RIGHT); // 8 column. Right Format.
    vformat[0] = DYNSOCC_TEXT_ALIGN_LEFT;
    vcolsize[0] = 34;

        
    vector<string> vcol;
    vector<string> vObs;
    vector<string> vPBino;
    vector<string> vExp;
    vector<vector<string>> vtable;
    
    vcol.push_back("Output");
    int_to_string(aObsX, 5, vcol, 8);
    vcol.push_back(">=5");

    vObs.push_back("Observed");
    int_to_string(aObsY, 6, vObs, 8);

    vPBino.push_back("H0-Probobability Binomal(8,0.5)");
    double_to_string(aP, 6, vPBino, 8);

    vExp.push_back("Expected");
    double_to_string(aExp, 6, vExp, 8);

    vtable.push_back(vcol);
    vtable.push_back(vObs);
    vtable.push_back(vPBino);
    vtable.push_back(vExp);
    */

    // std::cout << dynsocc::TextView::table_print(vcolsize, vformat, vtable, '|') << std::endl;

    return 0;
}

int main(int argc, char const *argv[])
{
    if (argc < 2) {
        printf(
            "Usage: "
            "    catalogue distrubition [create|load]");
    }
    else {
        /*============================================================*/
        return 0;
    }
    
    test1();
    printf("Chi-square Goodness of fit Test\r\n");
    // return assess_standard_distribution("C:\\Temp\\NormalDistribution_Range3_Bin1000.csv");
    // return build_normal_distribution_table("C:\\Temp\\NormalDistribution_Range3_Bin1000.csv", -3.0, 3.0, 1000);
    // return assess_standard_distribution("C:\\Temp\\NormalDistribution.csv");
    
    return assess_standard_distribution("C:\\Temp\\NormalDistribution_Range_10.csv");
    // return build_normal_distribution_table("C:\\Temp\\NormalDistribution_Range_10.csv");
    // return example_chi_square_goodness_fit();
    // return 0;

    return print_binomial_table(8, 0.5, FILE_BINOMIAL_8_0_5);
	return assess_standard_distribution("C:\\Temp\\NormalDistribution.csv");
	return critical_value_normal_distribution();

	// Choice(N,k) = Choice(N,k+1) * 
	// Double. P(X>=5)
    return 0;
}