#ifndef DTABLE_DATALOGUE_DECLARATION_H_
#define DTABLE_DATALOGUE_DECLARATION_H_

#include <string>
#include <vector>
#include <map>

int MAX_N 10000

#define LB_SUM "sum"
#define LB_ELEMENT "element"
#define LB_RANGE "range"
#define LB_FUNCTION "function"
#define LB_NUMBER "number"
#define LB_CONSTANT "constant" // Number - exchange value 
#define LB_VARIABLE "variable"


#define RL_PROPORTION "propotion"
#define RL_INVPROPORTION "inverse propotion"

#define LB_M_POINT_EST "point estimation"
#define LB_M_INTERVAL_EST "interval estimation"

#define LB_STATISTIC "statistic" // A calculation from data 
#define LB_SAMPLE "sample"
#define LB_DISCRETE "discrete"
#define LB_CONTINUOUS "continuous"

unsigned char component_variable;   // 
unsigned char integral;             // Integral 

std::map<std::string, std::string> mpGraph;

/* Setup everything */
int Initialize();
/*=====  =====*/
std::string azMethod[MAX_N];
int nMethod;

std::string azContrastL[MAX_N];
std::string azContrastR[MAX_N];
std::string azContrast[MAX_N];
int nConstrast;

#endif