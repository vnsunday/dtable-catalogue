#ifndef DTABLE_DATALOGUE_DECLARATION_H_
#define DTABLE_DATALOGUE_DECLARATION_H_

#include <string>
#include <vector>
#include <map>


#define LB_SUM "sum"
#define LB_ELEMENT "element"
#define LB_RANGE "range"
#define LB_FUNCTION "function"
#define LB_NUMBER "number"
#define LB_CONSTANT "constant" // Number - exchange value 
#define LB_VARIABLE "variable" 

#define LB_STATISTIC "statistic" // A calculation from data 
#define LB_SAMPLE "sample"
#define LB_DISCRETE "discrete"
#define LB_CONTINUOUS "continuous"

unsigned char component_variable;   // 
unsigned char integral;             // Integral 

std::map<std::string, std::string> mpGraph;

/* Setup everything */
int Initialize();

#endif