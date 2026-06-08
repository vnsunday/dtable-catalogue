#include "declarations.hpp"
#include <dynsocc/fundamental/stdex/algorithm.hpp>
#include <dynsocc/fundamental/stdex/tree.h>

int Initialize() {
	mpGraph["name"] = "dtable-catalogue";

	nMethod = 0;
	azMethod[nMethod++] = LB_M_POINT_EST;
	azMethod[nMethod++] = LB_M_INTERVAL_EST;

	
	dynsocc::TreeAdj tree;

	int nRoot;
	int nNewID;

	tree.add_node(-1, "root", nRoot);
	tree.add_node(nRoot, "confidence_interval", nNewID);
	

	return 0;
}

