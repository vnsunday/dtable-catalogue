#include "declarations.hpp"
#include <dynsocc/fundamental/stdex/algorithm.hpp>
#include <dynsocc/fundamental/stdex/tree.h>
#include <cmath>

double GammaFunction(double x) {
    double dt = 0.001;
	double f = 0.0;
	double x1 = x - 1;
	for (double t=0.0; t <= x; t++) {
		f += std::pow(t, x1) * std::exp(-t)*dt;
	}
	return f;
}

double BetaFunction(double z1, double z2) {
	double f = 0.0;
	double dt = 0.001;
	double z11 = z1 - 1;
	double z21 = z2 - 1;

	for (double t=0.0; t<1.0; t++) {
		f += std::pow(t, z11) * std::pow(1 - t, z21) * dt;
	}
	return f;
}

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

