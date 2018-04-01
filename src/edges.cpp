#include "Edge.h"

RegionEdge::RegionEdge(int i, int j, double normal, double area, double perimeter) {
	region_i = i;
	region_j = j;
	normal_weight = normal;
	area_weight = area;
	perimeter_weight = perimeter;
	weight = area * (normal*normal) * perimeter;
		//area*0.0001 + perimeter) / 3.0; // check later
};
