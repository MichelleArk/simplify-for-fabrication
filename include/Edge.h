 #ifndef EDGES_H
 #define EDGES_H
 #include <Eigen/Core>
class RegionEdge {
public:
	RegionEdge(int i, int j, double normal_weight, double area_weight, double perimeter_weight);

	bool operator < (const RegionEdge &other) const {
		return weight < other.weight;
	}

	int region_i;
	int region_j;
	double normal_weight;
	double area_weight;
	double perimeter_weight;
	double weight;
};
 #endif
