#pragma once
class Kmean
{
public:
	Kmean();
	~Kmean();

	double dist(MatrixXd x1, MatrixXd x2);
	MatrixXd randcent(MatrixXd data, int k);
	VectorXi kmeans(MatrixXd data, int k);
};

