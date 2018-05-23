#include "stdafx.h"
#include "Kmean.h"


Kmean::Kmean()
{
}


Kmean::~Kmean()
{
}


double Kmean::dist(MatrixXd x1, MatrixXd x2)
{
	double dict = 0;
	for (int i=0;i < x1.rows(); i++)
	{
		for (int j = 0; j < x1.cols(); j++)
		{
			dict += abs(pow(x1(i, j) - x2(i, j), 2));
		}
	}
	return dict;
}


MatrixXd Kmean::randcent(MatrixXd data,int k)
{
	MatrixXd centroids(k,data.cols());
	MatrixXd min = data.colwise().minCoeff();
	MatrixXd max = data.colwise().maxCoeff();
	for (int j = 0; j < data.cols(); j++)
	{
		double minj = min(0, j);
		double maxj = max(0, j);
		double range = maxj - minj;
		for (int i = 0; i < k; i++)
		{
			centroids(i, j) = minj + range*(rand() % 1000)*0.001;
		}
	}
	//cout << data.colwise().minCoeff() << endl;
	return centroids;
}


//MatrixXd Kmean::kmeans(MatrixXd data,int k)
//{
//	int m = data.rows();
//	VectorXi cluster(m);
//	cluster.setZero(m);
//	//cout << data.row(0) << endl;
//	//cout << cluster << endl;
//	MatrixXd centroids = randcent(data, k);
//	bool clusterchanged = true;
//	//cout << centroids << endl;
//	while (clusterchanged)
//	{
//		clusterchanged = false;
//		for (int i = 0; i < m; i++)
//		{
//			double mindist = MAX_VALUE;
//			int minindex = 0;
//			for (int j = 0; j < k; j++)
//			{
//				double distji = dist(centroids.row(j), data.row(i));
//				if (distji < mindist)
//				{
//					mindist = distji;
//					minindex = j;
//				}
//			}
//			if (cluster(i) != minindex)
//			{
//				clusterchanged = true;
//				cluster(i) = minindex;
//			}
//		}
//	//	cout << cluster << endl;
//		cout << centroids << endl;
//		cout << "-----------" << endl;
//		for (int i = 0; i < k; i++)
//		{
//			MatrixXd sum= MatrixXd::Zero(1, data.cols());
//			int n = 0;
//			for (int j = 0; j < m; j++)
//			{
//				if (cluster(j) == i)
//				{
//					//cout << "yes" << endl;
//					//cout << data.row(j) << endl;
//					sum+= data.row(j);
//					n++;
//				}
//			}
//			//cout << sum << endl;
//			if (n > 0)
//			{
//				centroids.row(i) = sum / n;
//			}
//			else
//			{
//				//centroids.row(k) = sum / 0.001;
//			}
//					
//		}
//
//	}
//	return centroids;
//}

//MatrixXd Kmean::kmeans(MatrixXd data, int k)
//{
//	int m = data.rows();
//	int N = 2;
//	//cout << data.row(0) << endl;
//	//cout << cluster << endl;
//	MatrixXd centroids = randcent(data, k);
//	bool clusterchanged = true;
//	//cout << centroids << endl;
//	while (clusterchanged)
//	{
//		clusterchanged = false;
//		VectorXi datastate(m);
//		datastate.setOnes(m);
//		VectorXi cluster(m);
//		cluster.setZero(m);
//		for (int j = 0; j < k; j++)
//		{
//			//MatrixXd tmp_data = data;
//			for (int n = 0; n < N; n++)
//			{
//				double mindist = MAX_VALUE;
//				int minindex = 0;
//				for (int i = 0; i < m; i++)
//				{
//					if (datastate(i))
//					{
//						double distji = dist(centroids.row(j), data.row(i));
//						if (distji < mindist)
//						{
//							mindist = distji;
//							minindex = i;
//						}
//					}
//					
//				}
//				
//					cluster(minindex) = j;
//					datastate(minindex) = 0;
//				
//				
//			}
//		}
//		cout << cluster << endl;
//		cout << centroids << endl;
//		cout << "-----------" << endl;
//		for (int i = 0; i < k; i++)
//		{
//			MatrixXd sum = MatrixXd::Zero(1, data.cols());
//			int n = 0;
//			for (int j = 0; j < m; j++)
//			{
//				if (cluster(j) == i)
//				{
//					//cout << "yes" << endl;
//					//cout << data.row(j) << endl;
//					sum += data.row(j);
//					n++;
//				}
//			}
//			//cout << sum << endl;
//			if (n > 0)
//			{
//				if (centroids.row(i)!=(sum/n))
//				{
//					clusterchanged = true;
//				}
//				centroids.row(i) = sum / n;
//			}
//			else
//			{
//				//centroids.row(k) = sum / 0.001;
//			}
//
//		}
//
//	}
//	return centroids;
//}

//VectorXi Kmean::kmeans(MatrixXd data,int k)
//{
//	int m = data.rows();
//	VectorXi cluster(m);
//	cluster.setZero(m);
//	MatrixXd centroids = randcent(data, k);
//	MatrixXd ijdist(m, m);
//	VectorXi cancluster(m);
//	cancluster.setOnes(m);
//	std::ptrdiff_t i, j;
//	//cout << data << endl;
//	//data.maxCoeff(&i, &j);
//	//cout << data.row(1).maxCoeff(&i) << endl;
//	//cout << data.col(1).maxCoeff(&j) << endl;
//	//cout << i << "  " << j << endl;
//	for (int I = 0; I < m; I++)
//	{
//		for (int J = 0; J < m; J++)
//		{
//			if (I == J)
//			{
//				ijdist(I, J) = MAX_VALUE;
//			}
//			else
//			{
//				ijdist(I, J) = dist(data.row(I), data.row(J));
//			}
//		}
//	}
//	//cout << ijdist << endl;
//	for (int K = 0; K < k; K++)
//	{
//		double mindist = ijdist.minCoeff(&i, &j);
//		centroids.row(K) = (data.row(i) + data.row(j)) / 2;
//		cluster(i) = K;
//		cluster(j) = K;
//		for (int I = 0; I < m; I++)
//		{
//			ijdist(I, j) = MAX_VALUE;
//			ijdist(I, i) = MAX_VALUE;
//			ijdist(i, I) = MAX_VALUE;
//			ijdist(j, I) = MAX_VALUE;
//
//		}
//
//	}
//	
//	//cout << cluster << endl;
//
//	
//	return cluster;
//}

VectorXi Kmean::kmeans(MatrixXd data, int k)
{
	int m = data.rows();
	VectorXi cluster(m);
	cluster.setZero(m);
	VectorXi cluster_tmp(m);
	cluster_tmp.setZero(m);
	VectorXi pointstate(m);
	pointstate.setOnes(m);
	MatrixXd centroids = randcent(data, k);
	MatrixXd centroids_tmp = randcent(data,8);
	//MatrixXd ijdist(m, m);
	VectorXi cancluster(m);
	cancluster.setOnes(m);
	std::ptrdiff_t i, j;
	for (int K = 0; K < 8; K++)
	{
		double mindist = MAX_VALUE;
		int index1, index2;
		for (int i = 0; i < m; i++)
		{
			for (int j = i + 1; j < m; j++)
			{
				if (pointstate(i) && pointstate(j))
				{
					double ijdist = dist(data.row(i), data.row(j));
					if (ijdist < mindist)
					{
						mindist = ijdist;
						index1 = i;
						index2 = j;
					}
				}
				
			}
		}
		cluster_tmp(index1) = K;
		cluster_tmp(index2) = K;
		pointstate(index1) = 0;
		pointstate(index2) = 0;
		centroids_tmp.row(K) = (data.row(index1) + data.row(index2)) / 2;
	}

	VectorXi pointstate_tmp(m/2);
	pointstate_tmp.setOnes(m/2);
	for (int K = 0; K < 4; K++)
	{
		double mindist = MAX_VALUE;
		int index1, index2;
		for (int i = 0; i < m/2; i++)
		{
			for (int j = i + 1; j < m/2; j++)
			{
				if (pointstate_tmp(i) && pointstate_tmp(j))
				{
					double ijdist = dist(centroids_tmp.row(i), centroids_tmp.row(j));
					if (ijdist < mindist)
					{
						mindist = ijdist;
						index1 = i;
						index2 = j;
					}
				}

			}
		}
		for (int i = 0; i < m; i++)
		{
			if (cluster_tmp(i) == index1 || cluster_tmp(i) == index2)
			{
				cluster(i) = K;
			}
		}
		
		pointstate_tmp(index1) = 0;
		pointstate_tmp(index2) = 0;

		centroids.row(K) = (centroids_tmp.row(index1) + centroids_tmp.row(index2)) / 2;
	}


	return cluster;
}