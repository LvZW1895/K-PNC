#include "stdafx.h"
#include "Demapping.h"
#include <map>
#include "Kmean.h"
#include <set>

Demapping::Demapping(int k)
{
	if (k == QAM16)
	{
		MatrixXd bit16(8, 256);
		Matrix<complex<double>, Dynamic, Dynamic> bit_16qam_tmp(2, 256);

		for (int i = 0; i < 256; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				if (i&(1 << j))
					bit16(7 - j, i) = 1;
				else
					bit16(7 - j, i) = 0;
			}
		}
		bit_16 = bit16;
		MatrixXd ncm16 = bit16.transpose();
		ncm_16qam = ncm16;
		//cout << ncm16 << endl;
		for (int i = 0; i < 256; i++)
		{
			int tmp1, tmp2;
			if (bit_16(0, i) == 0)
			{
				if (bit_16(2, i) == 0)
				{
					tmp1 = 1;
				}
				else
				{
					tmp1 = 3;
				}
			}
			else
			{
				if (bit_16(2, i) == 0)
				{
					tmp1 = -1;
				}
				else
				{
					tmp1 = -3;
				}
			}
			if (bit_16(1, i) == 0)
			{
				if (bit_16(3, i) == 0)
				{
					tmp2 = 1;
				}
				else
				{
					tmp2 = 3;
				}
			}
			else
			{
				if (bit_16(3, i) == 0)
				{
					tmp2 = -1;
				}
				else
				{
					tmp2 = -3;
				}
			}
			complex<double> c1(tmp1, tmp2);
			bit_16qam_tmp(0, i) = c1;
			if (bit_16(4, i) == 0)
			{
				if (bit_16(6, i) == 0)
				{
					tmp1 = 1;
				}
				else
				{
					tmp1 = 3;
				}
			}
			else
			{
				if (bit_16(6, i) == 0)
				{
					tmp1 = -1;
				}
				else
				{
					tmp1 = -3;
				}
			}
			if (bit_16(5, i) == 0)
			{
				if (bit_16(7, i) == 0)
				{
					tmp2 = 1;
				}
				else
				{
					tmp2 = 3;
				}
			}
			else
			{
				if (bit_16(7, i) == 0)
				{
					tmp2 = -1;
				}
				else
				{
					tmp2 = -3;
				}
			}
			complex<double> c2(tmp1, tmp2);
			bit_16qam_tmp(1, i) = c2;
		}
		bit_16qam = bit_16qam_tmp;
	}
	else if(k==QAM64)
	{
		MatrixXd bit64(12, 4096);
		Matrix<complex<double>, Dynamic, Dynamic> bit_64qam_tmp(2, 4096);
		for (int i = 0; i < 4096; i++)
		{
			for (int j = 0; j < 12; j++)
			{
				if (i&(1 << j))
					bit64(11 - j, i) = 1;
				else
					bit64(11 - j, i) = 0;
			}
		}
		bit_64 = bit64;
		for (int i = 0; i < 4096; i++)
		{
			int tmp1, tmp2;
			if (bit64(0,i) == 0)
			{
				if (bit64(1,i) == 0)
				{
					if ((bit64(2,i) == 0) && (bit64(4,i) == 1))
					{
						tmp1 = 1;
					}
					else if ((bit64(2,i) == 0) && (bit64(4,i) == 0))
					{
						tmp1 = 3;
					}
					else if ((bit64(2,i) == 1) && (bit64(4,i) == 0))
					{
						tmp1 = 5;
					}
					else if ((bit64(2,i) == 1) && (bit64(4,i) == 1))
					{
						tmp1 = 7;
					}
					if ((bit64(3,i) == 0) && (bit64(5,i) == 1))
					{
						tmp2 = 1;
					}
					else if ((bit64(3,i) == 0) && (bit64(5,i) == 0))
					{
						tmp2 = 3;
					}
					else if ((bit64(3,i) == 1) && (bit64(5,i) == 0))
					{
						tmp2 = 5;
					}
					else if ((bit64(3,i) == 1) && (bit64(5,i) == 1))
					{
						tmp2 = 7;
					}
				}
				else
				{
					if ((bit64(2,i) == 0) && (bit64(4,i) == 1))
					{
						tmp1 = 1;
					}
					else if ((bit64(2,i) == 0) && (bit64(4,i) == 0))
					{
						tmp1 = 3;
					}
					else if ((bit64(2,i) == 1) && (bit64(4,i) == 0))
					{
						tmp1 = 5;
					}
					else if ((bit64(2,i) == 1) && (bit64(4,i) == 1))
					{
						tmp1 = 7;
					}
					if ((bit64(3,i) == 0) && (bit64(5,i) == 1))
					{
						tmp2 = -1;
					}
					else if ((bit64(3,i) == 0) && (bit64(5,i) == 0))
					{
						tmp2 = -3;
					}
					else if ((bit64(3,i) == 1) && (bit64(5,i) == 0))
					{
						tmp2 = -5;
					}
					else if ((bit64(3,i) == 1) && (bit64(5,i) == 1))
					{
						tmp2 = -7;
					}
				}
			}
			else
			{
				if (bit64(1,i) == 0)
				{
					if ((bit64(2,i) == 0) && (bit64(4,i) == 1))
					{
						tmp1 = -1;
					}
					else if ((bit64(2,i) == 0) && (bit64(4,i) == 0))
					{
						tmp1 = -3;
					}
					else if ((bit64(2,i) == 1) && (bit64(4,i) == 0))
					{
						tmp1 = -5;
					}
					else if ((bit64(2,i) == 1) && (bit64(4,i) == 1))
					{
						tmp1 = -7;
					}
					if ((bit64(3,i) == 0) && (bit64(5,i) == 1))
					{
						tmp2 = 1;
					}
					else if ((bit64(3,i) == 0) && (bit64(5,i) == 0))
					{
						tmp2 = 3;
					}
					else if ((bit64(3,i) == 1) && (bit64(5,i) == 0))
					{
						tmp2 = 5;
					}
					else if ((bit64(3,i) == 1) && (bit64(5,i) == 1))
					{
						tmp2 = 7;
					}
				}
				else
				{
					if ((bit64(2,i) == 0) && (bit64(4,i) == 1))
					{
						tmp1 = -1;
					}
					else if ((bit64(2,i) == 0) && (bit64(4,i) == 0))
					{
						tmp1 = -3;
					}
					else if ((bit64(2,i) == 1) && (bit64(4,i) == 0))
					{
						tmp1 = -5;
					}
					else if ((bit64(2,i) == 1) && (bit64(4,i) == 1))
					{
						tmp1 = -7;
					}
					if ((bit64(3,i) == 0) && (bit64(5,i) == 1))
					{
						tmp2 = -1;
					}
					else if ((bit64(3,i) == 0) && (bit64(5,i) == 0))
					{
						tmp2 = -3;
					}
					else if ((bit64(3,i) == 1) && (bit64(5,i) == 0))
					{
						tmp2 = -5;
					}
					else if ((bit64(3,i) == 1) && (bit64(5,i) == 1))
					{
						tmp2 = -7;
					}
				}
			}
			complex<double> c1(tmp1, tmp2);
			bit_64qam_tmp(0, i) = c1;
			if (bit64(6, i) == 0)
			{
				if (bit64(7, i) == 0)
				{
					if ((bit64(8, i) == 0) && (bit64(10, i) == 1))
					{
						tmp1 = 1;
					}
					else if ((bit64(8, i) == 0) && (bit64(10, i) == 0))
					{
						tmp1 = 3;
					}
					else if ((bit64(8, i) == 1) && (bit64(10, i) == 0))
					{
						tmp1 = 5;
					}
					else if ((bit64(8, i) == 1) && (bit64(10, i) == 1))
					{
						tmp1 = 7;
					}
					if ((bit64(9, i) == 0) && (bit64(11, i) == 1))
					{
						tmp2 = 1;
					}
					else if ((bit64(9, i) == 0) && (bit64(11, i) == 0))
					{
						tmp2 = 3;
					}
					else if ((bit64(9, i) == 1) && (bit64(11, i) == 0))
					{
						tmp2 = 5;
					}
					else if ((bit64(9, i) == 1) && (bit64(11, i) == 1))
					{
						tmp2 = 7;
					}
				}
				else
				{
					if ((bit64(8, i) == 0) && (bit64(10, i) == 1))
					{
						tmp1 = 1;
					}
					else if ((bit64(8, i) == 0) && (bit64(10, i) == 0))
					{
						tmp1 = 3;
					}
					else if ((bit64(8, i) == 1) && (bit64(10, i) == 0))
					{
						tmp1 = 5;
					}
					else if ((bit64(8, i) == 1) && (bit64(10, i) == 1))
					{
						tmp1 = 7;
					}
					if ((bit64(9, i) == 0) && (bit64(11, i) == 1))
					{
						tmp2 = -1;
					}
					else if ((bit64(9, i) == 0) && (bit64(11, i) == 0))
					{
						tmp2 = -3;
					}
					else if ((bit64(9, i) == 1) && (bit64(11, i) == 0))
					{
						tmp2 = -5;
					}
					else if ((bit64(9, i) == 1) && (bit64(11, i) == 1))
					{
						tmp2 = -7;
					}
				}
			}
			else
			{
				if (bit64(7, i) == 0)
				{
					if ((bit64(8, i) == 0) && (bit64(10, i) == 1))
					{
						tmp1 = -1;
					}
					else if ((bit64(8, i) == 0) && (bit64(10, i) == 0))
					{
						tmp1 = -3;
					}
					else if ((bit64(8, i) == 1) && (bit64(10, i) == 0))
					{
						tmp1 = -5;
					}
					else if ((bit64(8, i) == 1) && (bit64(10, i) == 1))
					{
						tmp1 = -7;
					}
					if ((bit64(9, i) == 0) && (bit64(11, i) == 1))
					{
						tmp2 = 1;
					}
					else if ((bit64(9, i) == 0) && (bit64(11, i) == 0))
					{
						tmp2 = 3;
					}
					else if ((bit64(9, i) == 1) && (bit64(11, i) == 0))
					{
						tmp2 = 5;
					}
					else if ((bit64(9, i) == 1) && (bit64(11, i) == 1))
					{
						tmp2 = 7;
					}
				}
				else
				{
					if ((bit64(8, i) == 0) && (bit64(10, i) == 1))
					{
						tmp1 = -1;
					}
					else if ((bit64(8, i) == 0) && (bit64(10, i) == 0))
					{
						tmp1 = -3;
					}
					else if ((bit64(8, i) == 1) && (bit64(10, i) == 0))
					{
						tmp1 = -5;
					}
					else if ((bit64(8, i) == 1) && (bit64(10, i) == 1))
					{
						tmp1 = -7;
					}
					if ((bit64(9, i) == 0) && (bit64(11, i) == 1))
					{
						tmp2 = -1;
					}
					else if ((bit64(9, i) == 0) && (bit64(11, i) == 0))
					{
						tmp2 = -3;
					}
					else if ((bit64(9, i) == 1) && (bit64(11, i) == 0))
					{
						tmp2 = -5;
					}
					else if ((bit64(9, i) == 1) && (bit64(11, i) == 1))
					{
						tmp2 = -7;
					}
				}
			}
			complex<double> c2(tmp1, tmp2);
			bit_64qam_tmp(1, i) = c2;
		}
		bit_64qam = bit_64qam_tmp;
	}
	
	
}


Demapping::~Demapping()
{
}



MatrixXi Demapping::dp_pnc_bpsk(MatrixXd Rx_sig)
{
	MatrixXi res(Rx_sig.rows(), Rx_sig.cols());
	VectorXi s1xors2(Rx_sig.rows());
	VectorXi s1(Rx_sig.rows());
	VectorXi s2(Rx_sig.rows());
	//======================================
	// R1 Demapping
	//======================================
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		if (Rx_sig(i, 0) > -1.0&&Rx_sig(i,0)<=1.0)
		{
			s1xors2(i) = 1;
		}
		else
		{
			s1xors2(i) = 0;
		}
	}
	//======================================
	// R2 Demapping
	//======================================
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		if (Rx_sig(i, 1) >= 0.0)
		{
			s2(i) = 1;
		}
		else
		{
			s2(i) = 0;
		}
	}
	for (int i = 0; i < s1.size(); i++)
	{
		s1(i) = s1xors2(i) ^ s2(i);
	}
	
	res << s1, s2;
	//cout << s1<< endl;
	//cout << "--------" << endl;
	return res;
}


MatrixXi Demapping::dp_comp_bpsk(MatrixXd Rx_sig, MatrixXi A)
{
	MatrixXi res(Rx_sig.rows(), Rx_sig.cols());
	MatrixXi bit(2, 4);
	MatrixXi bit_bpsk(2, 4);
	bit << 0, 0, 1, 1,
		   0, 1, 0, 1;	  
	bit_bpsk << -1, -1, 1, 1,
	         	-1, 1, -1, 1;
	MatrixXi res_tmp=A*bit_bpsk;
	
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE*1.0;
		double a1 = Rx_sig(i, 0);
		double a2 = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(a1 - res_tmp(0, j), 2) + pow(a2 - res_tmp(1, j), 2);
			if (dis < min)
			{
				min = dis;
				res.row(i) = bit.col(j);
			}
		}
	}
	//cout << Rx_sig<<endl;
	//cout << "-------------" << endl;
	//cout << res<<endl;
	return res;
}

MatrixXd Demapping::dp_comp_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A)
{
	MatrixXd res(Rx_sig.rows()*2, 2);
	MatrixXd bit(4, 16);
	Matrix<complex<double>, Dynamic, Dynamic> bit_bpsk(2,16);
	bit << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		   0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
		   0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
		   0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
	
	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
	bit_bpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
		        c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;

	/*Matrix<complex<double>, Dynamic, Dynamic> A_temp(2, 2);
	A_temp << A(0, 0), A(0, 1),
		      A(1, 0), A(1, 1);*/
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_bpsk;
	//cout << res_tmp << endl;
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE*1.0;
		complex<double> a1 = Rx_sig(i, 0);
		complex<double> a2 = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(a1.real() - res_tmp(0, j).real(), 2) 
				+ pow(a1.imag() - res_tmp(0, j).imag(), 2)
				+ pow(a2.real() - res_tmp(1, j).real(), 2)
				+ pow(a2.imag() - res_tmp(1, j).imag(), 2);
			/*double dis = abs(a1.real() - res_tmp(0, j).real())
				+ abs(a1.imag() - res_tmp(0, j).imag())
				+ abs(a2.real() - res_tmp(1, j).real())
				+ abs(a2.imag() - res_tmp(1, j).imag());*/
			if (dis < min)
			{
				min = dis;
				res(2*i,0) = bit(0,j);
				res(2*i+1,0) = bit(1,j);
				res(2 * i, 1) = bit(2, j);
				res(2 * i + 1, 1) = bit(3, j);
			}
		}
	}
	//cout << Rx_sig<<endl;
	//cout << "-------------" << endl;
	//cout << res<<endl;
	return res;
}

MatrixXd Demapping::dp_comp_16qam(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A)
{
	MatrixXd res(Rx_sig.rows() * 4, 2);
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_16qam;
	//cout << res_tmp << endl;
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE*1.0;
		complex<double> a1 = Rx_sig(i, 0);
		complex<double> a2 = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(a1.real() - res_tmp(0, j).real(), 2)
				+ pow(a1.imag() - res_tmp(0, j).imag(), 2)
				+ pow(a2.real() - res_tmp(1, j).real(), 2)
				+ pow(a2.imag() - res_tmp(1, j).imag(), 2);
			/*double dis = abs(a1.real() - res_tmp(0, j).real())
			+ abs(a1.imag() - res_tmp(0, j).imag())
			+ abs(a2.real() - res_tmp(1, j).real())
			+ abs(a2.imag() - res_tmp(1, j).imag());*/
			if (dis < min)
			{
				min = dis;
				res(4 * i, 0) = bit_16(0, j);
				res(4 * i + 1, 0) = bit_16(1, j);
				res(4 * i + 2, 0) = bit_16(2, j);
				res(4 * i + 3, 0) = bit_16(3, j);
				res(4 * i, 1) = bit_16(4, j);
				res(4 * i + 1, 1) = bit_16(5, j);
				res(4 * i + 2, 1) = bit_16(6, j);
				res(4 * i + 3, 1) = bit_16(7, j);
			}
		}
	}
	//cout << Rx_sig<<endl;
	//cout << "-------------" << endl;
	//cout << res<<endl;
	return res;
}


MatrixXd Demapping::dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig)
{
	MatrixXd res(Rx_sig.rows() * 2, 2);
	VectorXd s1xors2(Rx_sig.rows()*2);
	VectorXd s1(Rx_sig.rows()*2);
	VectorXd s2(Rx_sig.rows()*2);
	//======================================
	// R1 Demapping
	//======================================
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		if (Rx_sig(i, 0).real() <= 1 && Rx_sig(i, 0).real() > -1)
		{
			s1xors2(2 * i) = 1;
		}
		else
		{
			s1xors2(2 * i) = 0;
		}
		if (Rx_sig(i, 0).imag() <= 1 && Rx_sig(i, 0).imag() > -1)
		{
			s1xors2(2 * i+1) = 1;
		}
		else
		{
			s1xors2(2 * i+1) = 0;
		}
	}
	//======================================
	// R2 Demapping
	//======================================
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		if (Rx_sig(i, 1).real() <= 0)
		{
			s2(2 * i) = 0;
		}
		else
		{
			s2(2 * i) = 1;
		}
		if (Rx_sig(i, 1).imag() <= 0)
		{
			s2(2 * i+1) = 0;
		}
		else
		{
			s2(2 * i+1) = 1;
		}

	}
	for (int i = 0; i < s1.size(); i++)
	{
		s1(i) = (int)s1xors2(i) ^(int)s2(i);
	}

	res << s1, s2;
	//cout << res<< endl;
	//cout << "--------" << endl;
	return res;

}

//MatrixXd Demapping::dp_llr_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A,double sigma,int bits)
//{
//	MatrixXd res(Rx_sig.rows() * 2, 2);
//	MatrixXd llr1(Rx_sig.rows() * 2, 2),llr2(Rx_sig.rows() * 2, 2),llr(Rx_sig.rows() * 2, 2);
//	Matrix<complex<double>, Dynamic, Dynamic> bit11_qpsk(2, 16), bit12_qpsk(2, 16), bit21_qpsk(2, 16), bit22_qpsk(2, 16);
//	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
//	bit11_qpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
//		          c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
//	bit12_qpsk << c1, c1, c1, c1, c3, c3, c3, c3, c2, c2, c2, c2, c4, c4, c4, c4,
//		          c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
//	bit21_qpsk << c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4,
//		          c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4;
//	bit22_qpsk << c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4,
//		          c1, c1, c1, c1, c3, c3, c3, c3, c2, c2, c2, c2, c4, c4, c4, c4;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp11 = A*bit11_qpsk;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp12 = A*bit12_qpsk;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp21 = A*bit21_qpsk;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp22 = A*bit22_qpsk;
//	for (int i = 0; i < Rx_sig.rows(); i++)
//	{
//		complex<double> a1 = Rx_sig(i, 0);
//		complex<double> a2 = Rx_sig(i, 1);
//		//==========
//		//LLR1
//		//==========
//		double sum1 = 0, sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2) / 2 * abs(pow(a1.real() - res_tmp11(0, j).real(), 2)+pow(a1.imag() - res_tmp11(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2) / 2 * abs(pow(a1.real() - res_tmp11(0, j).real(), 2) + pow(a1.imag() - res_tmp11(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr1(2 * i, 0) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp12(0, j).real(), 2)+pow(a1.imag() - res_tmp12(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp12(0, j).real(), 2) + pow(a1.imag() - res_tmp12(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr1(2 * i + 1, 0) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp21(0, j).real(), 2)+ pow(a1.imag() - res_tmp21(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp21(0, j).real(), 2) + pow(a1.imag() - res_tmp21(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr1(2 * i, 1) = log(sum1 / sum0);
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp22(0, j).real(), 2)+pow(a1.imag() - res_tmp22(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp22(0, j).real(), 2) + pow(a1.imag() - res_tmp22(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr1(2 * i + 1, 1) = log(sum1 / sum0);
//		//======
//		//LLR2
//		//======
//		sum1 = 0, sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp11(1, j).real(), 2) + pow(a2.imag() - res_tmp11(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp11(1, j).real(), 2) + pow(a2.imag() - res_tmp11(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr2(2 * i, 0) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp12(1, j).real(), 2) + pow(a2.imag() - res_tmp12(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp12(1, j).real(), 2) + pow(a2.imag() - res_tmp12(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr2(2 * i + 1, 0) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp21(1, j).real(), 2) + pow(a2.imag() - res_tmp21(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp21(1, j).real(), 2) + pow(a2.imag() - res_tmp21(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr2(2 * i, 1) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp22(1, j).real(), 2) + pow(a2.imag() - res_tmp22(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp22(1, j).real(), 2) + pow(a2.imag() - res_tmp22(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr2(2 * i + 1, 1) = log(sum1 / sum0);
//
//	}
//	//
//	//quantize llr1
//	//
//		MatrixXd llr1_q = llr1;
//		int qnum = pow(2, bits);
//		double maxvalue = 400.0, minvalue = -400.0;
//		double deta = (maxvalue - minvalue) / qnum;
//		for (int i = 0; i < llr1.rows(); i++)
//		{
//			llr1_q(i, 0) = llr1_q(i, 0) > maxvalue ? maxvalue - deta / 2 : llr1_q(i, 0);
//			llr1_q(i, 0) = llr1_q(i, 0) < minvalue ? minvalue + deta / 2 : llr1_q(i, 0);
//			llr1_q(i, 1) = llr1_q(i, 1) > maxvalue ? maxvalue - deta / 2 : llr1_q(i, 1);
//			llr1_q(i, 1) = llr1_q(i, 1) < minvalue ? minvalue + deta / 2 : llr1_q(i, 1);
//			for (int j = 1; j <= qnum; j++)
//			{
//				if (llr1_q(i, 0) >(minvalue + (j - 1)*deta) && llr1_q(i, 0) <= (minvalue + j*deta))
//				{
//					llr1_q(i, 0) = minvalue + j*deta - deta / 2;
//					break;
//				}
//			}
//			for (int j = 1; j <= qnum; j++)
//			{
//				if (llr1_q(i, 1) > (minvalue + (j - 1)*deta) && llr1_q(i, 1) <= (minvalue + j*deta))
//				{
//					llr1_q(i, 1) = minvalue + j*deta - deta / 2;
//					break;
//				}
//			}
//		}
//	//
//	//decode 
//	//
//	llr = llr1_q+llr2;
//	for (int i = 0; i < llr.rows(); i++)
//	{
//		res(i, 0) = llr(i, 0) >= 0 ? 1 : 0;
//		res(i, 1) = llr(i, 1) >= 0 ? 1 : 0;
//	}
//	
//	//cout << llr1_q << endl;
//	//cout << "---------"<< endl;
//	return res;
//}
MatrixXd Demapping::dp_llr_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, double sigma, int bits)
{
	MatrixXd res(Rx_sig.rows() * 2, 2);
	MatrixXd llr1(Rx_sig.rows() * 2, 2), llr2(Rx_sig.rows() * 2, 2), llr(Rx_sig.rows() * 2, 2);
	Matrix<complex<double>, Dynamic, Dynamic> bit11_qpsk(2, 16), bit12_qpsk(2, 16), bit21_qpsk(2, 16), bit22_qpsk(2, 16);
	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
	bit11_qpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
		c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
	bit12_qpsk << c1, c1, c1, c1, c3, c3, c3, c3, c2, c2, c2, c2, c4, c4, c4, c4,
		c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
	bit21_qpsk << c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4,
		c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4;
	bit22_qpsk << c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4,
		c1, c1, c1, c1, c3, c3, c3, c3, c2, c2, c2, c2, c4, c4, c4, c4;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp11 = A*bit11_qpsk;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp12 = A*bit12_qpsk;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp21 = A*bit21_qpsk;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp22 = A*bit22_qpsk;
	//double maxvalue = MIN_VALUE, minvalue = MAX_VALUE;
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		complex<double> a1 = Rx_sig(i, 0);
		complex<double> a2 = Rx_sig(i, 1);
		//==========
		//LLR1
		//==========
		double sum1 = 0, sum0 = 0,min=MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp11(0, j).real(), 2) + pow(a1.imag() - res_tmp11(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s= abs(pow(a1.real() - res_tmp11(0, j).real(), 2) + pow(a1.imag() - res_tmp11(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr1(2 * i, 0) = (sum0-sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp12(0, j).real(), 2) + pow(a1.imag() - res_tmp12(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp12(0, j).real(), 2) + pow(a1.imag() - res_tmp12(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr1(2 * i + 1, 0) = (sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp21(0, j).real(), 2) + pow(a1.imag() - res_tmp21(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp21(0, j).real(), 2) + pow(a1.imag() - res_tmp21(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr1(2 * i , 1) = (sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp22(0, j).real(), 2) + pow(a1.imag() - res_tmp22(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp22(0, j).real(), 2) + pow(a1.imag() - res_tmp22(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr1(2 * i + 1, 1) = (sum0 - sum1);

		//======
		//LLR2
		//======
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s =abs(pow(a2.real() - res_tmp11(1, j).real(), 2) + pow(a2.imag() - res_tmp11(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp11(1, j).real(), 2) + pow(a2.imag() - res_tmp11(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr2(2 * i, 0) = 1 / sigma*(sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp12(1, j).real(), 2) + pow(a2.imag() - res_tmp12(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp12(1, j).real(), 2) + pow(a2.imag() - res_tmp12(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr2(2 * i + 1, 0) = 1 / sigma*(sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp21(1, j).real(), 2) + pow(a2.imag() - res_tmp21(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp21(1, j).real(), 2) + pow(a2.imag() - res_tmp21(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr2(2 * i, 1) = 1 / sigma*(sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp22(1, j).real(), 2) + pow(a2.imag() - res_tmp22(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp22(1, j).real(), 2) + pow(a2.imag() - res_tmp22(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr2(2 * i + 1, 1) = 1 / sigma*(sum0 - sum1);

	}
	//
	//quantizer llr1
	//
	MatrixXd llr1_q = llr1;
	int qnum = pow(2, bits);
	double maxvalue = 5.0, minvalue = -5.0;
	double deta = (maxvalue - minvalue) / qnum;
	for (int i = 0; i < llr1.rows(); i++)
	{
		llr1_q(i, 0) = llr1_q(i, 0) > maxvalue ? maxvalue - deta / 2 : llr1_q(i, 0);
		llr1_q(i, 0) = llr1_q(i, 0) < minvalue ? minvalue + deta / 2 : llr1_q(i, 0);
		llr1_q(i, 1) = llr1_q(i, 1) > maxvalue ? maxvalue - deta / 2 : llr1_q(i, 1);
		llr1_q(i, 1) = llr1_q(i, 1) < minvalue ? minvalue + deta / 2 : llr1_q(i, 1);
		for (int j = 1; j <= qnum; j++)
		{
			if (llr1_q(i, 0) >(minvalue + (j - 1)*deta) && llr1_q(i, 0) <= (minvalue + j*deta))
			{
				llr1_q(i, 0) = minvalue + j*deta - deta / 2;
				break;
			}
		}
		for (int j = 1; j <= qnum; j++)
		{
			if (llr1_q(i, 1) > (minvalue + (j - 1)*deta) && llr1_q(i, 1) <= (minvalue + j*deta))
			{
				llr1_q(i, 1) = minvalue + j*deta - deta / 2;
				break;
			}
		}
	}






	//
	//decode
	//
	llr = llr1_q+llr2;
	for (int i = 0; i < llr.rows(); i++)
	{
		res(i, 0) = llr(i, 0) >= 0 ? 1 : 0;
		res(i, 1) = llr(i, 1) >= 0 ? 1 : 0;
	}

	return res;
}

MatrixXd Demapping::dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A)
{
	MatrixXd res(Rx_sig.rows() * 2, 2);
	VectorXd r1_s1(Rx_sig.rows() * 2);
	VectorXd r1_s2(Rx_sig.rows() * 2);
	VectorXd r2_s1(Rx_sig.rows() * 2);
	VectorXd r2_s2(Rx_sig.rows() * 2);
	MatrixXd bit(4, 16);
	Matrix<complex<double>, Dynamic, Dynamic> bit_qpsk(2, 16);
	MatrixXd ncm(9, 4);
	bit << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
		0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
		0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
	bit_qpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
		c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_qpsk;


	//
	//Computational network coding matrix
	//
	ncm << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1,
		1, 0, 1, 0,
		0, 1, 0, 1,
		1, 0, 0, 1,
		0, 1, 1, 0,
		1, 1, 0, 1;
	Matrix<double, 2, 4> nc[72];
	int n = 0;
	for (int i = 0; i < ncm.rows(); i++)
	{
		for (int j = 0; j < ncm.rows(); j++)
		{
			if (i != j)
			{
				nc[n].row(0) = ncm.row(i);
				nc[n].row(1) = ncm.row(j);
				n++;
			}
		}
	}
	//r1 matrix
	double maxdis = 0;
	int index1 = 0, index2 = 0;
	for (int k = 0; k < 72; k++)
	{
		double min = MAX_VALUE;
		MatrixXd pnc = nc[k] * bit;
		for (int i = 0; i < pnc.rows(); i++)
		{
			for (int j = 0; j < pnc.cols(); j++)
			{

				pnc(i, j) = (int)pnc(i, j) % 2;

			}
		}
		for (int i = 0; i < 16; i++)
		{
			for (int j = 0; j < 16; j++)
			{
				if (i != j)
				{
					double dic = pow(res_tmp(0, i).real() - res_tmp(0, j).real(), 2) + pow(res_tmp(0, i).imag() - res_tmp(0, j).imag(), 2);
					if (pnc.col(i) != pnc.col(j))
					{
						if (dic < min)
						{
							min = dic;
						}
					}
				}
			}
		}
		if (maxdis < min)
		{
			maxdis = min;
			index1 = k;
		}

	}
	//r2 matrix
	maxdis = 0;
	for (int k = 0; k < 72; k++)
	{
		if (k == index1)
		{
			continue;
		}
		double min = MAX_VALUE;
		MatrixXd pnc = nc[k] * bit;
		for (int i = 0; i < pnc.rows(); i++)
		{
			for (int j = 0; j < pnc.cols(); j++)
			{

				pnc(i, j) = (int)pnc(i, j) % 2;

			}
		}
		for (int i = 0; i < 16; i++)
		{
			for (int j = 0; j < 16; j++)
			{
				if (i != j)
				{
					double dic = pow(res_tmp(1, i).real() - res_tmp(1, j).real(), 2) + pow(res_tmp(1, i).imag() - res_tmp(1, j).imag(), 2);
					if (pnc.col(i) != pnc.col(j))
					{
						if (dic < min)
						{
							min = dic;
						}
					}
				}
			}
		}
		MatrixXd pnc_tmp(4, 4);
		pnc_tmp << nc[index1],
			nc[k];
		if (maxdis < min&&pnc_tmp.determinant()!=0)
		{
			maxdis = min;
			index2 = k;
		}

	}
	//
	//decode s1,s2
	//
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 0);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(0, j).real() - a.real(), 2) + pow(res_tmp(0, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				r1_s1(i * 2) = bit(0, j);
				r1_s1(i * 2 + 1) = bit(1, j);
				r1_s2(i * 2) = bit(2, j);
				r1_s2(i * 2 + 1) = bit(3, j);
			}
		}
	}
	//cout << r1_s1 << endl;
	//cout << "---------" << endl;
	//cout << r1_s2 << endl;
	//cout << "---------" << endl;
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(1, j).real() - a.real(), 2) + pow(res_tmp(1, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				r2_s1(i * 2) = bit(0, j);
				r2_s1(i * 2 + 1) = bit(1, j);
				r2_s2(i * 2) = bit(2, j);
				r2_s2(i * 2 + 1) = bit(3, j);
			}
		}
	}

	//
	//network coding
	//
	MatrixXd pnc_a(4, 4);
	pnc_a << nc[index1],
		nc[index2];
	//cout << pnc_a << endl;
	MatrixXd divpnc_a_tmp = pnc_a.inverse();
	MatrixXd divpnc_a = divpnc_a_tmp.array().abs();
	//cout << divpnc_a << endl;
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		VectorXd r1(4), r2(4);
		r1(0) = r1_s1(2 * i);
		r1(1) = r1_s1(2 * i + 1);
		r1(2) = r1_s2(2 * i);
		r1(3) = r1_s2(2 * i + 1);
		r2(0) = r2_s1(2 * i);
		r2(1) = r2_s1(2 * i + 1);
		r2(2) = r2_s2(2 * i);
		r2(3) = r2_s2(2 * i + 1);
		VectorXd tmp1(2), tmp2(2), tmp(4);
		tmp1 = nc[index1] * r1;
		tmp2 = nc[index2] * r2;
		tmp1(0) = (int)tmp1(0) % 2;
		tmp1(1) = (int)tmp1(1) % 2;
		tmp2(0) = (int)tmp2(0) % 2;
		tmp2(1) = (int)tmp2(1) % 2;
		tmp << tmp1,
			tmp2;
		VectorXd res_v = divpnc_a*tmp;
		res(2 * i, 0) = (int)res_v(0)%2;
		res(2 * i + 1, 0) = (int)res_v(1)%2;
		res(2 * i, 1) = (int)res_v(2)%2;
		res(2 * i + 1, 1) = (int)res_v(3)%2;
		/*res(2 * i, 0) = (int)res_v(0);
		res(2 * i + 1, 0) = (int)res_v(1);
		res(2 * i, 1) = (int)res_v(2);
		res(2 * i + 1, 1) = (int)res_v(3);*/
	}


	return res;
}

//MatrixXd Demapping::dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, int k)
//{
//	MatrixXd res(Rx_sig.rows() * 2, 2);
//	VectorXi pncr1(Rx_sig.rows()), pncr2(Rx_sig.rows());
//	Kmean km;
//	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
//	MatrixXd bit(4, 16);
//	Matrix<complex<double>, Dynamic, Dynamic> bit_qpsk(2, 16);
//	bit << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
//		0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
//		0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
//		0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
//	bit_qpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
//		c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_qpsk;
//	/*for (int i = 0; i < res_tmp.cols(); i++)
//	{
//		cout << res_tmp(0, i).real() << "   " << res_tmp(0, i).imag() << "			" << res_tmp(1, i).real() << "   " << res_tmp(1, i).imag() << endl;
//	}*/
//	MatrixXd data1(16, 2), data2(16, 2);
//	for (int i = 0; i < 16; i++)
//	{
//		data1(i, 0) = res_tmp(0, i).real();
//		data1(i, 1) = res_tmp(0, i).imag();
//		data2(i, 0) = res_tmp(1, i).real();
//		data2(i, 1) = res_tmp(1, i).imag();
//	}
//	VectorXi clures1(16), clures2(16);
//	clures1 = km.kmeans(data1, k);
//	clures2 = km.kmeans(data2, k);
//
//	for (int i = 0; i < Rx_sig.rows(); i++)
//	{
//		double min = MAX_VALUE;
//		complex<double> a = Rx_sig(i, 0);
//		for (int j = 0; j < res_tmp.cols(); j++)
//		{
//			double dis = pow(res_tmp(0, j).real() - a.real(), 2) + pow(res_tmp(0, j).imag() - a.imag(), 2);
//			if (dis < min)
//			{
//				min = dis;
//				pncr1(i) = clures1(j);
//			}
//		}
//	}
//	for (int i = 0; i < Rx_sig.rows(); i++)
//	{
//		double min = MAX_VALUE;
//		complex<double> a = Rx_sig(i, 1);
//		for (int j = 0; j < res_tmp.cols(); j++)
//		{
//			double dis = pow(res_tmp(1, j).real() - a.real(), 2) + pow(res_tmp(1, j).imag() - a.imag(), 2);
//			if (dis < min)
//			{
//				min = dis;
//				pncr2(i) = clures2(j);
//				res(i * 2,0) = bit(0, j);
//				res(i * 2 + 1,0) = bit(1, j);
//				res(i * 2,1) = bit(2, j);
//				res(i * 2 + 1,1) = bit(3, j);
//			}
//		}
//	}
//	for (int i = 0; i < pncr1.size(); i++)
//	{
//		int r11index = 0, r12index = 0, r21index = 0, r22index = 0;
//		for (int j = 0; j < clures1.size(); j++)
//		{
//			if (clures1(j) == pncr1(i))
//			{
//				r11index = r12index;
//				r12index = j;
//			}
//			if (clures2(j) == pncr2(i))
//			{
//				r21index = r22index;
//				r22index = j;
//			}
//		}
//		if (r11index==r21index&&r12index==r22index) 
//		{
//
//		}
//		else if (r11index == r21index)
//		{
//			res(2 * i, 0) = bit(0, r11index);
//			res(2 * i + 1, 0) = bit(1, r11index);
//			res(2 * i, 1) = bit(2, r11index);
//			res(2 * i + 1, 1) = bit(3, r11index);
//		}
//		else if (r11index == r22index)
//		{
//			res(2 * i, 0) = bit(0, r11index);
//			res(2 * i + 1, 0) = bit(1, r11index);
//			res(2 * i, 1) = bit(2, r11index);
//			res(2 * i + 1, 1) = bit(3, r11index);
//		}
//		else if (r12index == r21index)
//		{
//			res(2 * i, 0) = bit(0, r12index);
//			res(2 * i + 1, 0) = bit(1, r12index);
//			res(2 * i, 1) = bit(2, r12index);
//			res(2 * i + 1, 1) = bit(3, r12index);
//		}
//		else if (r12index == r22index)
//		{
//			res(2 * i, 0) = bit(0, r12index);
//			res(2 * i + 1, 0) = bit(1, r12index);
//			res(2 * i, 1) = bit(2, r12index);
//			res(2 * i + 1, 1) = bit(3, r12index);
//		}
//	}
//
//	//cout << clures1 << endl << clures2 << endl;
//	return res;
//}


MatrixXd Demapping::dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, int k)
{
	MatrixXd res(Rx_sig.rows() * 2, 2);
	VectorXi pncr1(Rx_sig.rows()), pncr2(Rx_sig.rows());
	Kmean km;
	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
	MatrixXd bit(4, 16);
	Matrix<complex<double>, Dynamic, Dynamic> bit_qpsk(2, 16);
	bit << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
		0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
		0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
	bit_qpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
		c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_qpsk;
	/*for (int i = 0; i < res_tmp.cols(); i++)
	{
	cout << res_tmp(0, i).real() << "   " << res_tmp(0, i).imag() << "			" << res_tmp(1, i).real() << "   " << res_tmp(1, i).imag() << endl;
	}*/
	MatrixXd data1(16, 2), data2(16, 2);
	for (int i = 0; i < 16; i++)
	{
		data1(i, 0) = res_tmp(0, i).real();
		data1(i, 1) = res_tmp(0, i).imag();
		data2(i, 0) = res_tmp(1, i).real();
		data2(i, 1) = res_tmp(1, i).imag();
	}
	VectorXi clures1(16), clures2(16);
	clures1 = km.kmeans(data1, k);
	clures2 = km.kmeans(data2, k);
	
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 0);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(0, j).real() - a.real(), 2) + pow(res_tmp(0, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				pncr1(i) = clures1(j);
				res(i * 2, 0) = bit(0, j);
				res(i * 2 + 1, 0) = bit(1, j);
				res(i *2, 1) = bit(2, j);
				res(i * 2 + 1, 1) = bit(3, j);
				
			}
		}
	}
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(1, j).real() - a.real(), 2) + pow(res_tmp(1, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				pncr2(i) = clures2(j);
				
			}
		}
	}
	
	
	for (int i = 0; i < pncr1.rows(); i++)
	{
		set<int> set1;
		set<int> set2;
		set<int> set1and2;
		for (int j = 0; j < clures1.size(); j++)
		{
			if (clures1(j) == pncr1(i))
			{
				set1.insert(j);
			}
			if (clures2(j) == pncr2(i))
			{
				set2.insert(j);
			}
		}
		set<int>::iterator it;
		for (it = set1.begin(); it != set1.end(); it++)
		{
			if (set2.find(*it) != set2.end())
			{
				set1and2.insert(*it);
			}
		}
		if (set1and2.size() == 0)
		{
			continue;
		}
		else if (set1and2.size() == 1)
		{
			int n = *set1and2.begin();
			res(i * 2, 0) = bit(0, n);
			res(i * 2 + 1, 0) = bit(1, n);
			res(i * 2, 1) = bit(2, n);
			res(i * 2 + 1, 1) = bit(3, n);
			

		}
		else
		{
			double min = MAX_VALUE;
			complex<double> a = Rx_sig(i, 0);
			for (it = set1.begin(); it != set1.end(); it++)
			{
				int n = *it;
				double dis = pow(res_tmp(0, n).real() - a.real(), 2) + pow(res_tmp(0, n).imag() - a.imag(), 2);
				if (dis < min)
				{
					min = dis;
					res(i * 2, 0) = bit(0, n);
					res(i * 2 + 1, 0) = bit(1, n);
					res(i * 2, 1) = bit(2, n);
					res(i * 2 + 1, 1) = bit(3, n);
				}

			}
		}
	}
	
	
	
	
	
	cout << data1 << endl;
	cout << "----------" << endl;
	cout << clures1 << endl;
	cout << "----------" << endl;
	cout << data2 << endl;
	cout << "----------" << endl;
	cout<< clures2 << endl;
	return res;
}

MatrixXd Demapping::dp_pnc_16qam(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, MatrixXd pnc_a)
{
	MatrixXd res(Rx_sig.rows() * 4, 2);
	VectorXd r1_s1(Rx_sig.rows() * 4);
	VectorXd r1_s2(Rx_sig.rows() * 4);
	VectorXd r2_s1(Rx_sig.rows() * 4);
	VectorXd r2_s2(Rx_sig.rows() * 4);
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_16qam;
	////
	////decode s1,s2
	////
	for (int i = 0; i < Rx_sig.rows(); ++i)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 0);
		for (int j = 0; j < res_tmp.cols(); ++j)
		{
			double dis = pow(res_tmp(0, j).real() - a.real(), 2) + pow(res_tmp(0, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				r1_s1(i * 4) = bit_16(0, j);
				r1_s1(i * 4 + 1) = bit_16(1, j);
				r1_s1(i * 4 + 2) = bit_16(2, j);
				r1_s1(i * 4 + 3) = bit_16(3, j);
				r1_s2(i * 4) = bit_16(4, j);
				r1_s2(i * 4 + 1) = bit_16(5, j);
				r1_s2(i * 4 + 2) = bit_16(6, j);
				r1_s2(i * 4 + 3) = bit_16(7, j);
			}
		}
	}
	////cout << r1_s1 << endl;
	////cout << "---------" << endl;
	////cout << r1_s2 << endl;
	////cout << "---------" << endl;
	for (int i = 0; i < Rx_sig.rows(); ++i)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); ++j)
		{
			double dis = pow(res_tmp(1, j).real() - a.real(), 2) + pow(res_tmp(1, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				r2_s1(i * 4) = bit_16(0, j);
				r2_s1(i * 4 + 1) = bit_16(1, j);
				r2_s1(i * 4 + 2) = bit_16(2, j);
				r2_s1(i * 4 + 3) = bit_16(3, j);
				r2_s2(i * 4) = bit_16(4, j);
				r2_s2(i * 4 + 1) = bit_16(5, j);
				r2_s2(i * 4 + 2) = bit_16(6, j);
				r2_s2(i * 4 + 3) = bit_16(7, j);
			}
		}
	}

	////
	////network coding
	////
	MatrixXd nc1(4, 8), nc2(4, 8);
	nc1 = pnc_a.topRows(4);
	nc2 = pnc_a.bottomRows(4);
	MatrixXd divpnc_a_tmp = pnc_a.inverse();
	MatrixXd divpnc_a = divpnc_a_tmp.array().abs();
	////cout << divpnc_a << endl;
	for (int i = 0; i < Rx_sig.rows(); ++i)
	{
		VectorXd r1(8), r2(8);
		r1(0) = r1_s1(4 * i);
		r1(1) = r1_s1(4 * i + 1);
		r1(2) = r1_s1(4 * i + 2);
		r1(3) = r1_s1(4 * i + 3);
		r1(4) = r1_s2(4 * i);
		r1(5) = r1_s2(4 * i + 1);
		r1(6) = r1_s2(4 * i + 2);
		r1(7) = r1_s2(4 * i + 3);
		r2(0) = r2_s1(4 * i);
		r2(1) = r2_s1(4 * i + 1);
		r2(2) = r2_s1(4 * i + 2);
		r2(3) = r2_s1(4 * i + 3);
		r2(4) = r2_s2(4 * i);
		r2(5) = r2_s2(4 * i + 1);
		r2(6) = r2_s2(4 * i + 2);
		r2(7) = r2_s2(4 * i + 3);
		VectorXd tmp1(4), tmp2(4), tmp(8);
		tmp1 = nc1 * r1;
		tmp2 = nc2 * r2;
		tmp1(0) = (int)tmp1(0) % 2;
		tmp1(1) = (int)tmp1(1) % 2;
		tmp1(2) = (int)tmp1(2) % 2;
		tmp1(3) = (int)tmp1(3) % 2;
		tmp2(0) = (int)tmp2(0) % 2;
		tmp2(1) = (int)tmp2(1) % 2;
		tmp2(2) = (int)tmp2(2) % 2;
		tmp2(3) = (int)tmp2(3) % 2;
		tmp << tmp1,
			tmp2;
		VectorXd res_v = divpnc_a*tmp;
		res(4 * i, 0) = (int)res_v(0) % 2;
		res(4 * i + 1, 0) = (int)res_v(1) % 2;
		res(4 * i + 2, 0) = (int)res_v(2) % 2;
		res(4 * i + 3, 0) = (int)res_v(3) % 2;
		res(4 * i, 1) = (int)res_v(4) % 2;
		res(4 * i + 1, 1) = (int)res_v(5) % 2;
		res(4 * i + 2, 1) = (int)res_v(6) % 2;
		res(4 * i + 3, 1) = (int)res_v(7) % 2;
		/*res(2 * i, 0) = (int)res_v(0);
		res(2 * i + 1, 0) = (int)res_v(1);
		res(2 * i, 1) = (int)res_v(2);
		res(2 * i + 1, 1) = (int)res_v(3);*/
	}


	return res;
}

MatrixXd Demapping::dp_pnc_16qam(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, VectorXi clures1, VectorXi clures2)
{
	MatrixXd res(Rx_sig.rows() * 4, 2);
	VectorXi pncr1(Rx_sig.rows()), pncr2(Rx_sig.rows());
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_16qam;
	
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 0);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(0, j).real() - a.real(), 2) + pow(res_tmp(0, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				pncr1(i) = clures1(j);
				res(i * 4, 0) = bit_16(0, j);
				res(i * 4 + 1, 0) = bit_16(1, j);
				res(i * 4 + 2, 0) = bit_16(2, j);
				res(i * 4 + 3, 0) = bit_16(3, j);
				res(i * 4, 1) = bit_16(4, j);
				res(i * 4 + 1, 1) = bit_16(5, j);
				res(i * 4 + 2, 1) = bit_16(6, j);
				res(i * 4 + 3, 1) = bit_16(7, j);
			}
		}
	}
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(1, j).real() - a.real(), 2) + pow(res_tmp(1, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				pncr2(i) = clures2(j);

			}
		}
	}


	for (int i = 0; i < pncr1.rows(); i++)
	{
		set<int> set1;
		set<int> set2;
		set<int> set1and2;
		for (int j = 0; j < clures1.size(); j++)
		{
			if (clures1(j) == pncr1(i))
			{
				set1.insert(j);
			}
			if (clures2(j) == pncr2(i))
			{
				set2.insert(j);
			}
		}
		set<int>::iterator it;
		for (it = set1.begin(); it != set1.end(); it++)
		{
			if (set2.find(*it) != set2.end())
			{
				set1and2.insert(*it);
			}
		}
		if (set1and2.size() == 0)
		{
			continue;
		}
		else if (set1and2.size() == 1)
		{
			int n = *set1and2.begin();
			res(i * 4, 0) = bit_16(0, n);
			res(i * 4 + 1, 0) = bit_16(1, n);
			res(i * 4 + 2, 0) = bit_16(2, n);
			res(i * 4 + 3, 0) = bit_16(3, n);
			res(i * 4, 1) = bit_16(4, n);
			res(i * 4 + 1, 1) = bit_16(5, n);
			res(i * 4 + 2, 1) = bit_16(6, n);
			res(i * 4 + 3, 1) = bit_16(7, n);


		}
		else
		{
			double min = MAX_VALUE;
			complex<double> a = Rx_sig(i, 0);
			for (it = set1.begin(); it != set1.end(); it++)
			{
				int n = *it;
				double dis = pow(res_tmp(0, n).real() - a.real(), 2) + pow(res_tmp(0, n).imag() - a.imag(), 2);
				if (dis < min)
				{
					min = dis;
					res(i * 4, 0) = bit_16(0, n);
					res(i * 4 + 1, 0) = bit_16(1, n);
					res(i * 4 + 2, 0) = bit_16(2, n);
					res(i * 4 + 3, 0) = bit_16(3, n);
					res(i * 4, 1) = bit_16(4, n);
					res(i * 4 + 1, 1) = bit_16(5, n);
					res(i * 4 + 2, 1) = bit_16(6, n);
					res(i * 4 + 3, 1) = bit_16(7, n);
				}

			}
		}
	}


	
	//cout << data1 << endl;
	//cout << "----------" << endl;
	//cout << clures1 << endl;
	//cout << "----------" << endl;
	//cout << data2 << endl;
	//cout << "----------" << endl;
	//cout << clures2 << endl;
	return res;
}

MatrixXd Demapping::dp_pnc_m16(Matrix<complex<double>, Dynamic, Dynamic> A)
{
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_16qam;
	MatrixXd nc1(4, 8), nc2(4, 8);
	double maxdis = 0;
	#pragma omp parallel for 
	for (int r1 = 0; r1 < ncm_16qam.rows(); ++r1)
	{
		MatrixXd nc1_tmp(4, 8);
		for (int r2 = 0; r2 < ncm_16qam.rows(); ++r2)
		{
			if (r2 == r1)
				continue;
			for (int r3 = 0; r3 < ncm_16qam.rows(); ++r3)
			{
				if (r3 == r2 || r3 == r1)
					continue;
				for (int r4 = 0; r4 < ncm_16qam.rows(); ++r4)
				{
					if (r4 == r1 || r4 == r2 || r4 == r3)
						continue;
					nc1_tmp.row(0) = ncm_16qam.row(r1);
					nc1_tmp.row(1) = ncm_16qam.row(r2);
					nc1_tmp.row(2) = ncm_16qam.row(r3);
					nc1_tmp.row(3) = ncm_16qam.row(r4);
					if (1)
					{
						double min = MAX_VALUE;
						MatrixXd pnc = nc1_tmp * bit_16;
						for (int i = 0; i < pnc.rows(); ++i)
						{
							for (int j = 0; j < pnc.cols(); ++j)
							{

								pnc(i, j) = (int)pnc(i, j) % 2;

							}
						}
						for (int i = 0; i < res_tmp.cols(); ++i)
						{
							for (int j = 0; j < res_tmp.cols(); ++j)
							{
								if (i != j)
								{
									double dic = pow(res_tmp(0, i).real() - res_tmp(0, j).real(), 2) + pow(res_tmp(0, i).imag() - res_tmp(0, j).imag(), 2);
									if (pnc.col(i) != pnc.col(j))
									{
										if (dic < min)
										{
											min = dic;
										}
									}
								}
							}
						}
						if (maxdis < min)
						{
							maxdis = min;
							nc1 = nc1_tmp;
						}
					}

				}
			}
		}
	}
	////r2 matrix
	maxdis = 0;
	#pragma omp parallel for
	for (int r1 = 0; r1 < ncm_16qam.rows(); ++r1)
	{
		Matrix<double, 4, 8> nc2_tmp;
		for (int r2 = 0; r2 < ncm_16qam.rows(); ++r2)
		{
			if (r2 == r1)
				continue;
			for (int r3 = 0; r3 < ncm_16qam.rows(); ++r3)
			{
				if (r3 == r2 || r3 == r1)
					continue;
				for (int r4 = 0; r4 < ncm_16qam.rows(); ++r4)
				{
					if (r4 == r1 || r4 == r2 || r4 == r3)
						continue;
					nc2_tmp.row(0) = ncm_16qam.row(r1);
					nc2_tmp.row(1) = ncm_16qam.row(r2);
					nc2_tmp.row(2) = ncm_16qam.row(r3);
					nc2_tmp.row(3) = ncm_16qam.row(r4);
					if (nc1.determinant() != 0)
					{
						double min = MAX_VALUE;
						MatrixXd pnc = nc2_tmp * bit_16;
						for (int i = 0; i < pnc.rows(); ++i)
						{
							for (int j = 0; j < pnc.cols(); ++j)
							{

								pnc(i, j) = (int)pnc(i, j) % 2;

							}
						}
						for (int i = 0; i < res_tmp.cols(); ++i)
						{
							for (int j = 0; j < res_tmp.cols(); ++j)
							{
								if (i != j)
								{
									double dic = pow(res_tmp(0, i).real() - res_tmp(0, j).real(), 2) + pow(res_tmp(0, i).imag() - res_tmp(0, j).imag(), 2);
									if (pnc.col(i) != pnc.col(j))
									{
										if (dic < min)
										{
											min = dic;
										}
									}
								}
							}
						}
						MatrixXd pnc_tmp(8, 8);
						pnc_tmp << nc1,
							nc2_tmp;
						if (maxdis < min&&pnc_tmp.determinant() != 0)
						{
							maxdis = min;
							nc2 = nc2_tmp;
						}
					}

				}
			}
		}
	}
	MatrixXd pnc_a(8, 8);
	pnc_a << nc1,
		nc2;
	return pnc_a;
}

VectorXi Demapping::dp_creat_cluster_16qam(Matrix<complex<double>, Dynamic, Dynamic> A, int k,int n)
{
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_16qam;
	Kmean km;
	MatrixXd data1(256, 2), data2(256, 2);	
	for (int i = 0; i < 256; i++)
	{
		data1(i, 0) = res_tmp(0, i).real();
		data1(i, 1) = res_tmp(0, i).imag();
		data2(i, 0) = res_tmp(1, i).real();
		data2(i, 1) = res_tmp(1, i).imag();
	}
	VectorXi clures, clures2;
	if (n == 1)
	{
		clures = km.kmean(data1, k);
	}
	else
	{
		clures = km.kmean(data2, k);
	}
	return clures;
}

MatrixXd Demapping::dp_pnc_64qam(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, VectorXi clures1, VectorXi clures2)
{
	MatrixXd res(Rx_sig.rows() * 6, 2);
	VectorXi pncr1(Rx_sig.rows()), pncr2(Rx_sig.rows());
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_64qam;

	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 0);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(0, j).real() - a.real(), 2) + pow(res_tmp(0, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				pncr1(i) = clures1(j);
				res(i * 6, 0) = bit_64(0, j);
				res(i * 6 + 1, 0) = bit_64(1, j);
				res(i * 6 + 2, 0) = bit_64(2, j);
				res(i * 6 + 3, 0) = bit_64(3, j);
				res(i * 6 + 4, 0) = bit_64(4, j);
				res(i * 6 + 5, 0) = bit_64(5, j);
				res(i * 6, 1) = bit_64(6, j);
				res(i * 6 + 1, 1) = bit_64(7, j);
				res(i * 6 + 2, 1) = bit_64(8, j);
				res(i * 6 + 3, 1) = bit_64(9, j);
				res(i * 6 + 4, 1) = bit_64(10, j);
				res(i * 6 + 5, 1) = bit_64(11, j);
			}
		}
	}
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(1, j).real() - a.real(), 2) + pow(res_tmp(1, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				pncr2(i) = clures2(j);

			}
		}
	}


	for (int i = 0; i < pncr1.rows(); i++)
	{
		set<int> set1;
		set<int> set2;
		set<int> set1and2;
		for (int j = 0; j < clures1.size(); j++)
		{
			if (clures1(j) == pncr1(i))
			{
				set1.insert(j);
			}
			if (clures2(j) == pncr2(i))
			{
				set2.insert(j);
			}
		}
		set<int>::iterator it;
		for (it = set1.begin(); it != set1.end(); it++)
		{
			if (set2.find(*it) != set2.end())
			{
				set1and2.insert(*it);
			}
		}
		if (set1and2.size() == 0)
		{
			continue;
		}
		else if (set1and2.size() == 1)
		{
			int n = *set1and2.begin();
			res(i * 6, 0) = bit_64(0, n);
			res(i * 6 + 1, 0) = bit_64(1, n);
			res(i * 6 + 2, 0) = bit_64(2, n);
			res(i * 6 + 3, 0) = bit_64(3, n);
			res(i * 6 + 4, 0) = bit_64(4, n);
			res(i * 6 + 5, 0) = bit_64(5, n);
			res(i * 6, 1) = bit_64(6, n);
			res(i * 6 + 1, 1) = bit_64(7, n);
			res(i * 6 + 2, 1) = bit_64(8, n);
			res(i * 6 + 3, 1) = bit_64(9, n);
			res(i * 6 + 4, 1) = bit_64(10, n);
			res(i * 6 + 5, 1) = bit_64(11, n);


		}
		else
		{
			double min = MAX_VALUE;
			complex<double> a = Rx_sig(i, 0);
			for (it = set1.begin(); it != set1.end(); it++)
			{
				int n = *it;
				double dis = pow(res_tmp(0, n).real() - a.real(), 2) + pow(res_tmp(0, n).imag() - a.imag(), 2);
				if (dis < min)
				{
					min = dis;
					res(i * 6, 0) = bit_64(0, n);
					res(i * 6 + 1, 0) = bit_64(1, n);
					res(i * 6 + 2, 0) = bit_64(2, n);
					res(i * 6 + 3, 0) = bit_64(3, n);
					res(i * 6 + 4, 0) = bit_64(4, n);
					res(i * 6 + 5, 0) = bit_64(5, n);
					res(i * 6, 1) = bit_64(6, n);
					res(i * 6 + 1, 1) = bit_64(7, n);
					res(i * 6 + 2, 1) = bit_64(8, n);
					res(i * 6 + 3, 1) = bit_64(9, n);
					res(i * 6 + 4, 1) = bit_64(10, n);
					res(i * 6 + 5, 1) = bit_64(11, n);
				}

			}
		}
	}



	//cout << data1 << endl;
	//cout << "----------" << endl;
	//cout << clures1 << endl;
	//cout << "----------" << endl;
	//cout << data2 << endl;
	//cout << "----------" << endl;
	//cout << clures2 << endl;
	return res;
}

VectorXi Demapping::dp_creat_cluster_64qam(Matrix<complex<double>, Dynamic, Dynamic> A, int k, int n)
{
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_64qam;
	Kmean km;
	ofstream log1, log2;
	log1.open("data1.txt", ios::trunc | ios::out);
	log2.open("data2.txt", ios::trunc | ios::out);
	MatrixXd data1(4096, 2), data2(4096, 2);
	for (int i = 0; i < 4096; i++)
	{
		data1(i, 0) = res_tmp(0, i).real();
		data1(i, 1) = res_tmp(0, i).imag();
		log1 << data1(i, 0) << "  " << data1(i, 1) << endl;
		data2(i, 0) = res_tmp(1, i).real();
		data2(i, 1) = res_tmp(1, i).imag();
		log2 << data2(i, 0) << "  " << data2(i, 1) << endl;
	}
	VectorXi clures, clures2;
	if (n == 1)
	{
		clures = km.kmean(data1, k);
	}
	else
	{
		clures = km.kmean(data2, k);
	}
	log1.close();
	log2.close();
	return clures;
}
