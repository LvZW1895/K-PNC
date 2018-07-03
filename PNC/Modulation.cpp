#include "stdafx.h"
#include "Modulation.h"


Modulation::Modulation()
{
}




Modulation::~Modulation()
{
}


VectorXi Modulation::mod_BPSK(VectorXi msg_bits)
{
	int len = msg_bits.size();
	VectorXi res(len);
	for (int i = 0; i < len; i++)
	{
		if (msg_bits(i) == 0)
			res(i) = -1;
		else if (msg_bits(i) == 1)
			res(i) = 1;
	}
	return res;
}


Matrix<complex<double>, Dynamic, Dynamic> Modulation::mod_QPSK(VectorXi msg_bits)
{
	int len = msg_bits.size()/2;
	Matrix<complex<double>, Dynamic, Dynamic> res(len,1);
	for (int i = 0; i < len; i++)
	{
		int tmp1 = msg_bits(2*i) == 0 ? -1.0: 1.0;
		int tmp2 = msg_bits(2 * i + 1) == 0 ? -1.0 : 1.0;
		complex<int> c(tmp1, tmp2);
		res(i, 0) = c;

	}
	return res;
}

Matrix<complex<double>, Dynamic, Dynamic> Modulation::mod_16QAM(VectorXi msg_bits)
{
	int len = msg_bits.size() / 4;
	Matrix<complex<double>, Dynamic, Dynamic> res(len, 1);
	for (int i = 0; i < len; i++)
	{
		int tmp1, tmp2;
		if (msg_bits(4 * i) == 0)
		{
			if (msg_bits(4 * i+2) == 0)
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
			if (msg_bits(4 * i + 2) == 0)
			{
				tmp1 = -1;
			}
			else
			{
				tmp1 = -3;
			}
		}
		if (msg_bits(4 * i + 1) == 0)
		{
			if (msg_bits(4 * i + 3) == 0)
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
			if (msg_bits(4 * i + 3) == 0)
			{
				tmp2 = -1;
			}
			else
			{
				tmp2 = -3;
			}
		}
		complex<double> c(tmp1, tmp2);
		res(i, 0) = c;
	}
	return res;
}

Matrix<complex<double>, Dynamic, Dynamic> Modulation::mod_64QAM(VectorXi msg_bits)
{
	int len = msg_bits.size() / 6;
	Matrix<complex<double>, Dynamic, Dynamic> res(len, 1);
	for (int i = 0; i < len; i++)
	{
		int tmp1, tmp2;
		if (msg_bits(6 * i) == 0)
		{
			if (msg_bits(6 * i + 1) == 0)
			{
				if ((msg_bits(6 * i + 2) == 0) && (msg_bits(6 * i + 4) == 1))
				{
					tmp1 = 1;
				}
				else if ((msg_bits(6 * i + 2) == 0) && (msg_bits(6 * i + 4) == 0))
				{
					tmp1 = 3;
				}
				else if ((msg_bits(6 * i + 2) == 1) && (msg_bits(6 * i + 4) == 0))
				{
					tmp1 = 5;
				}
				else if ((msg_bits(6 * i + 2) == 1) && (msg_bits(6 * i + 4) == 1))
				{
					tmp1 = 7;
				}
				if ((msg_bits(6 * i + 3) == 0) && (msg_bits(6 * i + 5) == 1))
				{
					tmp2 = 1;
				}
				else if ((msg_bits(6 * i + 3) == 0) && (msg_bits(6 * i + 5) == 0))
				{
					tmp2 = 3;
				}
				else if ((msg_bits(6 * i + 3) == 1) && (msg_bits(6 * i + 5) == 0))
				{
					tmp2 = 5;
				}
				else if ((msg_bits(6 * i + 3) == 1) && (msg_bits(6 * i + 5) == 1))
				{
					tmp2 = 7;
				}
			}
			else
			{
				if ((msg_bits(6 * i + 2) == 0) && (msg_bits(6 * i + 4) == 1))
				{
					tmp1 = 1;
				}
				else if ((msg_bits(6 * i + 2) == 0) && (msg_bits(6 * i + 4) == 0))
				{
					tmp1 = 3;
				}
				else if ((msg_bits(6 * i + 2) == 1) && (msg_bits(6 * i + 4) == 0))
				{
					tmp1 = 5;
				}
				else if ((msg_bits(6 * i + 2) == 1) && (msg_bits(6 * i + 4) == 1))
				{
					tmp1 = 7;
				}
				if ((msg_bits(6 * i + 3) == 0) && (msg_bits(6 * i + 5) == 1))
				{
					tmp2 = -1;
				}
				else if ((msg_bits(6 * i + 3) == 0) && (msg_bits(6 * i + 5) == 0))
				{
					tmp2 = -3;
				}
				else if ((msg_bits(6 * i + 3) == 1) && (msg_bits(6 * i + 5) == 0))
				{
					tmp2 = -5;
				}
				else if ((msg_bits(6 * i + 3) == 1) && (msg_bits(6 * i + 5) == 1))
				{
					tmp2 = -7;
				}
			}
		}
		else
		{
			if (msg_bits(6 * i + 1) == 0)
			{
				if ((msg_bits(6 * i + 2) == 0) && (msg_bits(6 * i + 4) == 1))
				{
					tmp1 = -1;
				}
				else if ((msg_bits(6 * i + 2) == 0) && (msg_bits(6 * i + 4) == 0))
				{
					tmp1 = -3;
				}
				else if ((msg_bits(6 * i + 2) == 1) && (msg_bits(6 * i + 4) == 0))
				{
					tmp1 = -5;
				}
				else if ((msg_bits(6 * i + 2) == 1) && (msg_bits(6 * i + 4) == 1))
				{
					tmp1 = -7;
				}
				if ((msg_bits(6 * i + 3) == 0) && (msg_bits(6 * i + 5) == 1))
				{
					tmp2 = 1;
				}
				else if ((msg_bits(6 * i + 3) == 0) && (msg_bits(6 * i + 5) == 0))
				{
					tmp2 = 3;
				}
				else if ((msg_bits(6 * i + 3) == 1) && (msg_bits(6 * i + 5) == 0))
				{
					tmp2 = 5;
				}
				else if ((msg_bits(6 * i + 3) == 1) && (msg_bits(6 * i + 5) == 1))
				{
					tmp2 = 7;
				}
			}
			else
			{
				if ((msg_bits(6 * i + 2) == 0) && (msg_bits(6 * i + 4) == 1))
				{
					tmp1 = -1;
				}
				else if ((msg_bits(6 * i + 2) == 0) && (msg_bits(6 * i + 4) == 0))
				{
					tmp1 = -3;
				}
				else if ((msg_bits(6 * i + 2) == 1) && (msg_bits(6 * i + 4) == 0))
				{
					tmp1 = -5;
				}
				else if ((msg_bits(6 * i + 2) == 1) && (msg_bits(6 * i + 4) == 1))
				{
					tmp1 = -7;
				}
				if ((msg_bits(6 * i + 3) == 0) && (msg_bits(6 * i + 5) == 1))
				{
					tmp2 = -1;
				}
				else if ((msg_bits(6 * i + 3) == 0) && (msg_bits(6 * i + 5) == 0))
				{
					tmp2 = -3;
				}
				else if ((msg_bits(6 * i + 3) == 1) && (msg_bits(6 * i + 5) == 0))
				{
					tmp2 = -5;
				}
				else if ((msg_bits(6 * i + 3) == 1) && (msg_bits(6 * i + 5) == 1))
				{
					tmp2 = -7;
				}
			}
		}
		complex<double> c(tmp1, tmp2);
		res(i, 0) = c;
	}
	return res;
}
