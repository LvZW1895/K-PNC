#pragma once
class Demapping
{
public:
	MatrixXd bit_16,bit_64;
	Matrix<complex<double>, Dynamic, Dynamic> bit_16qam, bit_64qam;
	MatrixXd ncm_16qam;
public:
	Demapping(int k);
	~Demapping();
	MatrixXi dp_pnc_bpsk(MatrixXd Rx_sig);
	MatrixXi dp_comp_bpsk(MatrixXd Rx_sig, MatrixXi A);
	MatrixXd dp_comp_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A);
	MatrixXd dp_comp_16qam(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A);

	MatrixXd dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig);
	MatrixXd dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A);
	MatrixXd dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A,int k);

	MatrixXd dp_pnc_16qam(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, MatrixXd pnc_a);
	MatrixXd dp_pnc_16qam(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, VectorXi clures1, VectorXi clures2);
	MatrixXd dp_pnc_m16(Matrix<complex<double>, Dynamic, Dynamic> A);
	VectorXi dp_creat_cluster_16qam(Matrix<complex<double>, Dynamic, Dynamic> A, int k,int n);

	MatrixXd dp_pnc_64qam(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, VectorXi clures1, VectorXi clures2);
	VectorXi dp_creat_cluster_64qam(Matrix<complex<double>, Dynamic, Dynamic> A, int k, int n);


	MatrixXd dp_llr_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A,double sigma,int bits);
};

