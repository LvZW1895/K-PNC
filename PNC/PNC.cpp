// PNC.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Modulation.h"
#include "Channel.h"
#include "Demapping.h"
#include "Quantizer.h"

int main()
{
	//srand((unsigned int)time(NULL));
	srand(20180529);
	ofstream log;
	int type =PNC;
	int modtype = QAM64;
	/*switch (type)
	{
	case CoMP:
		log.open("comp_result_throughput.txt", ios::trunc | ios::out);
		break;
	case PNC:
		log.open("pnc_result_throughput.txt", ios::trunc | ios::out);
		break;
	case nidealCoMP:
		log.open("nicomp_result_throughput.txt", ios::trunc | ios::out);
		break;
	default:
		break;
	}*/
	
	double relaytime=4.0;
	double frametime = 1.0;
	int block_num = 1;
	int msg_len = 6000;
	double Es = 1;
	int total_bit = 0;
	int error = 0;
	int recivebit = 0;
	int qtbits = 2;
	Modulation mod;
	Channel channel;
	Quantizer qt;
	Demapping dp(modtype);
	VectorXi EsN0dB=VectorXi::LinSpaced(21,0,20);
	VectorXd EsN0(21);
	for (int i = 0; i < 21; i++)
	{
		EsN0(i) = pow(10.0,0.1*EsN0dB(i));
	}
	VectorXd N0 = Es*EsN0.array().pow(-1);
	VectorXd sqrt_N0 = N0.array().sqrt();
	VectorXd sigma2 = N0 / 2;
	VectorXd sigma = sigma2.array().sqrt();
	VectorXi msg_u1(msg_len), msg_u2(msg_len);
	Matrix<complex<double>, Dynamic, Dynamic> A[100];
	MatrixXd pnc_a[100];
	VectorXi clures1[100], clures2[100];
	//A = channel.CreatH();
	/*A << 1, 1,
		0, 1;*/
    //cout << A << endl;
	
	
	
	ostringstream oss;
	oss << "nidealcomp_result" << qtbits << "_ber.txt";
	ostringstream oss1;
	oss1 << "llrcomp_result" << qtbits << "_ber.txt";
	switch (type)
	{
	case CoMP:
		log.open("comp_result_ber.txt", ios::trunc | ios::out);
		break;
	case nidealCoMP:
		log.open(oss.str(), ios::trunc | ios::out);
		break;
	case llrCoMP:
		log.open(oss1.str(), ios::trunc | ios::out);
		break;
	case PNC:
		log.open("pnc_result_ber.txt", ios::trunc | ios::out);
		break;
	default:
		break;
	}
	//cout << "#EsN0dB       #err         THROUGHPUT" << endl;
	//cout << sigma << endl;
	//log << "#EsN0dB     SER" << endl;
	switch (type)
	{
	case CoMP:
		relaytime = 0;
		break;
	case PNC:
		relaytime = 2.0;
		break;
	case nidealCoMP:
		relaytime = 4.0;
		break;
	default:
		break;
	}
	
	for (int k = 0; k <block_num; ++k)
	{
		A[k]= channel.CreatH();
		//pnc_a[k] = dp.dp_pnc_m16(A[k]);
		if (modtype == QAM16)
		{
			clures1[k] = dp.dp_creat_cluster_16qam(A[k], 16, 1);
			clures2[k] = dp.dp_creat_cluster_16qam(A[k], 16, 2);
		}
		else if (modtype == QAM64)
		{
			clures1[k] = dp.dp_creat_cluster_64qam(A[k], 64, 1);
			clures2[k] = dp.dp_creat_cluster_64qam(A[k], 64, 2);
		}
		
		//cout << A[k] << endl;
	}
	ofstream log3, log4;
	log3.open("clures1.txt", ios::trunc | ios::out);
	log4.open("clures2.txt", ios::trunc | ios::out);
	for (int i = 0; i < clures1[0].size(); i++)
	{
		log3 << clures1[0](i) << endl;
		log4 << clures2[0](i) << endl;
	}

	log3.close();
	log4.close();

	cout << "#EsN0dB       #err         SER" << endl;
	for (int k = 0; k <EsN0dB.size(); ++k)
	{
		error = 0;
		total_bit = 0;
		recivebit = 0;
		for (int i = 0; i < block_num; ++i)
		{
			//A = channel.CreatH();
			//======================================
			// generate message bits
			//======================================
			
			VectorXd msg_u1_tmp = VectorXd::Random(msg_len) + VectorXd::Constant(msg_len, 1.0);
			VectorXd msg_u2_tmp = VectorXd::Random(msg_len) + VectorXd::Constant(msg_len, 1.0);
			for (int j = 0; j< msg_len; ++j)
			{
				msg_u1(j) = (int)msg_u1_tmp(j);
				msg_u2(j) = (int)msg_u2_tmp(j);
			}
			//======================================
			// modulation
			//======================================
			Matrix<complex<double>, Dynamic, Dynamic> cv_txsig_u1, cv_txsig_u2;
			switch (modtype)
			{
			case BPSK:
				// cv_txsug_u1 = mod.mod_BPSK(msg_u1);
				 //cv_txsug_u2 = mod.mod_BPSK(msg_u2);
				break;
			case QPSK:
				cv_txsig_u1 = mod.mod_QPSK(msg_u1);
				cv_txsig_u2 = mod.mod_QPSK(msg_u2);
				break;
			case QAM16:
				cv_txsig_u1 = mod.mod_16QAM(msg_u1);
				cv_txsig_u2 = mod.mod_16QAM(msg_u2); 
				break;
			case QAM64:
				cv_txsig_u1 = mod.mod_64QAM(msg_u1);
				cv_txsig_u2 = mod.mod_64QAM(msg_u2);
				break;
			default:
				break;
			}
			/*cout << msg_u1 << endl;
			cout << "----------" << endl;
			cout << cv_txsig_u1 << endl;*/
			//======================================
			// Channel
			//======================================
			Matrix<complex<double>, Dynamic, Dynamic> Rx_sig = channel.AWGN(sigma(k), cv_txsig_u1, cv_txsig_u2, A[i]);
			//======================================
			// PNC Demapping
			//======================================
			//cout << Rx_sig << endl;
			MatrixXd res;
			Matrix<complex<double>, Dynamic, Dynamic> qt_rx;
			if (type == CoMP)
			{
				switch (modtype)
				{
				case BPSK:
					//res = dp.dp_comp_bpsk(Rx_sig, A);
					break;
				case QPSK:
					res = dp.dp_comp_qpsk(Rx_sig, A[i]);
					break;
				case QAM16:
					res = dp.dp_comp_16qam(Rx_sig, A[i]);
					break;
				default:
					break;
				}
				
			}
			else if (type == nidealCoMP)
			{
				switch (modtype)
				{
				case BPSK:
					//res = dp.dp_pnc_bpsk(Rx_sig);
					break;
				case QPSK:
					qt_rx= qt.quantize(Rx_sig,qtbits,2.0*sqrt(2),-2.0*sqrt(2));
					//A = qt.quantizeH(A,4);
					//cout << A << endl;
					res = dp.dp_comp_qpsk(qt_rx, A[i]);
					break;
				default:
					break;
				}
			}
			else if (type == llrCoMP)
			{
				switch (modtype)
				{
				case BPSK:
					//res = dp.dp_pnc_bpsk(Rx_sig);
					break;
				case QPSK:	
					res = dp.dp_llr_qpsk(Rx_sig, A[i], sigma(k),qtbits);
					break;
				default:
					break;
				}
			}
			else if (type == PNC)
			{
				switch (modtype)
				{
				case BPSK:
					//res = dp.dp_pnc_bpsk(Rx_sig);
					break;
				case QPSK:
					res = dp.dp_pnc_qpsk(Rx_sig,A[i],4);
					break;
				case QAM16:
					res= dp.dp_pnc_16qam(Rx_sig, A[i],clures1[i],clures2[i]);
					break;
				case QAM64:
					res = dp.dp_pnc_64qam(Rx_sig, A[i], clures1[i], clures2[i]);
					break;
				default:
					break;
				}
			}
			/*cout << msg_u1 << endl;
			cout << "--------------" << endl;
			cout << msg_u2 << endl;
			cout << "--------------" << endl;
			cout << res << endl;*/
			
			VectorXd res_u1 = res.col(0);
			VectorXd res_u2 = res.col(1);
			
			//======================================
			// BER calculation
			//======================================
			for (int j = 0; j < res_u1.size(); j++)
			{
				if (res_u1(j) != msg_u1(j))
				{
					error++;
				}		
				if (res_u2(j) != msg_u2(j))
				{
					error++;
				}
					
			}
			total_bit += 2 * msg_len;
			//recivebit = total_bit - error;
			//cout << "--------------" << endl;
			//cout << error << endl;
		}		
		
		std::printf(" %3i        %6i      %1.3e\n", EsN0dB(k), error, ((double)error) / total_bit);
	    log << EsN0dB(k)<<"        " << ((double)error) / total_bit << endl;
		//std::printf(" %3i        %6i      %1.3e\n", EsN0dB(k), error, ((double)recivebit) / (frametime+relaytime)/ block_num);
		//log << EsN0dB(k) << "        " << ((double)recivebit) / (frametime + relaytime)/ block_num << endl;
	}
	log.close();
	
    return 0;
}

