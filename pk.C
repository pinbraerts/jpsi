#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TH1.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <array>

using Kinematics = std::array<float, 3 * 4>;

void enable_kinematics(TTree& tree) {
	tree.SetBranchStatus("*", 0);
	tree.SetBranchStatus("sB_mu1_px", 1);
	tree.SetBranchStatus("sB_mu1_py", 1);
	tree.SetBranchStatus("sB_mu1_pz", 1);
	tree.SetBranchStatus("sB_mu2_px", 1);
	tree.SetBranchStatus("sB_mu2_py", 1);
	tree.SetBranchStatus("sB_mu2_pz", 1);
	tree.SetBranchStatus("sB_trk1_px", 1);
	tree.SetBranchStatus("sB_trk1_py", 1);
	tree.SetBranchStatus("sB_trk1_pz", 1);
	tree.SetBranchStatus("sB_trk2_px", 1);
	tree.SetBranchStatus("sB_trk2_py", 1);
	tree.SetBranchStatus("sB_trk2_pz", 1);
}

void setup_kinematics(TTree& tree, Kinematics& kinematics) {
	tree.SetBranchAddress( "sB_mu1_px", &kinematics[0]);
	tree.SetBranchAddress( "sB_mu1_py", &kinematics[1]);
	tree.SetBranchAddress( "sB_mu1_pz", &kinematics[2]);
	tree.SetBranchAddress( "sB_mu2_px", &kinematics[3]);
	tree.SetBranchAddress( "sB_mu2_py", &kinematics[4]);
	tree.SetBranchAddress( "sB_mu2_pz", &kinematics[5]);
	tree.SetBranchAddress("sB_trk1_px", &kinematics[6]);
	tree.SetBranchAddress("sB_trk1_py", &kinematics[7]);
	tree.SetBranchAddress("sB_trk1_pz", &kinematics[8]);
	tree.SetBranchAddress("sB_trk2_px", &kinematics[9]);
	tree.SetBranchAddress("sB_trk2_py", &kinematics[10]);
	tree.SetBranchAddress("sB_trk2_pz", &kinematics[11]);
}

void read_kinematics(Kinematics& kinematics,
					 TLorentzVector& m1, TLorentzVector& m2,
					 TLorentzVector& k1, TLorentzVector& k2,
					 TLorentzVector& p1, TLorentzVector& p2,
					 TLorentzVector& i1, TLorentzVector& i2,
					 TLorentzVector& jpsi) {
	float const m_muon = 0.1056583755;	// +- 0.0000000023	GeV
	float const m_kaon = 0.493677;		// +- 0.000016		GeV
	float const m_pion = 0.13957039;	// +- 0.00000018	GeV
	float const m_prot = 0.93827208816; // +- 0.00000000029 GeV
	float const m_jpsi = 3.096916;		// +- 0.000011		GeV

	for (auto& v: kinematics) {
		v /= 1000; // GeV
	}

	m1.SetXYZM(kinematics[0 + 0], kinematics[0 + 1], kinematics[0 + 2], m_muon);
	m2.SetXYZM(kinematics[3 + 0], kinematics[3 + 1], kinematics[3 + 2], m_muon);
	k1.SetXYZM(kinematics[6 + 0], kinematics[6 + 1], kinematics[6 + 2], m_kaon);
	k2.SetXYZM(kinematics[9 + 0], kinematics[9 + 1], kinematics[9 + 2], m_kaon);
	p1.SetXYZM(kinematics[6 + 0], kinematics[6 + 1], kinematics[6 + 2], m_prot);
	p2.SetXYZM(kinematics[9 + 0], kinematics[9 + 1], kinematics[9 + 2], m_prot);
	i1.SetXYZM(kinematics[6 + 0], kinematics[6 + 1], kinematics[6 + 2], m_pion);
	i2.SetXYZM(kinematics[9 + 0], kinematics[9 + 1], kinematics[9 + 2], m_pion);
	jpsi = m1 + m2;
}

void pk() {
	std::vector<float> lambda { 1380, 1405, 1520, 1600, 1670, 1690, 1710, 1800, 1810, 1820, 1830, 1890, 2000, 2050, 2070, 2080, 2085, 2100, 2110, 2325, 2350, 2585 };
	for (auto& l: lambda) {
		l /= 1000; // GeV
		std::cout << l << std::endl;
	}
	Kinematics kinematics;
	float c1 = 0, c2 = 0;
	float w = 0;
	float chi2 = 0;
	float lxy = 0;
	float pt = 0, spt = 0;

	TFile f_in("datasets/run2_1quarter.root");
	TFile f_out("masses/jkpk.root", "recreate");
	TTree& background = *(TTree*)f_in.Get("stree");

	enable_kinematics(background);
	setup_kinematics(background, kinematics);
	background.SetBranchStatus ("sB_trk1_charge", 1);
	background.SetBranchAddress("sB_trk1_charge", &c1);
	background.SetBranchStatus ("sB_trk2_charge", 1);
	background.SetBranchAddress("sB_trk2_charge", &c2);
	background.SetBranchStatus ("sB_Bs_chi2_ndof", 1);
	background.SetBranchAddress("sB_Bs_chi2_ndof", &chi2);
	background.SetBranchStatus ("sB_Bs_chi2_ndof", 1);
	background.SetBranchAddress("sB_Bs_chi2_ndof", &chi2);
	background.SetBranchStatus ("sB_Lxy_MinA0", 1);
	background.SetBranchAddress("sB_Lxy_MinA0", &lxy);
	background.SetBranchStatus ("sB_Bs_pt", 1);
	background.SetBranchAddress("sB_Bs_pt", &pt);
	background.SetBranchStatus ("sB_MinA0SumPt", 1);
	background.SetBranchAddress("sB_MinA0SumPt", &spt);

	TLorentzVector m1, m2, k1, k2, p1, p2, pi1, pi2, jpsi, pk, jpsik;

	auto const n_selections = 100;
	auto const n_background = background.GetEntries();
	TH1F h_pk("pk", "pk", n_selections, 1.1, 4);
	TH1F h_jk("jk", "jk", n_selections, 3.2, 5.8);
	TH1F h_jkp("jkp", "jkp", n_selections, 3.2, 5.8);
	TH1F h_jkm("jkm", "jkm", n_selections, 3.2, 5.8);
	for (size_t i = 0; i < n_background; ++i) {
		background.GetEntry(i);
		read_kinematics(kinematics, m1, m2, k1, k2, p1, p2, pi1, pi2, jpsi);
		const auto w = c1 * c2 > 0 ? -1 : 1;
		if (lxy < 0.85) continue;
		if (pt / (spt + pt) < 0.2) continue;
		if (chi2 < 1.5) continue;
		auto lb1 = jpsi + k1 + p2;
		auto lb2 = jpsi + k2 + p1;
		// if (lb1.M() < 5.590 || lb1.M() > 5.620) continue;
		// if (lb2.M() < 5.590 || lb2.M() > 5.620) continue;
		// if (c1 * c2 < 0) continue;
		jpsik = jpsi + k1;
		h_jk.Fill(jpsik.M(), w);
		h_jkp.Fill(jpsik.M(), w);
		jpsik = jpsi + k2;
		h_jk.Fill(jpsik.M(), w);
		h_jkm.Fill(jpsik.M(), w);
		pk = k1 + p2; h_pk.Fill(pk.M());
		if (i % 100000 == 0) {
			std::cout << '#';
			std::cout.flush();
		}
	}
	std::cout << std::endl;

	h_pk.Write();
	h_jk.Write();
	h_jkp.Write();
	h_jkm.Write();

	TCanvas canvas;

	h_pk.Draw();
	for (auto l: lambda) {
		TLine line(l, 0, l, n_background / n_selections);
		line.Draw();
	}
	canvas.SaveAs("masses/PKs.png");
	canvas.SaveAs("masses/PKs.C");

	h_jk.Draw();
	for (auto l: lambda) {
		TLine line(l, 0, l, n_background / n_selections);
		line.Draw();
	}
	canvas.SaveAs("masses/JK.png");
	canvas.SaveAs("masses/JK.C");

	h_jkp.Draw();
	for (auto l: lambda) {
		TLine line(l, 0, l, n_background / n_selections);
		line.Draw();
	}
	canvas.SaveAs("masses/JKPlus.png");
	canvas.SaveAs("masses/JKPlus.C");

	h_jkm.Draw();
	for (auto l: lambda) {
		TLine line(l, 0, l, n_background / n_selections);
		line.Draw();
	}
	canvas.SaveAs("masses/JKMinus.png");
	canvas.SaveAs("masses/JKMinus.C");
}
