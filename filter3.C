#define _USE_MATH_DEFINES
#include <cmath>
#include <array>

#include <TVectorD.h>
#include <TMath.h>
#include <TFormula.h>
#include <TF2.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMinuitMinimizer.h>
#include <Math/Functor.h>
#include <Math/IntegratorMultiDim.h>
#include <TStyle.h>
#include <TTree.h>
#include <TLorentzVector.h>

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

void filter3(const char* file = "datasets/run2_1quarter.root") {
	Kinematics kinematics;
	float c1 = 0, c2 = 0;
	float w = 0;
	float chi2 = 0;
	float lxy = 0;
	float pt = 0, spt = 0;
	float pt1 = 0, pt2 = 0;

	TFile f_background(file);
	auto obj = f_background.Get("stree");
	if (obj == nullptr) {
		obj = f_background.Get("BsAllCandidates");
	}
	if (obj == nullptr) {
		std::cout << "no tree in file " << file << std::endl;
		return;
	}
	TTree& background = *(TTree*)obj;
	TFile f_output((std::string("filtered/") + file).c_str(), "recreate");

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
	background.SetBranchStatus ("sB_trk1_pt", 1);
	background.SetBranchAddress("sB_trk1_pt", &pt1);
	background.SetBranchStatus ("sB_trk2_pt", 1);
	background.SetBranchAddress("sB_trk2_pt", &pt2);

	TLorentzVector m1, m2, k1, k2, p1, p2, pi1, pi2, jpsi;

	auto const n_background = background.GetEntries();
	auto const n_selections = 40;
	TH2F h_Jpk("jpk2d", "jpk2d", n_selections, 5.3, 7, n_selections, 5.3, 7);
	TH1F h_jpk("jpk", "jpk", n_selections, 5.3, 7);
	TH1F h_jkp("jkp", "jkp", n_selections, 5.3, 7);
	TH1F h_JKP("jkp_s", "jkp_s", n_selections, 5.3, 7);
	TH1F h_jki("jki", "jki", n_selections, 3.6, 6.6);
	TH1F h_jik("jik", "jik", n_selections, 3.6, 6.6);
	TH1F h_JKI("jki_s", "jki_s", n_selections, 3.6, 6.6);
	TH2F h_Jki("jki2d", "jki2d", n_selections, 3.6, 6.6, n_selections, 3.6, 6.6);
	TH1F h_jkk("jkk", "jkk", n_selections, 4.2, 6.6);
	TH1F h_jii("jii", "jii", n_selections, 3.4, 6.4);
	TH1F h_background("background", "background", n_selections, 0, n_selections);
	for (size_t i = 0; i < n_background; ++i) {
		if (i % 100000 == 0) {
			std::cout << '#';
			std::cout.flush();
		}
		h_background.Fill("total", 1);
		background.GetEntry(i);

		if (chi2 > 1.7) continue;
		h_background.Fill("chi2 < 1.7", 1);
		auto w = c1 * c2 < 0 ? 1 : -1;
		h_background.Fill("different charge", 1);
		read_kinematics(kinematics, m1, m2, k1, k2, p1, p2, pi1, pi2, jpsi);

		if (pt1 < 2000 || pt2 < 2000) continue;
		h_background.Fill("pt > 2", w);

		if (pt / (pt + spt) < 0.2) continue;
		h_background.Fill("pt/spt > 0.2", w);

		if (lxy < 0.85) continue;
		h_background.Fill("Blxy > 0.85", w);

		auto jk1 = jpsi + k1;
		auto jk2 = jpsi + k2;
		if (jk1.M() > 5.1 || jk2.M() > 5.1) continue;
		h_background.Fill("M(jpsi + k) < 5.1", w);

		auto pk = p1 + k2;
		auto kp = p2 + k1;
		if (pk.M() < 2 || kp.M() < 2) continue;
		h_background.Fill("M(p + k) > 2", w);

		auto kpi = k1 + pi2;
		auto pik = pi1 + k2;
		if (kpi.M() < 1.55 || pik.M() < 1.55) continue;
		h_background.Fill("M(pi + k) > 1.55", w);

		auto jpk = jpsi + p1 + k2;
		auto jkp = jpsi + k1 + p2;
		auto jkk = jpsi + k1 + k2;
		auto jii = jpsi + pi1 + pi2;
		auto jki = jpsi + k1 + pi2;
		auto jik = jpsi + pi1 + k2;
		if (jpk.M() < 5.3 || jkp.M() < 5.3) continue;
		h_background.Fill("M(jpsi + p + k) > 5.3", w);

		h_Jki.Fill(jki.M(), jik.M(), w);
		h_Jpk.Fill(jpk.M(), jkp.M(), w);
		h_jpk.Fill(jpk.M(), w);
		h_jkp.Fill(jkp.M(), w);
		h_JKP.Fill(jpk.M(), w);
		h_JKP.Fill(jkp.M(), w);
		h_jkk.Fill(jkk.M(), w);
		h_jii.Fill(jii.M(), w);
		h_jki.Fill(jki.M(), w);
		h_jik.Fill(jik.M(), w);
		h_JKI.Fill(jki.M(), w);
		h_JKI.Fill(jik.M(), w);
	}
	std::cout << std::endl;

	h_background.Write();
	h_Jpk.Write();
	h_jpk.Write();
	h_jkp.Write();
	h_JKP.Write();
	h_jkk.Write();
	h_jii.Write();
	h_jki.Write();
	h_jik.Write();
	h_JKI.Write();
	h_Jki.Write();
	f_output.Close();

	std::cout.precision(3);
	auto const background_total = h_background.GetBinContent(1);
	for (size_t i = 1; i < n_selections; ++i) {
		auto const b = h_background.GetBinContent(i);
		if (b == 0) break;
		auto const label = h_background.GetXaxis()->GetBinLabel(i);
		std::cout << label << ':' << b	<< '/' << background_total << '=' << (b * 100 / background_total) << '%' << std::endl;
	}

}
