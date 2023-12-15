
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

void selections(TH1& hist, TLorentzVector& mu1, TLorentzVector& mu2, float chi2, float w, float Bl) {
	hist.Fill("total", w);

	// if (mu1.Pt() < 4 || mu2.Pt() < 4) return;
	// hist.Fill("muon_pt", 1);

	// if (mu1.Eta() > 2.5 || mu2.Eta() > 2.5) return;
	// hist.Fill("muon_eta", 1);

	// for (float chi2_r = 1.9; chi2_r >= 0.5; chi2_r -= 0.05) {
	//	if (chi2 > chi2_r) return;
	//	string name = "chi2 < " + to_string(chi2_r);
	//	hist.Fill(name.c_str(), w);
	// }

	for (float l = 0.1; l < 10; l += 0.1) {
		if (Bl < l) return;
		char name[] = "B_Lxy > x.x";
		name[8]  = 0 + char(l);
		name[10] = 0 + char((l - int(l)) * 10);
		hist.Fill(name, w);
	}

}

void select() {
	Kinematics kinematics;
	float c1 = 0, c2 = 0;
	float w = 0;
	float chi2 = 0;

	TFile f_background("datasets/run2_1quarter.root");
	TTree& background = *(TTree*)f_background.Get("stree");
	TFile f_filtered("filtered/charge.root");
	TTree t_charge;

	background.SetBranchStatus ("*", 0);
	background.SetBranchStatus ("sB_trk1_charge", 1);
	background.SetBranchAddress("sB_trk1_charge", &c1);
	background.SetBranchStatus ("sB_trk2_charge", 1);
	background.SetBranchAddress("sB_trk2_charge", &c2);

	auto const n_signal = signal.GetEntries();
	for (size_t i = 0; i < n_signal; ++i) {
		signal.GetEntry(i);
		if (c1 * c2 < 0) return;
	}
	cout << endl;

	auto const n_background = background.GetEntries();
	TH1F h_background("background", "background", n_selections, 0, n_selections);
	for (size_t i = 0; i < n_background; ++i) {
		background.GetEntry(i);
		if (c1 * c2 > 0) continue;
		read_kinematics(kinematics, m1, m2, k1, k2, p1, p2, pi1, pi2, jpsi);
		selections(h_background, m1, m2, chi2, 1);
		if (i % 100000 == 0) {
			cout << '#';
			cout.flush();
		}
	}
	cout << endl;

	std::vector<float> background_supression, signal_acceptance;
	auto const background_total = h_background.GetBinContent(1);
	auto const signal_total = h_signal.GetBinContent(1);
	for (size_t i = 1; i < n_selections; ++i) {
		auto const b = h_background.GetBinContent(i);
		auto const s = h_signal.GetBinContent(i);
		if (s == 0 || b == 0) break;
		auto const bs = 1 - b / background_total;
		auto const sa = s / signal_total;
		background_supression.push_back(bs);
		signal_acceptance.push_back(sa);
		auto const label = h_background.GetXaxis()->GetBinLabel(i);
		cout << label << ':' << s << ' ' << b  << ' ' << sa << ' ' << bs << endl;
	}

	TCanvas canvas;
	TGraph graph(background_supression.size(), signal_acceptance.data(), background_supression.data());
	graph.Draw();
	canvas.SaveAs("acceptance.C");
}
