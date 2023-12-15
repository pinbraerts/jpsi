#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1F.h"

void filter() {
	float const m_muon = 0.1056583755;	// +- 0.0000000023	GeV
	float const m_kaon = 0.493677;		// +- 0.000016		GeV
	float const m_pion = 0.13957039;	// +- 0.00000018	GeV
	float const m_prot = 0.93827208816; // +- 0.00000000029 GeV
	float const m_jpsi = 3.096916;		// +- 0.000011		GeV

	TLorentzVector muon1, muon2, jpsi;
	TLorentzVector part1, part2, prot, k, pi, pipi, pik, pip, ppi, kpi, kk, pp, kp, pk;

	TFile input("datasets/run2_1quarter.root");
	TTree& in = *(TTree*)input.Get("stree");
	TFile output("filtered/run.root", "recreate");
	unsigned M = 100;
	TH1F h_jpsi("jpsi", "jpsi", M, 0, 10);
	TH1F h_pipi("pipi", "pipi", M, 0, 10);
	TH1F h_pik ( "pik",  "pik", M, 0, 10);
	TH1F h_pip ( "pip",  "pip", M, 0, 10);
	TH1F h_kpi ( "kpi",  "kpi", M, 0, 10);
	TH1F h_kk  (  "kk",   "kk", M, 0, 10);
	TH1F h_kp  (  "kp",   "kp", M, 0, 10);
	TH1F h_ppi ( "ppi",  "ppi", M, 0, 10);
	TH1F h_pk  (  "pk",   "pk", M, 0, 10);
	TH1F h_pp  (  "pp",   "pp", M, 0, 10);
	TH1F w_jpsi("jpsiw", "jpsiw", M, 0, 10);
	TH1F w_pipi("pipiw", "pipiw", M, 0, 10);
	TH1F w_pik ( "pikw",  "pikw", M, 0, 10);
	TH1F w_pip ( "pipw",  "pipw", M, 0, 10);
	TH1F w_kpi ( "kpiw",  "kpiw", M, 0, 10);
	TH1F w_kk  (  "kkw",   "kkw", M, 0, 10);
	TH1F w_kp  (  "kpw",   "kpw", M, 0, 10);
	TH1F w_ppi ( "ppiw",  "ppiw", M, 0, 10);
	TH1F w_pk  (  "pkw",   "pkw", M, 0, 10);
	TH1F w_pp  (  "ppw",   "ppw", M, 0, 10);
	// TH2F pipi_pik("pipi_pik", "pipi_pik", M, 0, m_pion + m_pion, M, 0, m_pion + m_kaon);
	// TH2F pipi_kpi("pipi_kpi", "pipi_kpi", M, 0, m_pion + m_pion, M, 0, m_pion + m_kaon);
	// TH2F pipi_kk ("pipi_kk" , "pipi_kk" , M, 0, m_pion + m_pion, M, 0, m_kaon + m_kaon);

	in.SetBranchStatus("*", 0);
	in.SetBranchStatus("sB_mu1_px", 1);
	in.SetBranchStatus("sB_mu1_py", 1);
	in.SetBranchStatus("sB_mu1_pz", 1);
	in.SetBranchStatus("sB_mu2_px", 1);
	in.SetBranchStatus("sB_mu2_py", 1);
	in.SetBranchStatus("sB_mu2_pz", 1);
	in.SetBranchStatus("sB_trk1_px", 1);
	in.SetBranchStatus("sB_trk1_py", 1);
	in.SetBranchStatus("sB_trk1_pz", 1);
	in.SetBranchStatus("sB_trk2_px", 1);
	in.SetBranchStatus("sB_trk2_py", 1);
	in.SetBranchStatus("sB_trk2_pz", 1);
	in.SetBranchStatus("sB_trk1_charge", 1);
	in.SetBranchStatus("sB_trk2_charge", 1);

	Float_t m1x, m1y, m1z, m2x, m2y, m2z, t1x, t1y, t1z, t2x, t2y, t2z;
	Float_t c1 = 0, c2 = 0;
	in.SetBranchAddress( "sB_mu1_px", &m1x);
	in.SetBranchAddress( "sB_mu1_py", &m1y);
	in.SetBranchAddress( "sB_mu1_pz", &m1z);
	in.SetBranchAddress( "sB_mu2_px", &m2x);
	in.SetBranchAddress( "sB_mu2_py", &m2y);
	in.SetBranchAddress( "sB_mu2_pz", &m2z);
	in.SetBranchAddress("sB_trk1_px", &t1x);
	in.SetBranchAddress("sB_trk1_py", &t1y);
	in.SetBranchAddress("sB_trk1_pz", &t1z);
	in.SetBranchAddress("sB_trk2_px", &t2x);
	in.SetBranchAddress("sB_trk2_py", &t2y);
	in.SetBranchAddress("sB_trk2_pz", &t2z);
	in.SetBranchAddress("sB_trk1_charge", &c1);
	in.SetBranchAddress("sB_trk2_charge", &c2);

	auto const N = in.GetEntries();
	for (size_t i = 0; i < N; ++i) {
		in.GetEntry(i);

		m1x /= 1000;
		m1y /= 1000;
		m1z /= 1000;
		m2x /= 1000;
		m2y /= 1000;
		m2z /= 1000;
		t1x /= 1000;
		t1y /= 1000;
		t1z /= 1000;
		t2x /= 1000;
		t2y /= 1000;
		t2z /= 1000;
		muon1.SetXYZM(m1x, m1y, m1z, m_muon); // GeV
		muon2.SetXYZM(m2x, m2y, m2z, m_muon); // GeV
		jpsi = muon1 + muon2;
		h_jpsi.Fill(jpsi.M());

		part1.SetXYZM(t1x, t1y, t1z, m_pion); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_pion); // GeV
		pipi = part1 + part2 + jpsi;

		part1.SetXYZM(t1x, t1y, t1z, m_pion); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_kaon); // GeV
		pik = part1 + part2 + jpsi;

		part1.SetXYZM(t1x, t1y, t1z, m_pion); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_prot); // GeV
		pip = part1 + part2 + jpsi;

		part1.SetXYZM(t1x, t1y, t1z, m_kaon); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_pion); // GeV
		kpi = part1 + part2 + jpsi;

		part1.SetXYZM(t1x, t1y, t1z, m_kaon); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_kaon); // GeV
		kk = part1 + part2 + jpsi;

		part1.SetXYZM(t1x, t1y, t1z, m_kaon); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_prot); // GeV
		kp = part1 + part2 + jpsi;

		part1.SetXYZM(t1x, t1y, t1z, m_prot); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_pion); // GeV
		ppi = part1 + part2 + jpsi;

		part1.SetXYZM(t1x, t1y, t1z, m_prot); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_kaon); // GeV
		pk = part1 + part2 + jpsi;

		part1.SetXYZM(t1x, t1y, t1z, m_prot); // GeV
		part2.SetXYZM(t2x, t2y, t2z, m_prot); // GeV
		pp = part1 + part2 + jpsi;

		auto const w = c1 * c2 < 0 ? 1 : -1;
		h_pipi.Fill(pipi.M());
		h_pik .Fill(pik .M());
		h_pip .Fill(pip .M());
		h_kpi .Fill(kpi .M());
		h_kk  .Fill(kk	.M());
		h_kp  .Fill(kp	.M());
		h_ppi .Fill(ppi .M());
		h_pk  .Fill(pk	.M());
		h_pp  .Fill(pp	.M());
		h_jpsi.Fill(jpsi.M());
		w_pipi.Fill(pipi.M(), w);
		w_pik .Fill(pik .M(), w);
		w_pip .Fill(pip .M(), w);
		w_kpi .Fill(kpi .M(), w);
		w_kk  .Fill(kk	.M(), w);
		w_kp  .Fill(kp	.M(), w);
		w_ppi .Fill(ppi .M(), w);
		w_pk  .Fill(pk	.M(), w);
		w_pp  .Fill(pp	.M(), w);
		w_jpsi.Fill(jpsi.M(), w);

		if (i % 100000 == 0) {
			std::cout << '#';
			std::cout.flush();
		}
	}
	std::cout << std::endl;

	h_pipi.Write();
	h_pik .Write();
	h_pip .Write();
	h_kpi .Write();
	h_kk  .Write();
	h_kp  .Write();
	h_ppi .Write();
	h_pk  .Write();
	h_pp  .Write();
	w_pipi.Write();
	w_pik .Write();
	w_pip .Write();
	w_kpi .Write();
	w_kk  .Write();
	w_kp  .Write();
	w_ppi .Write();
	w_pk  .Write();
	w_pp  .Write();


}
