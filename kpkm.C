#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

void kpkm() {
	TFile f_in("datasets/mc_LAMBDA0B_to_JPSI_P_K.root");
	auto& tree = *(TTree*)f_in.Get("BsAllCandidates");
	int k_pdg = 0;
	float pt_1 = 0, pt_2 = 0;
	float w = 1;
	tree.SetBranchStatus("*", 0);
	tree.SetBranchStatus ("sB_Kminus_PDG", 1);
	tree.SetBranchAddress("sB_Kminus_PDG", &k_pdg);
	tree.SetBranchStatus ("sB_trk1_pt", 1);
	tree.SetBranchAddress("sB_trk1_pt", &pt_1);
	tree.SetBranchStatus ("sB_trk2_pt", 1);
	tree.SetBranchAddress("sB_trk2_pt", &pt_2);
	tree.SetBranchStatus ("sB_w_Bpt", 1);
	tree.SetBranchAddress("sB_w_Bpt", &w);
	auto const events = tree.GetEntries();
	size_t const bins = 100;

	TFile f_out("filtered/output.root", "recreate");
	auto const p_max = 20000;
	TH1F k_p("k_p", "k_p", bins, 0, p_max);
	TH1F k_n("k_n", "k_n", bins, 0, p_max);
	TH1F p_p("p_p", "p_p", bins, 0, p_max);
	TH1F p_n("p_n", "p_n", bins, 0, p_max);
	// TH2F k  (  "k",   "k", bins, 0, 40000, bins, 0, 4000);
	// TH2F k  (  "p",   "p", bins, 0, 40000, bins, 0, 4000);
	for (size_t i = 0; i < events; ++i) {
		tree.GetEntry(i);
		if (k_pdg > 0) {
			// trk1 = k+
			k_p.Fill(pt_1, w);
			p_n.Fill(pt_2, w);
		}
		else {
			// trk2 = k-
			k_n.Fill(pt_2, w);
			p_p.Fill(pt_1, w);
		}
	}

	k_p.Write();
	k_n.Write();
	p_p.Write();
	p_n.Write();
	f_out.Close();

	TCanvas canvas;
	k_p.SetLineColor(kBlue);
	k_n.SetLineColor(kRed);
	k_p.Draw();
	k_n.Draw("same");
	canvas.SaveAs("k_p_vs_k_n.C");
	canvas.SaveAs("k_p_vs_k_n.png");
	canvas.Clear();
	p_p.SetLineColor(kBlue);
	p_n.SetLineColor(kRed);
	p_p.Draw();
	p_n.Draw("same");
	canvas.SaveAs("p_p_vs_p_n.C");
	canvas.SaveAs("p_p_vs_p_n.png");
}
