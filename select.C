#include "TStyle.h"
#include "common.hh"
#include <string>

void select_step(TH1* hist, float w, int more, int start, int end, int step, float v) {
	hist->Fill(0.0, w);
	for (int i = start; i < end; i += step) {
		if (more) {
			if (v < i / 100.0f) continue;
		}
		else {
			if (v > i / 100.0f) continue;
		}
		hist->Fill(std::to_string(i).c_str(), w);
	}
}

void select() {
	float c1 = 0, c2 = 0;
	float w = 0;
	float chi2 = 0;
	float lxy = 0;
	float pt = 0, spt = 0;

	auto f_background = new TFile("datasets/run2_1quarter.root");
	auto f_signal	  = new TFile("datasets/mc_LAMBDA0B_to_JPSI_P_K.root");
	auto f_output     = new TFile("output/selections.root", "recreate");
	auto background = (TTree*)f_background->Get("stree");
	auto signal	  = (TTree*)f_signal->Get("BsAllCandidates");

	background->SetBranchStatus ("sB_trk1_charge", 1);
	background->SetBranchAddress("sB_trk1_charge", &c1);
	background->SetBranchStatus ("sB_trk2_charge", 1);
	background->SetBranchAddress("sB_trk2_charge", &c2);
	background->SetBranchStatus ("sB_Bs_chi2_ndof", 1);
	background->SetBranchAddress("sB_Bs_chi2_ndof", &chi2);
	background->SetBranchStatus ("sB_Lxy_MinA0", 1);
	background->SetBranchAddress("sB_Lxy_MinA0", &lxy);
	background->SetBranchStatus ("sB_Bs_pt", 1);
	background->SetBranchAddress("sB_Bs_pt", &pt);
	background->SetBranchStatus ("sB_MinA0SumPt", 1);
	background->SetBranchAddress("sB_MinA0SumPt", &spt);

	signal->SetBranchStatus ("sB_w_trigger", 1);
	signal->SetBranchAddress("sB_w_trigger", &w);
	signal->SetBranchStatus ("sB_Bs_chi2_ndof", 1);
	signal->SetBranchAddress("sB_Bs_chi2_ndof", &chi2);
	signal->SetBranchStatus ("sB_Lxy_MinA0", 1);
	signal->SetBranchAddress("sB_Lxy_MinA0", &lxy);
	signal->SetBranchStatus ("sB_Bs_pt", 1);
	signal->SetBranchAddress("sB_Bs_pt", &pt);
	signal->SetBranchStatus ("sB_MinA0SumPt", 1);
	signal->SetBranchAddress("sB_MinA0SumPt", &spt);

	const char* name[] {
		"\\chi^2(H_b)/N_{dof}",
		"L_{xy}(H_{b})",
		"p_t(H_b)/\\sum p_t(track)",
	};
	const int start[] { 100, 70, 17, };
	const int end[] { 191, 120, 30, };
	const int step[] { 10, 5, 1, };

	TH1F* sa[std::size(name)];
	const size_t n_selections = 30;
	for (size_t i = 0; i < std::size(sa); ++i) {
		std::string nm = name[i];
		nm += "\\ signal\\ acceptance";
		sa[i] = new TH1F(nm.c_str(), nm.c_str(), n_selections, start[i], end[i]);
	}
	auto const n_signal = signal->GetEntries();
	for (size_t i = 0; i < n_signal; ++i) {
		if (i % 2000 == 0) {
			std::cout << '#';
			std::cout.flush();
		}
		signal->GetEntry(i);
		select_step(sa[0], w, 0, start[0], end[0], step[0], chi2);
		select_step(sa[1], w, 1, start[1], end[1], step[1], lxy);
		select_step(sa[2], w, 1, start[2], end[2], step[2], pt / (spt + pt));
	}
	std::cout << std::endl;

	TH1F* bs[std::size(name)];
	for (size_t i = 0; i < std::size(bs); ++i) {
		std::string nm = name[i];
		nm += "\\ background\\ supression";
		bs[i] = new TH1F(nm.c_str(), nm.c_str(), n_selections, start[i], end[i]);
	}
	auto const n_background = background->GetEntries();
	for (size_t i = 0; i < n_signal; ++i) {
		if (i % 2000 == 0) {
			std::cout << '#';
			std::cout.flush();
		}
		background->GetEntry(i);
		if (c1 * c2 < 0) continue;
		select_step(bs[0], w, 0, start[0], end[0], step[0], chi2);
		select_step(bs[1], w, 1, start[1], end[1], step[1], lxy);
		select_step(bs[2], w, 1, start[2], end[2], step[2], pt / (spt + pt));
	}
	std::cout << std::endl;

	for (size_t j = 0; j < std::size(name); ++j) {
		auto const t_bs = bs[j]->GetBinContent(0);
		auto const t_sa = sa[j]->GetBinContent(0);
		auto x = bs[j]->GetXaxis();
		std::vector<float> v_bs, v_sa;
		std::vector<int> v_th;
		auto const N = x->GetNbins();
		for (size_t i = 1; i <= N; ++i) {
			auto const b = bs[j]->GetBinContent(i);
			auto const s = sa[j]->GetBinContent(i);
			auto const t = std::atoi(x->GetBinLabel(i));
			if (b == 0 || s == 0) continue;
			auto const bss = 1 - b / t_bs;
			auto const sas = s / t_sa;
			std::cout
				<< name[j] << " | " << (t / 100.0f)
				<< " : " << s << ' ' << b
				<< " = " << sas << ' ' << bss
				<< std::endl;
			v_bs.emplace_back(bss);
			v_sa.emplace_back(sas);
			v_th.emplace_back(t);
		}
		bs[j]->Write();
		sa[j]->Write();

		auto c = new TCanvas(name[j], name[j]);
		auto g = new TGraph(v_bs.size(), v_sa.data(), v_bs.data());
		g->SetTitle(name[j]);
		g->SetName(name[j]);
		auto gx = g->GetXaxis();
		// g->GetXaxis()->SetRangeUser(0, 1);
		gx->SetTitle("signal acceptance");
		auto gy = g->GetYaxis();
		gy->SetRangeUser(0, 1);
		gy->SetTitle("background supression");
		g->SetMarkerStyle(20);
		g->Draw();
		g->Write();

		auto l = new TText();
		l->SetTextSize(0.04);
		for (size_t i = 0; i < v_bs.size(); ++i) {
			char nm[] = "x.xx";
			nm[0] = '0' + v_th[i] / 100;
			nm[2] = '0' + (v_th[i] % 100) / 10;
			nm[3] = '0' + v_th[i] % 10;
			l->DrawText(v_sa[i] + 0.005, v_bs[i], nm);
		}

		auto line = new TLine(0.9, 0, 0.9, 1);
		line->Draw();
		c->Update();
		c->Draw();
		c->SaveAs((std::string("output/") + std::to_string(j) + ".eps").c_str());
		c->Write();
	}

}
