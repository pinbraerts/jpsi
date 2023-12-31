#include "common.hh"

using Graph = std::array<TH1F*, 9>;
using Fit   = std::array<Graph, 3>;
using Back  = std::array<TF1* , 3>;

Fit h;
Back f;

double background(double const* x, double const* a) {
	double const xx = x[0] - 5;
	double const b = (a[3] + a[4] * xx) / (1 + a[5] * xx);
	return exp(a[0] + xx * a[1] + xx * xx * a[2]) * exp(-b * b);
}

double l_step(Graph& h, TF1* f, size_t i, double const* c) {
	double const value = h[0]->GetBinContent(i);
	double const x = h[0]->GetBinCenter(i);
	double integral = f->Eval(x);
	for (size_t j = 1; j < std::size(h); ++j) {
		integral += h[j]->GetBinContent(i) * c[j - 1];
	}
	return value * log(integral) - integral;
}

double likelihood(double const* p) {
	double L = 0;
	for (size_t j = 0; j < std::size(h); ++j) {
		f[j]->SetParameters(p + 8 + j * f[j]->GetNpar());
		auto const N = h[j][0]->GetNbinsX();
		for (size_t i = 1; i <= N; ++i) {
			L += l_step(h[j], f[j], i, p);
		}
	}
	return -L;
}

void fit_jkk(bool fix = false) {
	TFile* file[] = {
		new TFile("filtered/datasets/run2_1quarter.root"),
		new TFile("filtered/datasets/mc_BS_to_JPSI_K_K.root"),
		new TFile("filtered/datasets/mc_BS_to_JPSI_K_PI.root"),
		new TFile("filtered/datasets/mc_BS_to_JPSI_PI_PI.root"),
		new TFile("filtered/datasets/mc_BD_to_JPSI_K_K.root"),
		new TFile("filtered/datasets/mc_BD_to_JPSI_PI_PI.root"),
		new TFile("filtered/datasets/mc_BD_to_JPSI_K_PI.root"),
		new TFile("filtered/datasets/mc_LAMBDA0B_to_JPSI_P_K.root"),
		new TFile("filtered/datasets/mc_LAMBDA0B_to_JPSI_P_PI.root"),
	};
	auto f_output = new TFile("output/results_jkk.root", "recreate");

	char const* label[] = {
		"Data",
		// "Combinatorial BG",
		"B_s \\rightarrow J/\\psi K K",
		"B_s \\rightarrow J/\\psi K \\pi",
		"B_s \\rightarrow J/\\psi \\pi \\pi",
		"B_d \\rightarrow J/\\psi K K",
		"B_d \\rightarrow J/\\psi \\pi \\pi",
		"B_d \\rightarrow J/\\psi K \\pi",
		"\\Lambda^0_b \\rightarrow J/\\psi p K",
		"\\Lambda^0_b \\rightarrow J/\\psi p \\pi",
	};
	std::string const obj_name[] {
		"jkk", "jki", "jik",
	};

	for (size_t j = 0; j < std::size(h); ++j) {
		f[j] = new TF1((obj_name[j] + "_comb_bg").c_str(), background, 5.0, 7.0, 6, 1);
		for (size_t i = 0; i < std::size(h[0]); ++i) {
			h[j][i] = (TH1F*)file[i]->Get(obj_name[j].c_str());
			h[j][i]->SetName((obj_name[j] + label[i]).c_str());
			if (i > 0) {
				auto const integral = h[j][i]->Integral();
				h[j][i]->Scale(1 / integral);
			}
		}
	}

	TMinuitMinimizer m;
	m.SetMaxFunctionCalls(1e7);
	m.SetMaxIterations(1e6);
	m.SetTolerance(1e-5);
	m.SetErrorDef(0.5);
	m.SetPrintLevel(1);

	double step[]  {
		1, 1, 1, 1, 1, 1, 1, 1,
		1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4,
		1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4,
		1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4,
	};
	double start[] {
3086.33, 1074.38, 3556.19, 500.677, 1.15471e-08, 17018.4, 588.073, 2.70347e-08, 
		20.7639, 139.776, -44.7429, 3.73766188, 21.78677276, 1.40820501,
		13.8807, 32.6371, -11.2272, -2.55324514, -10.01506379, 1.21720411,
		22.8181, 61.1985, -19.3603, -3.93336816, -12.39543926, 1.08315551,
	};
	char const* name[] {
		"c_skk", "c_ski", "c_sii", "c_dkk", "c_dii", "c_dki", "c_lpk", "c_lpi",
		"jkk_a", "jkk_a_x", "jkk_a_x^2", "jkk_b", "jkk_bx", "jkk_cx",
		"jki_a", "jki_a_x", "jki_a_x^2", "jki_b", "jki_bx", "jki_cx",
		"jik_a", "jik_a_x", "jik_a_x^2", "jik_b", "jik_bx", "jik_cx",
	};

	ROOT::Math::Functor functor(&likelihood, std::size(start));
	m.SetFunction(functor);
	for (size_t i = 0; i < std::size(start); ++i) {
		if (fix && i == 6) {
			m.SetFixedVariable(i, name[i], start[i]);
		}
		else if (i < 7 && i != 4) {
			m.SetLimitedVariable(i, name[i], start[i], step[i], 1e-8, 1e6);
		}
		else {
			m.SetVariable(i, name[i], start[i], step[i]);
		}
	}
	m.Minimize();

	auto const* result = m.X();
	if (result == nullptr) {
		result = start;
	}
	else {
		TVectorD errs(functor.NDim(), m.Errors());
		errs.Write("errors");
		TVectorD params(std::size(start), result);
		params.Write("constants");
	}

	TH1F* back[std::size(h)];
	TH1F* full[std::size(h)];
	for (size_t j = 0; j < std::size(h); ++j) {
		back[j] = (TH1F*)h[j][0]->Clone((obj_name[j] + "back").c_str());
		back[j]->Clear();
		full[j] = (TH1F*)h[j][0]->Clone((obj_name[j] + "full").c_str());
		full[j]->Clear();
		f[j]->SetParameters(result + 8 + j * f[j]->GetNpar());
		auto const N = h[j][0]->GetNbinsX();
		for (size_t k = 1; k < std::size(h[0]); ++k) {
			h[j][k]->Scale(result[k - 1]);
		}
		for (size_t i = 1; i <= N; ++i) {
			auto const x = h[j][0]->GetBinCenter(i);
			auto value = f[j]->Eval(x);
			back[j]->SetBinContent(i, value);
			for (size_t k = 1; k < std::size(h[0]); ++k) {
				value += h[j][k]->GetBinContent(i);
			}
			full[j]->SetBinContent(i, value);
		}
	}

	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(0);
	std::vector<double> chi2;
	for (size_t j = 0; j < std::size(h); ++j) {
		auto c = new TCanvas(obj_name[j].c_str(), obj_name[j].c_str());
		auto s = new THStack();
		auto l = new TLegend(0.6, 0.3, 0.9, 0.9);
		s->Add(back[j]);
		l->AddEntry(back[j], "Combinatorial BG");
		back[j]->Write();
		full[j]->Write();
		for (size_t k = 0; k < std::size(h[0]); ++k) {
			auto x = h[j][k];
			x->SetFillStyle(3001);
			if (k == 0) {
				x->SetMarkerStyle(20);
				x->SetMarkerSize(1);
				x->SetMarkerColor(kRed);
				x->SetLineColor(kRed);
				x->Draw("E");
				x->SetTitle(obj_name[j].c_str());
				chi2.emplace_back(full[j]->Chi2Test(x));
			}
			else {
				s->Add(x);
			}
			x->Write();
			l->AddEntry(x, label[k]);
		}

		s->Draw("noclear same hist pfc");
		l->Draw();
		c->SetTitle(obj_name[j].c_str());
		c->Draw();
		c->Write();
		c->SaveAs(("output/" + obj_name[j] + ".eps").c_str());
		c->SaveAs(("output/" + obj_name[j] + ".png").c_str());

	}

	TVectorD chi2s(chi2.size(), chi2.data());
	chi2s.Write("chi2s");

	if (m.Status() != 0) {
		std::cerr << "FAILED" << std::endl;
	}

	std::copy(
		result, result + 8,
		std::ostream_iterator<double>(std::cout, ", ")
	);
}
