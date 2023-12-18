#include "common.hh"

TH1F* h_jpk;
TH1F* h_skk;
TH1F* h_ski;
TH1F* h_sii;
TH1F* h_dkk;
TH1F* h_dki;
TH1F* h_dii;
TH1F* h_lpk;
TH1F* h_lpi;

double background(double const* x, double const* a) {
	double const xx = x[0] - 5;
	double const b = a[5] != 0 || a[6] != 0
		? (a[3] + a[4] * xx) / (a[5] + a[6] * xx)
		: 0
	;
	return exp(a[0] + xx * a[1] + xx * xx * a[2]) * exp(-b * b);
}

TF1 f_back("f_jkk", background, 5.0, 7.0, 7, 1);
double likelihood(double const* parameters) {
	f_back.SetParameters(parameters + 8);
	auto x = h_jpk->GetXaxis();
	const auto N = x->GetNbins();
	double L = 0;
	for (size_t i = 1; i <= N; ++i) {
		double const value = h_jpk->GetBinContent(i);
		double const v_skk = h_skk->GetBinContent(i) * parameters[0];
		double const v_ski = h_ski->GetBinContent(i) * parameters[1];
		double const v_sii = h_sii->GetBinContent(i) * parameters[2];
		double const v_dkk = h_dkk->GetBinContent(i) * parameters[3];
		double const v_dki = h_dki->GetBinContent(i) * parameters[4];
		double const v_dii = h_dii->GetBinContent(i) * parameters[5];
		double const v_lpk = h_lpk->GetBinContent(i) * parameters[6];
		double const v_lpi = h_lpi->GetBinContent(i) * parameters[7];
		double const xx = x->GetBinCenter(i);
		double const integral = f_back.Eval(xx) +
			v_skk + v_ski + v_sii +
			v_dkk + v_dki + v_dii +
			v_lpk + v_lpi;
		L += value * log(integral) - integral;
	}
	return -L;
}

void fit_jkk() {
	auto f_jpk = new TFile("filtered/datasets/run2_1quarter.root");
	auto f_skk = new TFile("filtered/datasets/mc_BS_to_JPSI_K_K.root");
	auto f_ski = new TFile("filtered/datasets/mc_BS_to_JPSI_K_PI.root");
	auto f_sii = new TFile("filtered/datasets/mc_BS_to_JPSI_PI_PI.root");
	auto f_dkk = new TFile("filtered/datasets/mc_BD_to_JPSI_K_K.root");
	auto f_dki = new TFile("filtered/datasets/mc_BD_to_JPSI_K_PI.root");
	auto f_dii = new TFile("filtered/datasets/mc_BD_to_JPSI_PI_PI.root");
	auto f_lpk = new TFile("filtered/datasets/mc_LAMBDA0B_to_JPSI_P_K.root");
	auto f_lpi = new TFile("filtered/datasets/mc_LAMBDA0B_to_JPSI_P_PI.root");
	auto f_output = new TFile("output/results_jkk.root", "recreate");

	h_jpk = (TH1F*)f_jpk->Get("jkk");
	h_skk = (TH1F*)f_skk->Get("jkk");
	h_ski = (TH1F*)f_ski->Get("jkk");
	h_sii = (TH1F*)f_sii->Get("jkk");
	h_dkk = (TH1F*)f_dkk->Get("jkk");
	h_dki = (TH1F*)f_dki->Get("jkk");
	h_dii = (TH1F*)f_dii->Get("jkk");
	h_lpk = (TH1F*)f_lpk->Get("jkk");
	h_lpi = (TH1F*)f_lpi->Get("jkk");
	h_skk->SetName("skk");
	h_ski->SetName("ski");
	h_sii->SetName("sii");
	h_dkk->SetName("dkk");
	h_dki->SetName("dki");
	h_dii->SetName("dii");
	h_lpk->SetName("lpk");
	h_lpi->SetName("lpi");

	// h_jpk->Scale(1 / h_jpk->Integral());
	h_skk->Scale(1 / h_skk->Integral());
	h_ski->Scale(1 / h_ski->Integral());
	h_sii->Scale(1 / h_sii->Integral());
	h_dkk->Scale(1 / h_dkk->Integral());
	h_dki->Scale(1 / h_dki->Integral());
	h_dii->Scale(1 / h_dii->Integral());
	h_lpk->Scale(1 / h_lpk->Integral());
	h_lpi->Scale(1 / h_lpi->Integral());

	TMinuitMinimizer m;
	m.SetMaxFunctionCalls(1e7);
	m.SetMaxIterations(1e6);
	m.SetTolerance(1e-3);
	m.SetErrorDef(0.5);
	m.SetPrintLevel(1);

	double step[]  {
		1, 1, 1, 1, 1, 1, 1, 1,
		1e-4, 1e-4, 1e-4,
		1e-4, 1e-4, 1e-4, 1e-4,
		// 1e-1, 1e-1,
	};
	double start[] {
		1008, 210, 3116, 1, 11352, 1665, 3060, 1,
		9.6, 2, -19.56,
		-0.02, 0.02, 0.11, 0.22,
		// 1, 1,
	};
	char const* name[] {
		"c_skk",
		"c_ski",
		"c_sii",
		"c_dkk",
		"c_dki",
		"c_dii",
		"c_lpk",
		"c_lpi",
		"a", "a_x", "a_x^2",
		"b", "bx", "c", "cx",
		// "x", "x^",
	};

	ROOT::Math::Functor f(&likelihood, std::size(start));
	m.SetFunction(f);
	for (size_t i = 0; i < std::size(start); ++i) {
		if (i <= 7) {
			m.SetLimitedVariable(i, name[i], start[i], step[i], 0e0, 1e6);
			// m.SetFixedVariable(i, name[i], 2 * start[i]);
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
		TVectorD errs(f.NDim(), m.Errors());
		errs.Write("errors");
		TVectorD params(std::size(start), result);
		params.Write("constants");
	}
	f_back.SetParameters(result + 8);
	f_back.Write("f_background_combinatorial");

	h_skk->Scale(result[0]);
	h_ski->Scale(result[1]);
	h_sii->Scale(result[2]);
	h_dkk->Scale(result[3]);
	h_dki->Scale(result[4]);
	h_dii->Scale(result[5]);
	h_lpk->Scale(result[6]);
	h_lpi->Scale(result[7]);

	auto back = (TH1F*)h_skk->Clone("combinatorial BG");
	auto full = (TH1F*)h_skk->Clone("full BG");
	auto x = full->GetXaxis();
	auto const N = x->GetNbins();
	full->Clear();
	back->Clear();
	for (size_t i = 1; i <= N; ++i) {
		auto const xx = x->GetBinCenter(i);
		double const value = f_back.Eval(xx);
		double const v_skk = h_skk->GetBinContent(i);
		double const v_ski = h_ski->GetBinContent(i);
		double const v_sii = h_sii->GetBinContent(i);
		double const v_dkk = h_dkk->GetBinContent(i);
		double const v_dki = h_dki->GetBinContent(i);
		double const v_dii = h_dii->GetBinContent(i);
		double const v_lpk = h_lpk->GetBinContent(i);
		double const v_lpi = h_lpi->GetBinContent(i);
		double const integral = value + (
			v_skk + v_ski + v_sii +
			v_dkk + v_dki + v_dii +
			v_lpk + v_lpi
		);
		back->SetBinContent(i, value);
		full->SetBinContent(i, integral);
	}
	full->Write("background");
	back->Write("h_background_combinatorial");

	TH1F* histograms[] = {
		h_jpk, back,
		h_skk, h_ski, h_sii,
		h_dkk, h_dii, h_dki,
		h_lpk, h_lpi
	};
	char const* labels[] = {
		"Data",
		"Combinatorial BG",
		"B_s \\rightarrow J/\\psi K K",
		"B_s \\rightarrow J/\\psi K \\pi",
		"B_s \\rightarrow J/\\psi \\pi \\pi",
		"B_d \\rightarrow J/\\psi K K",
		"B_d \\rightarrow J/\\psi \\pi \\pi",
		"B_d \\rightarrow J/\\psi K \\pi",
		"\\Lambda^0_b \\rightarrow J/\\psi p K",
		"\\Lambda^0_b \\rightarrow J/\\psi p \\pi",
	};

	gStyle->SetOptTitle(0);
	// gStyle->SetOptStat(0);
	auto canvas = new TCanvas();
	auto legend = new TLegend(0.6, 0.25, 0.9, 0.9);
	auto stack = new THStack("stack", "stack");
	double chi2 = 0;
	for (size_t i = 0; i < std::size(histograms); ++i) {
		auto x = histograms[i];
		x->SetFillStyle(3001);
		if (i > 0) {
			stack->Add(x);
		}
		else {
			x->SetFillColor(kRed);
			x->SetLineColor(kRed);
			x->Draw("E");
			chi2 = full->Chi2Test(x);
		}
		x->Write(labels[i]);
		legend->AddEntry(x, labels[i]);
	}

	TVectorD chi2s(1, &chi2);
	chi2s.Write("chi2s");
	std::cout << "chi2: " << chi2 << std::endl;

	stack->Draw("noclear same hist pfc");
	legend->Draw();
	canvas->Draw();
	canvas->Write("plot");

}
