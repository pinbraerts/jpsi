#include "common.hh"

TH2F* h_jpk;
TH2F* h_skk;
TH2F* h_ski;
TH2F* h_sii;
TH2F* h_dkk;
TH2F* h_dki;
TH2F* h_dii;
TH2F* h_lpk;
TH2F* h_lpi;

double background(double const* x, double const* a) {
	double const xx = (x[0] + x[1] - 10) * M_SQRT1_2;
	double const yy = (x[1] - x[0]) * M_SQRT1_2;
	double const b = a[5] != 0 || a[6] != 0
		? (yy - (a[3] + a[4] * xx)) / (a[5] + a[6] * xx)
		: 0
	;
	return exp(a[0] + xx * a[1] + xx * xx * a[2]) * exp(-b * b);
}

TF2 f_2d("f_2d", background, 5.0, 7.0, 5.0, 7.0, 7, 2);
double likelihood(double const* parameters) {
	f_2d.SetParameters(parameters + 8);
	auto x = h_jpk->GetXaxis();
	auto y = h_jpk->GetYaxis();
	const auto N = x->GetNbins();
	const auto M = y->GetNbins();
	double L = 0;
	for (size_t i = 1; i <= N; ++i) {
		for (size_t j = 1; j <= M; ++j) {
			double const value = h_jpk->GetBinContent(i, j);
			double const v_skk = h_skk->GetBinContent(i, j) * parameters[0];
			double const v_ski = h_ski->GetBinContent(i, j) * parameters[1];
			double const v_sii = h_sii->GetBinContent(i, j) * parameters[2];
			double const v_dkk = h_dkk->GetBinContent(i, j) * parameters[3];
			double const v_dki = h_dki->GetBinContent(i, j) * parameters[4];
			double const v_dii = h_dii->GetBinContent(i, j) * parameters[5];
			double const v_lpk = h_lpk->GetBinContent(i, j) * parameters[6];
			double const v_lpi = h_lpi->GetBinContent(i, j) * parameters[7];
			double const xx = x->GetBinCenter(i);
			double const yy = y->GetBinCenter(j);
			double integral = f_2d.Eval(xx, yy);
			integral +=
				v_skk + v_ski + v_sii +
				v_dkk + v_dki + v_dii +
				v_lpk + v_lpi;
			L += value * log(integral) - integral;
		}
	}
	return -L;
}

void fit() {
	auto f_jpk = new TFile("filtered/datasets/run2_1quarter.root");
	auto f_skk = new TFile("filtered/datasets/mc_BS_to_JPSI_K_K.root");
	auto f_ski = new TFile("filtered/datasets/mc_BS_to_JPSI_K_PI.root");
	auto f_sii = new TFile("filtered/datasets/mc_BS_to_JPSI_PI_PI.root");
	auto f_dkk = new TFile("filtered/datasets/mc_BD_to_JPSI_K_K.root");
	auto f_dki = new TFile("filtered/datasets/mc_BD_to_JPSI_K_PI.root");
	auto f_dii = new TFile("filtered/datasets/mc_BD_to_JPSI_PI_PI.root");
	auto f_lpk = new TFile("filtered/datasets/mc_LAMBDA0B_to_JPSI_P_K.root");
	auto f_lpi = new TFile("filtered/datasets/mc_LAMBDA0B_to_JPSI_P_PI.root");
	auto f_output = new TFile("output/results.root", "recreate");

	h_jpk = (TH2F*)f_jpk->Get("jpk2d");
	h_skk = (TH2F*)f_skk->Get("jpk2d");
	h_ski = (TH2F*)f_ski->Get("jpk2d");
	h_sii = (TH2F*)f_sii->Get("jpk2d");
	h_dkk = (TH2F*)f_dkk->Get("jpk2d");
	h_dki = (TH2F*)f_dki->Get("jpk2d");
	h_dii = (TH2F*)f_dii->Get("jpk2d");
	h_lpk = (TH2F*)f_lpk->Get("jpk2d");
	h_lpi = (TH2F*)f_lpi->Get("jpk2d");
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
	m.SetMaxFunctionCalls(1e5);
	m.SetMaxIterations(1e4);
	m.SetTolerance(1e-5);
	m.SetErrorDef(0.5);
	m.SetPrintLevel(1);

	double step[]  {
		1e1, 1e1, 1e1, 1e1, 1e1, 1e1, 1e1, 1e1,
			   1,			1,	  1,
		1e-2, 1e-4, 1e-4, 1e-4,
		// 1e-1, 1e-1,
	};
	double start[] {
		6e4, 8e3, 1e4, 2e4, 2e4, 2e4, 2e4, 6e2,
		-140, -3, -1,
		-1e-1, 1e-1, 1e-2, 3e-1,
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
		m.SetVariable(i, name[i], start[i], step[i]);
		if (i <= 7) {
			m.SetVariableLimits(i, 1e0, 1e6);
			// m.FixVariable(i);
			// m.SetVariableLimits(i, 1e2, 1e5);
			// m.SetFixedVariable(i, name[i], 0);
		}
	}
	m.Minimize();

	TVectorD errs(f.NDim(), m.Errors());
	errs.Write("errors");

	auto const* result = m.X();
	f_2d.SetParameters(result + 8);
	f_2d.Write("f_background_combinatorial");

	TVectorD params(std::size(start), result);
	params.Write("constants");

	h_skk->Scale(result[0]);
	h_ski->Scale(result[1]);
	h_sii->Scale(result[2]);
	h_dkk->Scale(result[3]);
	h_dki->Scale(result[4]);
	h_dii->Scale(result[5]);
	h_lpk->Scale(result[6]);
	h_lpi->Scale(result[7]);

	auto back = (TH2F*)h_skk->Clone("combinatorial BG");
	auto full = (TH2F*)h_skk->Clone("full BG");
	auto x = back->GetXaxis();
	auto y = back->GetYaxis();
	auto const N = x->GetNbins();
	auto const M = y->GetNbins();
	back->Clear();
	for (size_t i = 1; i <= N; ++i) {
		for (size_t j = 1; j <= M; ++j) {
			auto const xx = x->GetBinCenter(i);
			auto const yy = y->GetBinCenter(j);
			double const value = f_2d.Eval(xx, yy);
			double const v_skk = h_skk->GetBinContent(i, j);
			double const v_ski = h_ski->GetBinContent(i, j);
			double const v_sii = h_sii->GetBinContent(i, j);
			double const v_dkk = h_dkk->GetBinContent(i, j);
			double const v_dki = h_dki->GetBinContent(i, j);
			double const v_dii = h_dii->GetBinContent(i, j);
			double const v_lpk = h_lpk->GetBinContent(i, j);
			double const v_lpi = h_lpi->GetBinContent(i, j);
			double const integral = value +(
				v_skk + v_ski + v_sii +
				v_dkk + v_dki + v_dii +
				v_lpk + v_lpi);
			full->SetBinContent(i, j, integral);
			back->SetBinContent(i, j, value);
		}
	}
	full->Write("background");
	// back->Write("h_background_combinatorial");

	TH2F* histograms[] = {
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

	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(0);
	auto canvas = new TCanvas("jkk", "jkk", 1000, 500);
	auto legendx = new TLegend(0.55, 0.25, 0.9, 0.9);
	auto legendy = new TLegend(0.55, 0.25, 0.9, 0.9);
	auto stackx = new THStack("stackx", "stackx");
	auto stacky = new THStack("stacky", "stacky");
	auto fullx = full->ProjectionX("", 1, N);
	auto fully = full->ProjectionY("", 1, M);
	double chi2[2];
	canvas->Divide(2);
	for (size_t i = 0; i < std::size(histograms); ++i) {
		histograms[i]->Write(labels[i]);
		canvas->cd(1);
		auto x = histograms[i]->ProjectionX("_px", 1, N);
		x->SetFillStyle(3001);
		if (i > 0) {
			stackx->Add(x);
		}
		else {
			x->SetFillColor(kRed);
			x->SetLineColor(kRed);
			x->Draw("E");
			x->GetYaxis()->SetRangeUser(0, 4100);
			x->SetTitle("M(J/\\psi p K)");
			chi2[0] = fullx->Chi2Test(x);
		}
		x->Write();
		legendx->AddEntry(x, labels[i]);

		canvas->cd(2);
		auto y = histograms[i]->ProjectionY("_py", 1, M);
		y->SetFillStyle(3001);
		if (i > 0) {
			stacky->Add(y);
		}
		else {
			y->SetFillColor(kRed);
			y->SetLineColor(kRed);
			y->Draw("E");
			y->SetTitle("M(J/\\psi K p)");
			chi2[1] = fully->Chi2Test(y);
		}
		y->Write();
		legendy->AddEntry(y, labels[i]);
	}

	TVectorD chi2s(std::size(chi2), chi2);
	chi2s.Write("chi2s");
	std::cout << "chi2_x: " << chi2[0] << ", chi2_y: " << chi2[1] << std::endl;

	canvas->cd(1)->SetTitle("x projection");
	stackx->Draw("noclear same hist pfc");
	legendx->Draw();
	canvas->cd(2)->SetName("y projection");
	stacky->Draw("noclear same hist pfc");
	legendy->Draw();
	canvas->Draw();
	canvas->Write("plot");
	canvas->SaveAs("output/jpk.eps");
	canvas->SaveAs("output/jpk.png");

	if (m.Status() != 0) {
		std::cerr << "FAILED" << std::endl;
	}

	std::copy(
		result, result + 8,
		std::ostream_iterator<double>(std::cout, ", ")
	);
}
