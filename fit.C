#define _USE_MATH_DEFINES
#include <cmath>
#include <TMath.h>
#include <TFormula.h>
#include <TF2.h>
#include <TH2F.h>

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
	double const xx = (x[0] + x[1]) * M_SQRT1_2;
	double const yy = (x[1] - x[0]) * M_SQRT1_2;
	double b = (yy - (a[3] + a[4] * xx)) / (a[5] + a[6] * xx);
	return exp(a[0] + a[1] * xx + a[2] * xx * xx) * exp(b * b);
}

TF2 f_2d("f_2d", background, 0.0, 2.0, 0.0, 2.0, 7, 2);
double likelihood(double const* parameters) {
	f_2d.SetParameters(parameters + 8);
	const auto N = h_jpk->GetXaxis()->GetNbins();
	const auto M = h_jpk->GetYaxis()->GetNbins();
	double L = 0;
	for (size_t i = 1; i < N; ++i) {
		for (size_t j = 1; j < M; ++j) {
			double const value = h_jpk->GetBinContent(i, j);
			double const v_skk = h_skk->GetBinContent(i, j) * parameters[0];
			double const v_ski = h_ski->GetBinContent(i, j) * parameters[1];
			double const v_sii = h_sii->GetBinContent(i, j) * parameters[2];
			double const v_dkk = h_dkk->GetBinContent(i, j) * parameters[3];
			double const v_dki = h_dki->GetBinContent(i, j) * parameters[4];
			double const v_dii = h_dii->GetBinContent(i, j) * parameters[5];
			double const v_lpk = h_lpk->GetBinContent(i, j) * parameters[6];
			double const v_lpi = h_lpi->GetBinContent(i, j) * parameters[7];
			double const x = h_jpk->GetXaxis()->GetBinCenter(i) - 5.0;
			double const y = h_jpk->GetYaxis()->GetBinCenter(j) - 5.0;
			double const integral = f_2d.Eval(x, y) +
				v_skk + v_ski + v_sii +
				v_dkk + v_dki + v_dii +
				v_lpk + v_lpi;
			L += value * log(integral) - integral;
		}
	}
	return -L;
}

void fit() {
	TFile f_jpk("filtered/datasets/run2_1quarter.root");
	TFile f_skk("filtered/datasets/mc_BS_to_JPSI_K_K.root");
	TFile f_ski("filtered/datasets/mc_BS_to_JPSI_K_PI.root");
	TFile f_sii("filtered/datasets/mc_BS_to_JPSI_PI_PI.root");
	TFile f_dkk("filtered/datasets/mc_BD_to_JPSI_K_K.root");
	TFile f_dki("filtered/datasets/mc_BD_to_JPSI_K_PI.root");
	TFile f_dii("filtered/datasets/mc_BD_to_JPSI_PI_PI.root");
	TFile f_lpk("filtered/datasets/mc_LAMBDA0B_to_JPSI_P_K.root");
	TFile f_lpi("filtered/datasets/mc_LAMBDA0B_to_JPSI_P_PI.root");

	h_jpk = (TH2F*)f_jpk.Get("jpk2d");
	h_skk = (TH2F*)f_skk.Get("jpk2d");
	h_ski = (TH2F*)f_ski.Get("jpk2d");
	h_sii = (TH2F*)f_sii.Get("jpk2d");
	h_dkk = (TH2F*)f_dkk.Get("jpk2d");
	h_dki = (TH2F*)f_dki.Get("jpk2d");
	h_dii = (TH2F*)f_dii.Get("jpk2d");
	h_lpk = (TH2F*)f_lpk.Get("jpk2d");
	h_lpi = (TH2F*)f_lpi.Get("jpk2d");
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
	m.SetMaxFunctionCalls(1e3);
	m.SetMaxIterations(1e3);
	m.SetTolerance(1e-1);
	m.SetErrorDef(0.5);
	m.SetPrintLevel(1);

	double step[]  { 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e-3, 1e-3,  1e-4, 1e-3,  1e-4, 1e-4, 1e-4 };
	double start[] { 4e4, 1e4, 3e4, 5e4, 6e4, 3e4, 5e4, 5e3,    6,   -2, -2e-2,    1, -3e-1,    1, 2e-1 };
	char const* name[] {
		"c_skk", "c_ski", "c_sii",
		"c_dkk", "c_dki", "c_dii",
		"c_lpk", "c_lpi",
		"a", "a_x", "a_x^2", "b", "bx", "c", "cx"
	};

	ROOT::Math::Functor f(&likelihood, 15);
	m.SetFunction(f);
	for (size_t i = 0; i < 15; ++i) {
		m.SetVariable(i, name[i], start[i], step[i]);
		if (i == 2 || i == 4 || i == 7) {
			// m.SetFixedVariable(i, name[i], 0);
		}
		if (i > 7) {
			// m.SetFixedVariable(i, name[i], 1);
		}
	}
	m.Minimize();

	auto const* result = m.X();
	// result = start;
	f_2d.SetParameters(result + 8);
	for (size_t i = 0; i < 15; ++i) {
		std::cout << (i + 1) << ") " << name[i] << ' ' << result[i] << std::endl;
	}

	h_skk->Scale(result[0]);
	h_ski->Scale(result[1]);
	h_sii->Scale(result[2]);
	h_dkk->Scale(result[3]);
	h_dki->Scale(result[4]);
	h_dii->Scale(result[5]);
	h_lpk->Scale(result[6]);
	h_lpi->Scale(result[7]);

	auto x_jpk = h_jpk->ProjectionX();
	auto x_skk = h_skk->ProjectionX();
	auto x_ski = h_ski->ProjectionX(); // 0
	auto x_sii = h_sii->ProjectionX();
	auto x_dkk = h_dkk->ProjectionX(); // 0
	auto x_dki = h_dki->ProjectionX();
	auto x_dii = h_dii->ProjectionX();
	auto x_lpk = h_lpk->ProjectionX();
	auto x_lpi = h_lpi->ProjectionX(); // 0

	auto y_jpk = h_jpk->ProjectionY();
	auto y_skk = h_skk->ProjectionY();
	auto y_ski = h_ski->ProjectionY();
	auto y_sii = h_sii->ProjectionY();
	auto y_dkk = h_dkk->ProjectionY();
	auto y_dki = h_dki->ProjectionY();
	auto y_dii = h_dii->ProjectionY();
	auto y_lpk = h_lpk->ProjectionY();
	auto y_lpi = h_lpi->ProjectionY();

	TCanvas c;
	x_jpk->Draw();
	x_skk->Draw("plc,same");
	x_ski->Draw("plc,same");
	x_sii->Draw("plc,same");
	x_dkk->Draw("plc,same");
	x_dki->Draw("plc,same");
	x_dii->Draw("plc,same");
	x_lpk->Draw("plc,same");
	x_lpi->Draw("plc,same");
	c.SaveAs("filtered/fit_x.png");
	c.SaveAs("filtered/fit_x.C");

	y_jpk->Draw();
	y_skk->Draw("plc,same");
	y_ski->Draw("plc,same");
	y_sii->Draw("plc,same");
	y_dkk->Draw("plc,same");
	y_dki->Draw("plc,same");
	y_dii->Draw("plc,same");
	y_lpk->Draw("plc,same");
	y_lpi->Draw("plc,same");
	c.SaveAs("filtered/fit_y.png");
	c.SaveAs("filtered/fit_y.C");

}
