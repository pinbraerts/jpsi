#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <TROOT.h>
#include <TCanvas.h>
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
#include <TGraph.h>
#include <TLine.h>

#if __cplusplus < 201703L

namespace std {

template<class T>
constexpr std::size_t size(T const& array) noexcept
{
	return array.size();
}

template<class T, std::size_t N>
constexpr std::size_t size(const T (&array)[N]) noexcept
{
	return N;
}

} // namespace std

#endif // __cplusplus < 201703L
