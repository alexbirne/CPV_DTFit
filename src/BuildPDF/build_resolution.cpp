#include "CPTimeFit.h"

void CPTimeFit::build_resolution(RooRealVar BeautyTime){

    RooCubicSplineFun* acceptance = NULL;
	if(use_acceptance_&&!second_fit_when_comparing_fits_) acceptance = (RooCubicSplineFun*)ws_.obj("acceptance");
    else if(use_acceptance_&&second_fit_when_comparing_fits_) acceptance = (RooCubicSplineFun*)ws_for_alt_pdf_.obj("acceptance");

	RooRealVar resolution_mean = RooRealVar("resolution_mean", "resolution_mean", 0.0);
    resolution_mean.setConstant(true);
    RooRealVar resolution_sigma = RooRealVar("resolution_sigma", "resolution_sigma", 0.05);
    resolution_sigma.setConstant(true);

    RooResolutionModel* resolution = NULL;
    if(use_acceptance_){
        resolution = new RooGaussEfficiencyModel("resolution",
                                                 "resolution",
                                                 BeautyTime,
                                                 dynamic_cast<RooAbsGaussModelEfficiency&>(*acceptance),
                                                 resolution_mean,
                                                 resolution_sigma);
    }
    else{
        resolution = new RooGaussModel("resolution",
                                       "resolution",
                                       BeautyTime,
                                       resolution_mean,
                                       resolution_sigma);
        // for testing purpose
        // resolution = new RooTruthModel("resolution", "resolution", BeautyTime);
    }

	// importing resolution function into workspace
    if(second_fit_when_comparing_fits_) ws_for_alt_pdf_.import(*resolution);
	else ws_.import(*resolution);

}

