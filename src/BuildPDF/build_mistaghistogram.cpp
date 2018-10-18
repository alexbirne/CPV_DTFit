#include "CPTimeFit.h"

void CPTimeFit::build_mistaghistogram(std::string tagger){

	TString pdfName = "MistagPDF_" + tagger;

	TString MistagVarName = "Mistag" + tagger;
	if(use_single_tag_) MistagVarName = single_mistag_var_;
	else if(tagger == "OS") MistagVarName = os_mistag_var_;
	else if(tagger == "SS") MistagVarName = ss_mistag_var_;

	if(generatetoy_){
		if(tagger == "OS" && use_os_mistag_template_){
			TFile mistagfile(os_mistag_file_.c_str(), "READ");
			RooWorkspace *mistagws = (RooWorkspace*)mistagfile.Get(os_mistag_ws_.c_str());
			RooHistPdf *mistagpdf = (RooHistPdf*)mistagws->pdf(os_mistag_pdf_.c_str());
			RooHistPdf mistagpdf_os(*mistagpdf, pdfName.Data());
			ws_.import(mistagpdf_os);
		}
		else if(tagger == "SS" && use_ss_mistag_template_){
			TFile mistagfile(ss_mistag_file_.c_str(), "READ");
			RooWorkspace *mistagws = (RooWorkspace*)mistagfile.Get(ss_mistag_ws_.c_str());
			RooHistPdf *mistagpdf = (RooHistPdf*)mistagws->pdf(ss_mistag_pdf_.c_str());
			RooHistPdf mistagpdf_ss(*mistagpdf, pdfName.Data());
			ws_.import(mistagpdf_ss);
		}
		else{
			TString mistag_mean_name = "mistag_mean_" + tagger;
			TString mistag_sigma_name = "mistag_sigma_" + tagger;
			RooRealVar mistag_mean = RooRealVar(mistag_mean_name.Data(), "#eta", 0.37, 0.0, 1.0);
			RooRealVar mistag_sigma = RooRealVar(mistag_sigma_name.Data(), "#eta", 0.02, 0.0, 1.0);
			if(tagger == "SS") {
				mistag_mean.setVal(0.45);
				mistag_sigma.setVal(0.01);
			}
			mistag_mean.setConstant(true);
			mistag_sigma.setConstant(true);
			RooRealVar *Mistag_var = dynamic_cast<RooRealVar*>(ws_.obj(MistagVarName.Data()));
			RooGaussian mistagPdf(pdfName.Data(), pdfName.Data(), *Mistag_var, mistag_mean, mistag_sigma);
			mistagPdf.Print();
			// import mistag pdf into workspace
			ws_.import(mistagPdf);
		}
	}
	else{
		int bin = 50;
		RooDataSet *data = dynamic_cast<RooDataSet*>(ws_.data("dataset"));
		if(bootstrap_) data = dynamic_cast<RooDataSet*>(ws_.data("bootstrapped_dataset"));

		const RooArgSet *obs = data->get();
		RooRealVar *Mistag_var = (RooRealVar*)obs->find(MistagVarName);

		RooDataSet* sliceData;
		// RooHistPdf* mistagPdf;
		TString tagdecision_var = "TagDec" + tagger + "_idx";
		if(tagger == "OS") tagdecision_var = os_tag_var_;
		else if(tagger == "SS") tagdecision_var = ss_tag_var_;
		if(use_single_tag_) tagdecision_var = single_tag_var_;
		sliceData = (RooDataSet*)data->reduce(*obs,Form("((" + tagdecision_var + " == %d) || (" + tagdecision_var + " == %d))", 1, -1));

		TString histoName = "MistagHisto_" + tagger;

		TH1* hist = sliceData->createHistogram(histoName.Data(), *Mistag_var, RooFit::Binning(bin));

		// correcting for numerical reasons some bin contents
		for (int i = 1; i< bin+1; i++)
	    {
    	  double c = hist->GetBinContent(i);
    	  if (c < 1e-37)
    	  {
    	    hist->SetBinContent(i, 1e-37);
    	  }
    	}

    	TString histDataName = "histData_" + tagger;
    	RooDataHist histData(histDataName.Data(), histDataName.Data(), RooArgList(*Mistag_var), hist);

    	RooHistPdf mistagPdf(pdfName.Data(), pdfName.Data(), RooArgSet(*Mistag_var), histData);
    	mistagPdf.Print();
    	// import mistag pdf into workspace
    	if(second_fit_when_comparing_fits_) ws_for_alt_pdf_.import(mistagPdf);
		else ws_.import(mistagPdf);
	}
}

