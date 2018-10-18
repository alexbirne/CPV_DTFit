#include "CPTimeFit.h"
#include "RooNLLVar.h"

void CPTimeFit::Fit(int argc, char *argv[], bool plot_only, bool compare_two_fits){

	if(generatetoy_){
		this->Generate_ToySample(argc, argv);
	}

	RooDataSet *data = dynamic_cast<RooDataSet*>(ws_.data("dataset"));
	if(bootstrap_) data = dynamic_cast<RooDataSet*>(ws_.data("bootstrapped_dataset"));
	RooAbsPdf *pdf = NULL;
	if(!second_fit_when_comparing_fits_) pdf = ws_.pdf("pdfSigTime_conditional");
	else pdf= ws_for_alt_pdf_.pdf("pdfSigTime_conditional");

	RooArgSet constraint_argset("constraint_argset");
	RooGaussian *gaussianConstraint;

	// building 4D gaussian constraints for FT calibrations
	RooMultiVarGaussian *multivargaussianConstraint;
	if(!use_GlM_){
		for(auto MultiVarconstraint : MultiVarconstraints_){
			this->build4DConstraint(MultiVarconstraint);
			TString MultiVarconstraint_name = "constraint";
			for(auto name : MultiVarconstraint){
				MultiVarconstraint_name = MultiVarconstraint_name + "_" + name;
			}
			multivargaussianConstraint = (RooMultiVarGaussian*)ws_.pdf(MultiVarconstraint_name);
			if(second_fit_when_comparing_fits_) multivargaussianConstraint = (RooMultiVarGaussian*)ws_for_alt_pdf_.pdf(MultiVarconstraint_name);
			constraint_argset.add(*multivargaussianConstraint);
		}
	}
	else if(use_GlM_ && constrained_tagging_){
		std::string taggers[2] = {"OS", "SS"};
		int i = 0;
		for(auto MultiVarconstraint : MultiVarconstraints_){
			this->build4DConstraint(taggers[i]);
			i += 1;
			TString MultiVarconstraint_name = "constraint";
			for(auto name : MultiVarconstraint){
				MultiVarconstraint_name = MultiVarconstraint_name + "_" + name;
			}
			multivargaussianConstraint = (RooMultiVarGaussian*)ws_.pdf(MultiVarconstraint_name);
			if(second_fit_when_comparing_fits_) multivargaussianConstraint = (RooMultiVarGaussian*)ws_for_alt_pdf_.pdf(MultiVarconstraint_name);
			constraint_argset.add(*multivargaussianConstraint);
		}
	}

	// building gaussian constraints
	for(auto constraint : constraints_){
		this->buildConstraint(constraint);
		TString constraint_name = "constraint_" + constraint;
		gaussianConstraint = (RooGaussian*)ws_.pdf(constraint_name);
		if(second_fit_when_comparing_fits_) gaussianConstraint = (RooGaussian*)ws_for_alt_pdf_.pdf(constraint_name);
		constraint_argset.add(*gaussianConstraint);
	}


	std::string parameter_val_file = config_dir_+"parameters.txt.new";
	pdf->getParameters(data)->writeToFile(parameter_val_file.c_str());
	parameter_val_file = config_dir_+"parameters.txt";
	pdf->getParameters(data)->readFromFile(parameter_val_file.c_str());
	doocore::io::sinfo << "-INFO- \t" << "Full PDF Tree:" << doocore::io::endmsg;
	pdf->printCompactTree();
	doocore::io::sinfo << "-INFO- \t" << "Full PDF:" << doocore::io::endmsg;
	pdf->Print();
	doocore::io::sinfo << "-INFO- \t" << "Full PDF (Verbose output):" << doocore::io::endmsg;
	pdf->Print("V");
	doocore::io::sinfo << "-INFO- \t" << "DataSet:" << doocore::io::endmsg;
	data->Print("V");

	constraint_argset.Print();

	RooArgSet cp_params = RooArgSet();
	if(toyfit_ || reduced_minos_){
		RooRealVar *S = (RooRealVar*)ws_.obj("S_f");
		RooRealVar *Sbar = (RooRealVar*)ws_.obj("S_fbar");
		cp_params.add(*S);
		cp_params.add(*Sbar);
	}

	// Fitting
	RooLinkedList fitting_args;
	// options to control construction of fit
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::NumCPU(numcpu_))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Offset(true))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Extended(false))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::ExternalConstraints(constraint_argset))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Minimizer("Minuit2","migrad"))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Optimize(true))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Hesse(true))));
	if(toyfit_ || reduced_minos_) fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Minos(cp_params))));
	else fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Minos(true))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Save(true))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Strategy(2))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::SumW2Error(false))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::PrintLevel(1))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::Warnings(false))));
	fitting_args.Add((TObject*)(new RooCmdArg(RooFit::PrintEvalErrors(-1))));

	if(!plot_only){
		// these lines prevent the warning of a too large likelihood value when running with constraints;
		RooMsgService::instance().getStream(0).removeTopic(RooFit::MsgTopic::Eval);
		RooMsgService::instance().getStream(1).removeTopic(RooFit::MsgTopic::Eval);
		const RooFitResult* fitresult = NULL;
		if(!DoContour_) fitresult = pdf->fitTo(*data, fitting_args);
		else{
			RooLinkedList nll_args;
			nll_args.Add((TObject*)(new RooCmdArg(RooFit::NumCPU(numcpu_))));
			nll_args.Add((TObject*)(new RooCmdArg(RooFit::Offset(true))));
			nll_args.Add((TObject*)(new RooCmdArg(RooFit::Extended(false))));
			nll_args.Add((TObject*)(new RooCmdArg(RooFit::ExternalConstraints(constraint_argset))));
			nll_args.Add((TObject*)(new RooCmdArg(RooFit::Optimize(true))));

			// creating negative log likelihood
			RooAbsReal* nll = pdf->createNLL(*data, nll_args);
			// creating minimizer/miuit object and configure it
			RooMinimizer m(*nll);
			m.setMinimizerType("Minuit2");
			m.setStrategy(2);
			m.setPrintEvalErrors(-1);
			m.migrad();
			m.hesse();

			RooRealVar *S = (RooRealVar*)ws_.obj("S_f");
			RooRealVar *Sbar = (RooRealVar*)ws_.obj("S_fbar");
			RooArgSet asyms = RooArgSet();
			asyms.add(*S);
			asyms.add(*Sbar);
			m.minos(cp_params);

			fitresult = m.save();
			RooPlot* p1 = m.contour(*S,*Sbar,1.52,2.49,3.44);
			TCanvas canv("canv", "canv", 2000, 1200);
			p1->Draw();
			std::string contour_plot_string = plot_dir_+"contour/contour.pdf";
			canv.SaveAs(contour_plot_string.c_str());
			contour_plot_string = plot_dir_+"contour/contour.tex";
			canv.SaveAs(contour_plot_string.c_str());
			contour_plot_string = plot_dir_+"contour/contour.C";
			canv.SaveAs(contour_plot_string.c_str());
		}
		did_Fit_ = true;

		// including the topics again into the message streams to not miss anything in the following
		RooMsgService::instance().getStream(0).addTopic(RooFit::MsgTopic::Eval);
		RooMsgService::instance().getStream(1).addTopic(RooFit::MsgTopic::Eval);

		doofit::fitter::easyfit::FitResultPrinter fit_result_printer(*fitresult);
		fit_result_printer.set_print_const_pars(true);
		fit_result_printer.set_full_precision_printout(true);
		fit_result_printer.Print();
		TMatrixTSym<double> corMatrix = fitresult->correlationMatrix();
		corMatrix.Print();

		if(toyfit_){
			parameter_val_file = output_dir_ + "FitResults_" + toyIdentifier_ + ".root";
			// std::cout << "Saving FitResult to file " << parameter_val_file << std::endl;
			// TFile fitresultFile(parameter_val_file.c_str(), "RECREATE");
			// fitresult->Write();
			// fitresultFile.Close();
		}
		else{
			parameter_val_file = output_dir_+"FitResults.txt";
			pdf->getParameters(data)->writeToFile(parameter_val_file.c_str());
		}

		if(save_toy_with_toystudy_&&!compare_two_fits){
			doofit::config::CommonConfig cfg_com("common");
			cfg_com.InitializeOptions(argc, argv);
			doofit::toy::ToyFactoryStdConfig cfg_tfac("toyfac");
			cfg_tfac.InitializeOptions(cfg_com);
			doofit::toy::ToyStudyStdConfig cfg_tstudy("toystudy");
			cfg_tstudy.InitializeOptions(cfg_tfac);
			doofit::plotting::PlotConfig cfg_plot("cfg_plot");

			doofit::toy::ToyStudyStd tstudy(cfg_com, cfg_tstudy, cfg_plot);
			tstudy.StoreFitResult(fitresult);
		}
		else if(save_toy_with_toystudy_&&compare_two_fits){
			doofit::config::CommonConfig cfg_com("common");
			cfg_com.InitializeOptions(argc, argv);
			doofit::toy::ToyFactoryStdConfig cfg_tfac("toyfac");
			cfg_tfac.InitializeOptions(cfg_com);
			doofit::toy::ToyStudyStdConfig cfg_tstudy("toystudy");
			cfg_tstudy.InitializeOptions(cfg_tfac);
			doofit::plotting::PlotConfig cfg_plot("cfg_plot");

			doofit::toy::ToyStudyStd tstudy(cfg_com, cfg_tstudy, cfg_plot);
			cfg_com.PrintAll();

			tstudy.StoreFitResult(fitresult, NULL, NULL, NULL, seed_, toy_study_run_id_);
			toy_study_run_id_++;
			second_fit_when_comparing_fits_=true;
		}
		else{
			// save fitresult into workspace
			RooFitResult *fitres = const_cast<RooFitResult*>(fitresult);
			ws_.import(*fitres);  // not sure if this works now, has to try in other functions
		}
	}
}

