#include "CPTimeFit.h"

void CPTimeFit::buildConstraint(std::string var){

	TString mean_name = "mean_" + var;
	TString sigma_name = "sigma_" + var;
	RooRealVar mean = RooRealVar(mean_name.Data(),mean_name.Data(), 1.0);
	RooRealVar sigma = RooRealVar(sigma_name.Data(),sigma_name.Data(), 1.0);

	RooRealVar *obs = (RooRealVar*)ws_.obj(var.c_str());
	if(second_fit_when_comparing_fits_) obs = (RooRealVar*)ws_for_alt_pdf_.obj(var.c_str());
	TString constraint_name = "constraint_" + var;
	RooGaussian constraint(constraint_name.Data(), constraint_name.Data(), *obs, mean, sigma);

	// get the parameter values
	RooArgSet for_params = RooArgSet(*obs);
	std::string constrain_dir_read = config_dir_ + "constraint_" + var + ".txt";
	std::string constrain_dir_write = config_dir_ + "constraint_" + var + ".txt.new";
	constraint.getParameters(for_params)->writeToFile(constrain_dir_write.c_str());
	constraint.getParameters(for_params)->readFromFile(constrain_dir_read.c_str());

	if(toyfit_){
		RooDataSet *temp;
		RooRandom::randomGenerator()->SetSeed(seed_);
		RooArgSet params(*obs);
		temp = constraint.generate(params, RooFit::NumEvents(1));
		temp->Print();
		std::cout << "Constraint before sampling the mean values with seed " << seed_ << std::endl;
		constraint.Print();
		obs = (RooRealVar*)temp->get()->find(var.c_str());
		std::cout << "Old value of " << mean.GetName() << ": " << mean.getVal() << std::endl;
		mean.setVal(obs->getVal());
		mean.setConstant(true);
		std::cout << "New value of " << mean.GetName() << ": " << mean.getVal() << std::endl;
	}

	// import to workspace
	constraint.Print();
	if(!second_fit_when_comparing_fits_) ws_.import(constraint);
	else ws_for_alt_pdf_.import(constraint);
}
