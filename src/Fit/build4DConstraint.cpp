#include "CPTimeFit.h"

void CPTimeFit::build4DConstraint(std::vector<std::string> vars){

	RooArgList means = RooArgList();
	RooArgList params = RooArgList();
	TString constraint_name = "constraint";
	std::string dir_name = "";
	for(auto var: vars){
		TString mean_name = "mean_" + var;
		RooRealVar *mean = new RooRealVar(mean_name.Data(), mean_name.Data(), 1.0);
		means.add(*mean);
		RooRealVar *obs = (RooRealVar*)ws_.obj(var.c_str());
		if(second_fit_when_comparing_fits_) obs = (RooRealVar*)ws_for_alt_pdf_.obj(var.c_str());
		params.add(*obs);
		constraint_name = constraint_name + "_" + var;
		dir_name = dir_name + "_" + var;
	}

	TMatrixDSym cov_matrix = TMatrixDSym(vars.size());
	for(unsigned int i=0; i<vars.size(); i++){
		for(unsigned int j=0; j<vars.size(); j++){
			for(unsigned int m=0; m<vars.size(); m++){
				for(unsigned int n=0; n<vars.size(); n++){
					if(i != j) cov_matrix[i][j] = 0.0;
					else cov_matrix[i][j] = 0.0003;
				}
			}
		}
	}

	RooMultiVarGaussian constraint(constraint_name.Data(), constraint_name.Data(), params, means, cov_matrix);

	constraint.Print();
	// get the parameter values for the means
	std::string constrain_dir_read = config_dir_ + "constraint" + dir_name + ".txt";
	std::string constrain_dir_write = config_dir_ + "constraint" + dir_name + ".txt.new";
	constraint.getParameters(params)->writeToFile(constrain_dir_write.c_str());
	constraint.getParameters(params)->readFromFile(constrain_dir_read.c_str());

	//next getting the covariance matrix
	std::string cov_dir_read = config_dir_ + "covMatrix" + dir_name + ".txt";
	std::string cov_dir_write = config_dir_ + "covMatrix" + dir_name + ".txt.new";
	std::ofstream writeMatrix(cov_dir_write);
	for(unsigned int i=0; i < vars.size(); i++){
		for(unsigned int j=0; j<vars.size(); j++){
			writeMatrix << cov_matrix[i][j] << "\t";
		}
		writeMatrix << "\n";
	}
	writeMatrix.close();

	std::ifstream readMatrix(cov_dir_read);
	for(unsigned int i=0; i<vars.size(); i++){
		for(unsigned int j=0; j<vars.size(); j++){
			readMatrix >> cov_matrix[i][j];
		}
	}
	readMatrix.close();
	doocore::io::sinfo << "-INFO-  \t" << "Covariance Matrix" << doocore::io::endmsg;
	cov_matrix.Print();

	if(toyfit_){
		RooDataSet *temp;
		RooRandom::randomGenerator()->SetSeed(seed_);
		temp = constraint.generate(params, RooFit::NumEvents(1));
		temp->Print();
		std::cout << "Constraint before sampling the mean values with seed " << seed_ << std::endl;
		constraint.Print();
		for(auto var : vars){
			TString mean_name = "mean_" + var;
			RooRealVar *obs = (RooRealVar*)temp->get()->find(var.c_str());
			RooRealVar *oldMeans = (RooRealVar*)means.find(mean_name);
			std::cout << "Old value of " << oldMeans->GetName() << ": " << oldMeans->getVal() << std::endl;
			oldMeans->setVal(obs->getVal());
			oldMeans->setConstant(true);
			std::cout << "New value of " << oldMeans->GetName() << ": " << oldMeans->getVal() << std::endl;
		}
	}

	constraint.Print();
	// import to workspace
	if(!second_fit_when_comparing_fits_) ws_.import(constraint);
	else ws_for_alt_pdf_.import(constraint);
}

//=======================================================================================================================

void CPTimeFit::build4DConstraint(std::string tagger){

	RooRealVar *Mistag = NULL;
	if(!second_fit_when_comparing_fits_){
      	if(use_single_tag_) Mistag = (RooRealVar*)ws_.obj(single_mistag_var_.c_str());
      	else if(tagger == "OS") Mistag = (RooRealVar*)ws_.obj(os_mistag_var_.c_str());
      	else Mistag = (RooRealVar*)ws_.obj(ss_mistag_var_.c_str());
    }
    else{
    	if(use_single_tag_) Mistag = (RooRealVar*)ws_for_alt_pdf_.obj(single_mistag_var_.c_str());
    	else if (tagger == "OS") Mistag = (RooRealVar*)ws_for_alt_pdf_.obj(os_mistag_var_.c_str());
    	else Mistag = (RooRealVar*)ws_for_alt_pdf_.obj(ss_mistag_var_.c_str());
    }

    TString constraint_name = "constraint";
    TString GLM_name = "Calibration_for_the_" + tagger;
    TString GLM_title = tagger + "_Calibration";

    RooMultiVarGaussian *constraint = NULL;
    RooArgList *params = NULL;
    if(tagger == "OS") {
    	constraint = new RooMultiVarGaussian(glm_os_->covariance_matrix());
    	params = new RooArgList(glm_os_->coefficients());
    	params->add(glm_os_->delta_coefficients());
    }
    else if(tagger == "SS"){
    	constraint = new RooMultiVarGaussian(glm_ss_->covariance_matrix());
    	params = new RooArgList(glm_ss_->coefficients());
    	params->add(glm_ss_->delta_coefficients());
    }
    else{
    	doocore::io::serr << "-ERROR-  \t" << "4D constraint for neither OS nor SS requested - Severe Malfunction!!" <<  doocore::io::endmsg;
    }
	std::vector<std::string> vars;
	for(int i = 0; i < params->getSize(); i++){
		vars.push_back(params->at(i)->GetName());
		constraint_name = constraint_name + "_" + params->at(i)->GetName();
	}

	TMatrixDSym cov_matrix_stat = constraint->covarianceMatrix();
	//extracting statistical uncertainties and writing to txt file
	std::string statUncert_dir_write = config_dir_ + "statUncert_GLM_" + tagger + ".txt";
	std::vector<double> stat_uncert;
	std::ofstream writeUncerts(statUncert_dir_write);
	for(unsigned int i=0; i < vars.size(); i++){
		for(unsigned int j=0; j<vars.size(); j++){
			if(i == j) {
				writeUncerts << sqrt(cov_matrix_stat[i][j]);
				stat_uncert.push_back(sqrt(cov_matrix_stat[i][j]));
			}
		}
		writeUncerts << "\n";
	}
	writeUncerts.close();

	// writing correlation matrix to txt file
	std::string corr_dir_write = config_dir_ + "corrMatrix_GLM_" + tagger + ".txt.new";
	std::ofstream writeCorrMatrix(corr_dir_write);
	for(unsigned int i=0; i < vars.size(); i++){
		for(unsigned int j=0; j<vars.size(); j++){
			writeCorrMatrix << cov_matrix_stat[i][j]/(stat_uncert[i]*stat_uncert[j]) << "\t";
		}
		writeCorrMatrix << "\n";
	}
	writeCorrMatrix.close();

	// reading systematic uncertainties from file
	std::string systUncert_dir_read = config_dir_ + "systUncert_GLM_" + tagger + ".txt";
	std::vector<double> syst_uncert;
	double syst_uncert_temp = 0;
	std::ifstream readsystUncert(systUncert_dir_read);
	for(unsigned int i=0; i<vars.size(); i++){
		readsystUncert >> syst_uncert_temp;
		syst_uncert.push_back(syst_uncert_temp);
	}
	readsystUncert.close();

	// creating covariance matrix from collected information
	TMatrixDSym cov_matrix = TMatrixDSym(vars.size());
	std::string corr_dir_read = config_dir_ + "corrMatrix_GLM_" + tagger + ".txt.new";
	std::ifstream readMatrix(corr_dir_read);
	double correlation = 0;
	for(unsigned int i=0; i<vars.size(); i++){
		for(unsigned int j=0; j<vars.size(); j++){
			readMatrix >> correlation;
			cov_matrix[i][j] = correlation * sqrt(pow(syst_uncert[i], 2) + pow(stat_uncert[i], 2)) * sqrt(pow(syst_uncert[j],2) + pow(stat_uncert[j], 2));
		}
	}
	readMatrix.close();

	doocore::io::sinfo << "-INFO-  \t" << "Covariance Matrix" << doocore::io::endmsg;
	cov_matrix.Print();
	// creating new RooMultiVarGaussian with correct covariance matrix
	TString constraint_name_temp = constraint_name + "_temp";
	RooArgList means = RooArgList();
	int i = 0;
	for(auto var : vars){
		TString mean_name = "mean_" + TString(params->at(i)->GetName());
		RooRealVar *mean = new RooRealVar(*(RooRealVar*)params->at(i), mean_name.Data());
		mean->setConstant(true);
		means.add(*mean);
		i++;
	}
	RooMultiVarGaussian *constraint_new = new RooMultiVarGaussian(constraint_name_temp, constraint_name_temp, *params, means, cov_matrix);
	// this ugly hack with a second RooMultiVarGaussian is needed because the first one from the GLM model is built using RooConstVars instead of
	// RooRealVars. RooConstVars cannot be accessed for some strange reason in the same way...
	RooMultiVarGaussian *gaussianConstraint = NULL;
	if(toyfit_){
		RooDataSet *temp;
		RooRandom::randomGenerator()->SetSeed(seed_);
		temp = constraint_new->generate(*params, RooFit::NumEvents(1));
		temp->Print();
		std::cout << "Constraint before sampling the mean values with seed " << seed_ << std::endl;
		constraint_new->Print();
		TMatrixDSym cov_matrix = constraint_new->covarianceMatrix();
		// RooArgList means = RooArgList();
		int i = 0;
		for(auto var : vars){
			TString mean_name = "mean_" + var;
			RooRealVar *obs = (RooRealVar*)temp->get()->find(var.c_str());
			RooRealVar *oldMeans = (RooRealVar*)means.find(mean_name);
			std::cout << "Old value of " << oldMeans->GetName() << ": " << oldMeans->getVal() << std::endl;
			oldMeans->setVal(obs->getVal());
			oldMeans->setConstant(true);
			std::cout << "New value of " << oldMeans->GetName() << ": " << oldMeans->getVal() << std::endl;
			i += 1;
		}
		gaussianConstraint = new RooMultiVarGaussian(constraint_name, constraint_name, *params, means, cov_matrix);
	}
	else{
		gaussianConstraint = new RooMultiVarGaussian(*constraint_new, constraint_name);
	}
	gaussianConstraint->Print();
	// import to workspace
	if(!second_fit_when_comparing_fits_) ws_.import(*gaussianConstraint, RooFit::RecycleConflictNodes(), RooFit::Silence());
	else ws_for_alt_pdf_.import(*gaussianConstraint, RooFit::RecycleConflictNodes(), RooFit::Silence());
}
