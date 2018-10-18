#include "CPTimeFit.h"

void CPTimeFit::DataHandler(std::string filename,
							std::string treeorwsname,
							std::string datasetname,
							std::string Cut,
							std::string Weight,
							bool DataSetInWorkspace,
							bool truetoytagger)
{
	if(second_fit_when_comparing_fits_){
		if(generatetoy_){
			RooArgSet observables;
			observables.add(*(dynamic_cast<RooRealVar*>(ws_.obj(time_var_.c_str()))));
			if(use_single_tag_){
				observables.add(*(dynamic_cast<RooRealVar*>(ws_.obj(single_mistag_var_.c_str()))));
				observables.add(*(dynamic_cast<RooCategory*>(ws_.obj(single_tag_var_.c_str()))));
			}
			else{
				observables.add(*(dynamic_cast<RooRealVar*>(ws_.obj(os_mistag_var_.c_str()))));
				observables.add(*(dynamic_cast<RooRealVar*>(ws_.obj(ss_mistag_var_.c_str()))));
				observables.add(*(dynamic_cast<RooCategory*>(ws_.obj(os_tag_var_.c_str()))));
				observables.add(*(dynamic_cast<RooCategory*>(ws_.obj(ss_tag_var_.c_str()))));
			}
			observables.add(*(dynamic_cast<RooCategory*>(ws_.obj(finalstate_var_.c_str()))));
			ws_for_alt_pdf_.import(observables);
		}
		else{
			RooDataSet *data = dynamic_cast<RooDataSet*>(ws_.data("dataset"));
			ws_for_alt_pdf_.import(*data, RooFit::Rename("dataset"));
		}
	}
	else if(DataSetInWorkspace){
		TFile *file = TFile::Open(filename.c_str(),"READ");
		if(file != NULL) std::cout << "Open file: "<< filename << std::endl;
		else {
			std::cout << "Cannot open file: " << filename  << ". Please Check!! (program will crash at some point)" << std::endl;
		}
		RooWorkspace *work = (RooWorkspace*) file->Get(treeorwsname.c_str());
		if( work != NULL ){
			std::cout << "Read workspace: " << treeorwsname << std::endl;
			work->Print();
		}
		else{
			std::cout << "Cannot read workspace. Please Check!! (program will crash at some point)" << std::endl;
		}
		RooDataSet *dataset = dynamic_cast<RooDataSet*>(work->data(datasetname.c_str()));
		doocore::io::sinfo << "-INFO-  \t" << "Dataset:" <<  doocore::io::endmsg;
		dataset->Print();
		RooDataSet *dataReady = NULL;
		if(truetoytagger && use_single_tag_){
			this->AddTrueToyTagger(*dataset);
			dataReady = dynamic_cast<RooDataSet*>(ws_.data("dataReady"));
			doocore::io::sinfo << "-INFO-  \t" << "Dataset with true toy tagger added:" <<  doocore::io::endmsg;
			dataReady->Print();
		}
		else if(truetoytagger && !use_single_tag_){
			doocore::io::swarn << "-INFO-  \t" << "Creation of toy tagger is requested while at the same time multiple taggers shall be used"
			<< doocore::io::endmsg;
			doocore::io::swarn << "-INFO-  \t" << "This does not make sense - therefore toy tagger is not created. Check your configuration!"
			<< doocore::io::endmsg;
			dataReady = new RooDataSet("dataReady", "dataReady", dataset, *(dataset->get()));
		}
		else{
			dataReady = new RooDataSet("dataReady", "dataReady", dataset, *(dataset->get()));
		}
		RooDataSet *reduced_dataset = NULL;
		if(Cut != ""){
			reduced_dataset = (RooDataSet*)dataReady->reduce(Cut.c_str());
			doocore::io::sinfo << "-INFO-  \t" << "Reduced dataset:" <<  doocore::io::endmsg;
			reduced_dataset->Print();
			ws_.import(*reduced_dataset, RooFit::Rename("dataset"));
		}
		else ws_.import(*dataReady,RooFit::Rename("dataset"));
	}
	else{
		// observables
		// time
		RooRealVar BeautyTime = RooRealVar(time_var_.c_str(),"B decay time",range_low_,range_high_,"ps");
		BeautyTime.SetTitle("t_{#kern[-0.1]{B}^{#kern[-0.1]{0}}}");
		BeautyTime.setUnit("ps");
		// finalstate/charge of bachelor particle
		RooCategory BacCharge = RooCategory(finalstate_var_.c_str(),"Finalstate category");
		BacCharge.defineType("D-pi+",+1);
		BacCharge.defineType("D+pi-",-1);
		// mixing state for testing purposes
		RooCategory MixState = RooCategory("catMixingTrue","Mixing category");
		MixState.defineType("unmixed",+1);
		MixState.defineType("mixed",-1);
		// OS tag decision and mistag
		RooCategory* TagDecOS = NULL;
		if(use_single_tag_){
			TagDecOS = new RooCategory(single_tag_var_.c_str(), "single tag decision");
		}
		else{
			TagDecOS = new RooCategory(os_tag_var_.c_str(), "OS tag decision");
		}
		TagDecOS->defineType("Bbar",-1);
		TagDecOS->defineType("B",+1);
		TagDecOS->defineType("untagged",0);

		RooRealVar* MistagOS = NULL;
		if(use_single_tag_){
			MistagOS = new RooRealVar(single_mistag_var_.c_str(), "single mistag",single_mistag_range_low_,single_mistag_range_high_);
		}
		else{
			MistagOS = new RooRealVar(os_mistag_var_.c_str(),"OS mistag",mistag_os_range_low_,mistag_os_range_high_);
		}

		// SS tag decision and mistag
		RooCategory TagDecSS = RooCategory(ss_tag_var_.c_str(),"SS tag decision");
		TagDecSS.defineType("Bbar",-1);
		TagDecSS.defineType("untagged",0);
		TagDecSS.defineType("B",+1);

		RooRealVar MistagSS = RooRealVar(ss_mistag_var_.c_str(),"SS mistag",mistag_ss_range_low_,mistag_ss_range_high_);

		RooRealVar weightVar = RooRealVar(Weight.c_str(), "sWeight", -2, 2);

		std::vector<std::string> cut_vars;
		std::string temp_substr = "";
		std::size_t found_start = 0;
		std::size_t found_end = 0;

		while(true){
			found_start = Cut.find_first_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ", found_end);
			if(found_start == std::string::npos) break;
			found_end = Cut.find_first_of("><=!", found_start);
			temp_substr = Cut.substr(found_start,found_end-found_start);
			// doocore::io::sdebug << temp_substr << doocore::io::endmsg;
			if(!(std::find(cut_vars.begin(), cut_vars.end(),temp_substr) != cut_vars.end())) cut_vars.push_back(temp_substr);
		}

		RooArgSet *observables = NULL;
		if(use_single_tag_){
			observables = new RooArgSet(BeautyTime,
										BacCharge,
										*TagDecOS,
										*MistagOS,
										// MixState,
										weightVar);
		}
		else{
			observables = new RooArgSet(BeautyTime,
										BacCharge,
										*TagDecOS,
										*MistagOS,
										TagDecSS,
										MistagSS,
										weightVar);
		}

		for(auto& cut_var : cut_vars){
			bool test_realVar_in_set = observables->getRealValue(cut_var.c_str(), -10000)==-10000;
			bool test_category_in_set = TString(observables->getCatLabel(cut_var.c_str(), "qazwsx"))==TString("qazwsx");
			if(test_realVar_in_set&&test_category_in_set){
				RooRealVar *temp = new RooRealVar(cut_var.c_str(), cut_var.c_str(), 0);
				observables->add(*temp);
			}
		}

		if(generatetoy_){
			doocore::io::sinfo << "-INFO-  \t" << "Toy Data will be generated - therefore skip reading an external data sample" << doocore::io::endmsg;
			ws_.import(*observables);
		}
		else{
			file_ = filename;
			tree_ = treeorwsname;
			doocore::io::sinfo << "-INFO-  \t" << "File: " << file_ <<  doocore::io::endmsg;
			doocore::io::sinfo << "-INFO-  \t" << "Tree: " << tree_ <<  doocore::io::endmsg;

			doocore::io::EasyTuple etuple(filename,treeorwsname,*observables);
			RooDataSet &data = etuple.ConvertToDataSet(RooFit::Cut(Cut.c_str()));
			doocore::io::sinfo << "-INFO-  \t" << "Unweighted dataset:" <<  doocore::io::endmsg;
			data.Print();
			RooDataSet *dataReady = NULL;
			if(truetoytagger && use_single_tag_){
				this->AddTrueToyTagger(data);
				dataReady = dynamic_cast<RooDataSet*>(ws_.data("dataReady"));
				doocore::io::sinfo << "-INFO-  \t" << "Dataset with true toy tagger added:" <<  doocore::io::endmsg;
				dataReady->Print();
			}
			else if(truetoytagger && !use_single_tag_){
				doocore::io::swarn << "-INFO-  \t" << "Creation of toy tagger is requested while at the same time multiple taggers shall be used"
				<< doocore::io::endmsg;
				doocore::io::swarn << "-INFO-  \t" << "This does not make sense - therefore toy tagger is not created. Check your configuration!"
				<< doocore::io::endmsg;
				dataReady = new RooDataSet("dataReady", "dataReady", &data, *(data.get()));
			}
			else{
				dataReady = new RooDataSet("dataReady", "dataReady", &data, *(data.get()));
			}
			if(Weight == "noWeight"){
				ws_.import(*dataReady,RooFit::Rename("dataset"));
			}
			else{
				doocore::io::sinfo << "-INFO-  \t" << "Weighted dataset (" << Weight << ") will be used" <<  doocore::io::endmsg;
				double correction = 1;
				if(swCorr_) correction = this->calculateSweightCorrection(*dataReady, Weight);

				RooRealVar sWeightCorr = RooRealVar("correctionFactor","correctionFactor",correction);
				RooFormulaVar SW_Corr = RooFormulaVar("sweight","sweight","@0*@1",RooArgList(weightVar,sWeightCorr));
				RooRealVar* weight_corrected = (RooRealVar*) dataReady->addColumn(SW_Corr);
				doocore::io::sinfo << "-INFO-  \t" << "Unweighted dataset" <<  doocore::io::endmsg;
				dataReady->Print();

				TString weightName = Weight.c_str();
				if(swCorr_) weightName = weight_corrected->GetName();

				RooDataSet data_sw(dataReady->GetName(), dataReady->GetTitle(), dataReady, *(dataReady->get()), 0, weightName);
				doocore::io::sinfo << "-INFO-  \t" << "Weighted dataset" <<  doocore::io::endmsg;
				data_sw.Print();

				// import everything into workspace
				ws_.import(data_sw,RooFit::Rename("dataset"));
			}
		}
	}
	if(bootstrap_){
		TString weightName = Weight.c_str();
		this->doBootstrapping(weightName);
	}
}

