#include "CPTimeFit.h"

void CPTimeFit::doBootstrapping(TString weightName){

	if(second_fit_when_comparing_fits_){
		RooDataSet *data = dynamic_cast<RooDataSet*>(ws_.data("bootstrapped_dataset"));
		ws_for_alt_pdf_.import(*data, RooFit::Rename("bootstrapped_dataset"));
	}
	else{
		RooDataSet *data = dynamic_cast<RooDataSet*>(ws_.data("dataset"));

		bool weightedData = data->isWeighted();

		RooArgSet *observables = const_cast<RooArgSet*>(data->get());
		RooDataSet *bootstrapped = NULL;

		if(weightedData){
			RooRealVar weightVar = RooRealVar(weightName.Data(),weightName.Data(),1);
			observables->add(weightVar);
			bootstrapped = new RooDataSet("dataset_Bootstrapped",
										  "dataset_Bootstrapped",
										  *observables,
										  weightName.Data());
		}
		else{
			bootstrapped = new RooDataSet("dataset_Bootstrapped",
										  "dataset_Bootstrapped",
										  *observables);
		}

		int numEntries = numBootstrapEntries_;
		int DataSetEntries = data->numEntries();
		if(numEntries == -1) numEntries = DataSetEntries;

		const RooArgSet *temp;
		TRandom3 rnd(0);

		for(int i=0; i < numEntries; i++){
			int randomNumber = rnd.Integer(DataSetEntries);
			temp = data->get(randomNumber);
			if(weightedData){
				double weight = data->weight();
				bootstrapped->add(*temp, weight, 0);
			}
			else{
				bootstrapped->add(*temp);
			}
		}
		bootstrapped->Print();
		ws_.import(*bootstrapped,RooFit::Rename("bootstrapped_dataset"));
	}




}

