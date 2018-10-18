#include "CPTimeFit.h"

void CPTimeFit::Generate_ToySample(int argc, char *argv[]){

	RooAbsPdf *pdf = ws_.pdf("pdfSigTime_conditional");
	RooRealVar SigYield = RooRealVar("SigYield","SigYield",540000,0,1000000);
    RooExtendPdf pdfSigTime_fit = RooExtendPdf("pdfSigTime_fit","pdfSigTime_fit",*pdf,SigYield);

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

	RooWorkspace toy_ws;

	toy_ws.import(pdfSigTime_fit);
	toy_ws.defineSet("observables",observables, true);

	doofit::config::CommonConfig cfg_com("common");
	cfg_com.InitializeOptions(argc, argv);
	doofit::toy::ToyFactoryStdConfig cfg_tfac("toyfac");
	cfg_tfac.InitializeOptions(cfg_com);

	cfg_tfac.set_workspace(&toy_ws);
	cfg_tfac.set_generation_pdf_workspace("pdfSigTime_fit");
	cfg_tfac.set_random_seed(seed_);
	cfg_tfac.set_argset_generation_observables_workspace("observables");

	cfg_com.CheckHelpFlagAndPrintHelp();

	cfg_com.PrintAll();

	doofit::toy::ToyFactoryStd tfac(cfg_com, cfg_tfac);

	RooDataSet* data = tfac.Generate();
	data->Print();

	// for debugging this currently commented lines might be helpful
	std::string seed_string = std::to_string(seed_);
	TString file_name = WriteGeneratedDataFile_;
	// TString file_name = "/net/nfshome/home/abirnkraut/Bd2JpsiKSTest_" + std::to_string(seed_) + ".root";
	TFile file(file_name,"RECREATE");
	RooTreeDataStore *treestore = new RooTreeDataStore(WriteGeneratedDataTree_.c_str(),WriteGeneratedDataTree_.c_str(),observables,*(data->store()));
	(treestore->tree()).Write();
	file.Close();

	// import everything into workspace
	ws_.import(*data,RooFit::Rename("dataset"));
}
