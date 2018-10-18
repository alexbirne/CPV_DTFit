#include "CPTimeFit.h"

void CPTimeFit::AddTrueToyTagger(RooDataSet& data){
	TRandom3 random_gen(seed_);
	// create new RooDataSet with columns of old dataset + cheated tagger
	const RooArgSet *old_obs = data.get();
	RooArgSet new_obs(*old_obs, "new_obs");
	RooCategory CheatedTrueTag = RooCategory("CheatedTrueTag","cheated true tag decision");
	CheatedTrueTag.defineType("Bbar",-1);
	CheatedTrueTag.defineType("untagged",0);
	CheatedTrueTag.defineType("B",+1);
	new_obs.add(CheatedTrueTag);
	RooDataSet dataReady = RooDataSet("dataReady", "dataReady", new_obs);
	for (int i=0; i<data.numEntries(); i++){
		// now getting each step the argset containing the current values/category indices and storing them in a new argset (fucking constness)
		const RooArgSet *temp = data.get(i);
		RooArgSet temp_new(*temp, "temp_new");
		temp_new.add(CheatedTrueTag);
		// based on the mistag and the true tag decision assigning a new toy tagger decision
		double single_mistag = temp_new.getRealValue(single_mistag_var_.c_str());
		int true_tag = temp_new.getCatIndex(single_tag_var_.c_str());
		bool flip_tag = random_gen.Uniform()<=single_mistag;
		temp_new.setCatIndex("CheatedTrueTag", flip_tag ? -1*true_tag : true_tag);
		// adding the new argset to the dataset
		dataReady.add(temp_new);
	}

	// from now on using the cheated true tag decision as single tag
	single_tag_var_ = "CheatedTrueTag";

	// the following lines could be needed for debugging
	// TFile file("/net/nfshome/home/abirnkraut/test.root","RECREATE");
	// RooTreeDataStore *treestore = new RooTreeDataStore("test","test",*dataReady.get(),*(dataReady.store()));
	// (treestore->tree()).Write();
	// file.Close();
	ws_.import(dataReady);
}
