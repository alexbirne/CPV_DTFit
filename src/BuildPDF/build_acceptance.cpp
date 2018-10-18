#include "CPTimeFit.h"

void CPTimeFit::build_acceptance(RooRealVar BeautyTime){

	if(!use_B2DXacceptance_){
		// own acceptance model
		doocore::io::sinfo << "-INFO-  \t" << "Using the default acceptance model" <<  doocore::io::endmsg;
		doocore::io::sinfo << "-INFO-  \t" << "Check if you're fixing one coefficient - e.g. to one" <<  doocore::io::endmsg;
		// empty arg list for coefficients
		int num_knots = acceptance_knots_.size();
		RooArgList* list = new RooArgList();
		RooRealVar *coefficient = NULL;
		for(int i=1; i<=num_knots; i++){
			std::string coefficient_num = std::to_string(i);
			TString coefficient_name = "parCSpline" + coefficient_num;
			TString coefficient_title = "v_{" + coefficient_num + "}";
			coefficient = new RooRealVar(coefficient_name.Data(), coefficient_title.Data(), 1.0, 0. ,3.);
			coefficient->setConstant(false);
			list->add(*coefficient);
			if(i == 1 || i == num_knots) list->add(*coefficient);
		}

		list->Print();
		RooCubicSplineFun spline("acceptance","acceptance",BeautyTime,acceptance_knots_,*list);

		// importing acceptance into workspace
		if(second_fit_when_comparing_fits_) ws_for_alt_pdf_.import(spline);
		else ws_.import(spline);
	}
	else{
		// B2DX acceptance model
		doocore::io::sinfo << "-INFO-  \t" << "Using B2DX acceptance model" <<  doocore::io::endmsg;
		// list of knots
		std::vector<double> knots;
		knots.push_back(0.5);
		knots.push_back(1.0);
		knots.push_back(1.5);
		knots.push_back(2.0);
		knots.push_back(2.3);
		knots.push_back(2.6);
		knots.push_back(3.0);
		knots.push_back(4.0);
		knots.push_back(10.0);

		RooBinning knotbinning = RooBinning(BeautyTime.getMin(), BeautyTime.getMax(), "knotbinning");
		for(auto knot : knots) knotbinning.addBoundary(knot);
		knotbinning.removeBoundary(BeautyTime.getMin());
		knotbinning.removeBoundary(BeautyTime.getMax());
		// not sure why the boundaries are removed twice - for the moment copying this from B2DXfitter,
		// see https://gitlab.cern.ch/lhcb/Urania/blob/dev_PyEspresso_merge/PhysFit/B2DXFitters/python/B2DXFitters/acceptanceutils.py)
		knotbinning.removeBoundary(BeautyTime.getMin());
		knotbinning.removeBoundary(BeautyTime.getMax());

		double lo = BeautyTime.getMin();
		double hi = BeautyTime.getMax();
		BeautyTime.setBinning(knotbinning, "knotbinning");
		BeautyTime.setRange(lo, hi);

		// empty arg list for coefficients
		RooArgList* list = new RooArgList();
		int num_knots = knots.size();
		RooRealVar *coefficient = NULL;
		for(int i=1; i<=num_knots; i++){
			std::string coefficient_num = std::to_string(i);
			TString coefficient_name = "parCSpline" + coefficient_num;
			TString coefficient_title = "v_{" + coefficient_num + "}";
			coefficient = new RooRealVar(coefficient_name.Data(), coefficient_title.Data(), 1.0, 0. ,3.);
			coefficient->setConstant(false);
			list->add(*coefficient);
		}

		RooConstVar *one = new RooConstVar("one", "1", 1.0);
		list->add ( *one );
		knots.push_back(BeautyTime.getMax());

		double fudge = (knots.back() - knots.at(knots.size()-2))/(knots.at(knots.size()-3) - knots.at(knots.size()-2));
		RooConstVar constVar1 = RooConstVar("parCSpline8_0","parCSpline8_0",1-fudge);
		RooConstVar constVar2 = RooConstVar("parCSpline8_1","parCSpline8_1",fudge);
		RooArgList *polylist = new RooArgList(constVar1,constVar2);
		RooPolyVar *polyvar = new RooPolyVar("parCSpline9","v_{9}", *(RooAbsReal*)(list->at(list->getSize() - 2)), *polylist);
		list->add( *polyvar );

		list->Print();
		RooCubicSplineFun spline("acceptance","acceptance",BeautyTime,"knotbinning",*list);

		// importing acceptance into workspace
		if(second_fit_when_comparing_fits_) ws_for_alt_pdf_.import(spline);
		else ws_.import(spline);
	}
}

