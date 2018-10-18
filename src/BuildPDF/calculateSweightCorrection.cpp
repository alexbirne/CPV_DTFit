#include "CPTimeFit.h"

double CPTimeFit::calculateSweightCorrection(RooDataSet &data, std::string Weight){

	double sum_sw = 0;
	double sum_sw_squared = 0;
	for (int i=0; i<data.numEntries(); i++){
		sum_sw += data.get(i)->getRealValue(Weight.c_str());
		sum_sw_squared += pow(data.get(i)->getRealValue(Weight.c_str()),2);
	}
	double correction = sum_sw/sum_sw_squared;
	doocore::io::sinfo << "-INFO-  \t" << "sweight correction factor = " << correction <<  doocore::io::endmsg;
	return correction;

}

