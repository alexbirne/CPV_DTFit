#include "CPTimeFit.h"

CPTimeFit::CPTimeFit(std::string identifier):
file_(""),
tree_(""),
output_dir_(""),
plot_dir_(""),
config_dir_(""),
time_var_("BeautyTime"),
range_low_(0.2),
range_high_(15),
os_tag_var_("TagDecOS_idx"),
os_mistag_var_("MistagOS"),
ss_tag_var_("TagDecSS_idx"),
ss_mistag_var_("MistagSS"),
mistag_ss_range_low_(0.0),
mistag_ss_range_high_(0.5),
finalstate_var_("BacCharge_idx"),
single_tag_var_("obsTagTrue"),
single_mistag_var_("obsTagOS"),
xml_file_os_(""),
constrained_tagging_(false),
glm_os_(NULL),
glm_ss_(NULL),
numcpu_(1),
make_plot_(true),
ResultBlind_(true),
swCorr_(true),
bootstrap_(false),
numBootstrapEntries_(-1),
toyfit_(false),
reduced_minos_(false),
generatetoy_(false),
saveGeneratedData_(false),
WriteGeneratedDataFile_("~/Test.root"),
WriteGeneratedDataTree_("DecayTree"),
use_acceptance_(true),
use_single_tag_(false),
save_toy_with_toystudy_(false),
use_GlM_(false),
second_fit_when_comparing_fits_(false),
did_Fit_(false),
DoContour_(false),
tag_var_to_cut_(""),
os_mistag_file_(""),
os_mistag_ws_(""),
os_mistag_pdf_(""),
ss_mistag_file_(""),
ss_mistag_ws_(""),
ss_mistag_pdf_(""),
seed_(0),
toy_study_run_id_(0),
toyIdentifier_(identifier),
ws_(RooWorkspace()),
ws_for_alt_pdf_(RooWorkspace())
{
    doocore::io::sinfo << "-INFO-  \t" << "Welcome to CPTimeFit - Tool for sweighted CP fits" << doocore::io::endmsg;

    // not sure if the following lines are necessary - for the moment adapting them from the B2DX fitter
    RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-7);
    RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-7);
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("extrapolation","Wynn-Epsilon");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("maxSteps","1000");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("minSteps","0");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","21Points");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 1000);
    // RooAbsReal::defaultIntegratorConfig()->Print("V");
}

double CPTimeFit::Round(double number, int digits){
    number *= pow(10, digits);
    if(number >= 0) number = floor(number+0.5);
    else number = ceil(number-0.5);
    number /= pow(10,digits);
    return number;
}

void CPTimeFit::BlindResult(bool ResultBlind){
    ResultBlind_ = ResultBlind;
    if(ResultBlind_) doocore::io::sinfo << "-INFO-  \t" << "Result will be blinded using a RooUnblindUniform" <<  doocore::io::endmsg;
    else doocore::io::sinfo << "-INFO-  \t" << "Result will not be blinded, use at own risk!!" <<  doocore::io::endmsg;
}

void CPTimeFit::sweightCorrection(bool swCorr){
    swCorr_ = swCorr;
}

void CPTimeFit::setOutputDir(std::string output_dir){
    output_dir_ = output_dir;
    doocore::io::sinfo << "-INFO-  \t" << "Output directory: " << output_dir_ <<  doocore::io::endmsg;
}

void CPTimeFit::setPlotDir(std::string plot_dir){
    plot_dir_ = plot_dir;
    doocore::io::sinfo << "-INFO-  \t" << "Plot directory: " << plot_dir_ <<  doocore::io::endmsg;
}

void CPTimeFit::setConfigDir(std::string config_dir){
    config_dir_ = config_dir;
    doocore::io::sinfo << "-INFO-  \t" << "Config directory: " << config_dir_ <<  doocore::io::endmsg;
}

void CPTimeFit::setFitNumCPU(int numcpu){
    numcpu_ = numcpu;
}

void CPTimeFit::setTimeVar(std::string time_var, double range_low, double range_high){
    time_var_ = time_var;
    range_low_ = range_low;
    range_high_ = range_high;
    doocore::io::sinfo << "-INFO-  \t" << "Decaytime variable: " << time_var_ << ", Range [" << range_low_ << "," << range_high_ << "]" << doocore::io::endmsg;
}

void CPTimeFit::setFinalstateVar(std::string finalstate_var){
    finalstate_var_ = finalstate_var;
    doocore::io::sinfo << "-INFO-  \t" << "Finalstate variable: " << finalstate_var_ <<  doocore::io::endmsg;
}

void CPTimeFit::setTagOSVar(std::string os_tag_var){
    os_tag_var_ = os_tag_var;
    doocore::io::sinfo << "-INFO-  \t" << "OS tagdecision variable: " << os_tag_var_ <<  doocore::io::endmsg;
}

void CPTimeFit::setMistagOSVar(std::string os_mistag_var, double range_low, double range_high){
    os_mistag_var_ = os_mistag_var;
    mistag_os_range_low_ = range_low;
    mistag_os_range_high_ = range_high;
    doocore::io::sinfo << "-INFO-  \t" << "OS mistag variable: " << os_mistag_var_ <<
    " in range [" << mistag_os_range_low_ << "," << mistag_os_range_high_ << "]" << doocore::io::endmsg;
}

void CPTimeFit::setTagSSVar(std::string ss_tag_var){
    ss_tag_var_ = ss_tag_var;
    doocore::io::sinfo << "-INFO-  \t" << "SS tagdecision variable: " << ss_tag_var_ <<  doocore::io::endmsg;
}

void CPTimeFit::setMistagSSVar(std::string ss_mistag_var, double range_low, double range_high){
    ss_mistag_var_ = ss_mistag_var;
    mistag_ss_range_low_ = range_low;
    mistag_ss_range_high_ = range_high;
    doocore::io::sinfo << "-INFO-  \t" << "SS mistag variable: " << ss_mistag_var_ <<
    " in range [" << mistag_ss_range_low_ << "," << mistag_ss_range_high_ << "]" <<  doocore::io::endmsg;
}

void CPTimeFit::setSingleTagVar(std::string single_tag_var){
    single_tag_var_ = single_tag_var;
    doocore::io::sinfo << "-INFO-  \t" << "Single tagdecision variable: " << single_tag_var_ <<  doocore::io::endmsg;
}

void CPTimeFit::setSingleMistagVar(std::string single_mistag_var, double range_low, double range_high){
    single_mistag_var_ = single_mistag_var;
    single_mistag_range_low_ = range_low;
    single_mistag_range_high_ = range_high;
    doocore::io::sinfo << "-INFO-  \t" << "Single mistag variable: " << single_mistag_var_ <<
    " in range [" << single_mistag_range_low_ << "," << single_mistag_range_high_ << "]" << doocore::io::endmsg;
}

void CPTimeFit::bootstrapData(bool bootstrap, int numBootstrapEntries){
    bootstrap_ = bootstrap;
    numBootstrapEntries_ = numBootstrapEntries;
    if(bootstrap_){
        doocore::io::sinfo << "-INFO-  \t" << "Dataset will be bootstrapped" <<  doocore::io::endmsg;
        doocore::io::sinfo << "-INFO-  \t" << "Using " << numBootstrapEntries_ << " candidates for bootstrapping" <<  doocore::io::endmsg;
    }
}

void CPTimeFit::useSingleTag(bool use_single_tag){
    use_single_tag_ = use_single_tag;
    if(use_single_tag_){
        doocore::io::sinfo << "-INFO-  \t" << "Using single tagger" <<  doocore::io::endmsg;
    }
}

void CPTimeFit::UseGLMCalibration(bool use_glm, bool constrained_tagging){
    use_GlM_ = use_glm;
    constrained_tagging_ = constrained_tagging;
    if(use_GlM_){
        doocore::io::sinfo << "-INFO-  \t" << "Using GLM models for calibration" <<  doocore::io::endmsg;
        doocore::io::sinfo << "-INFO-  \t" << "Make sure to set the appropriate XML files using the following functions:" << doocore::io::endmsg;
        doocore::io::sinfo << "-INFO-  \t" << "CPTimeFit::Set_XML_OS(...)\t CPTimeFit::Set_XML_SS(...)" << doocore::io::endmsg;
        if(constrained_tagging_) doocore::io::sinfo << "-INFO-  \t" << "The tagging parameters will be constrained and not fixed" <<  doocore::io::endmsg;
    }
}

void CPTimeFit::setXMLOS(std::string xml_file_os){
    xml_file_os_ = xml_file_os;
    doocore::io::sinfo << "-INFO-  \t" << "The following XML file is used for the OS tagging: " << xml_file_os_ << doocore::io::endmsg;
}

void CPTimeFit::setXMLSS(std::string xml_file_ss){
    xml_file_ss_ = xml_file_ss;
    doocore::io::sinfo << "-INFO-  \t" << "The following XML file is used for the SS tagging: " << xml_file_ss_ << doocore::io::endmsg;
}

void CPTimeFit::ToyFit(bool toyfit,int seed){
    toyfit_ = toyfit;
    seed_ = seed;
    if(toyfit_){
        doocore::io::sinfo << "-INFO-  \t" << "Fit is performed to toy data - therefore the constraint-means are sampled before fitting" <<  doocore::io::endmsg;
    }
}

void CPTimeFit::UseReducedMinos(bool reduced_minos){
    reduced_minos_ = reduced_minos;
}

void CPTimeFit::GenerateToy(bool generatetoy){
    generatetoy_ = generatetoy;
    if(generatetoy_){
        doocore::io::sinfo << "-INFO-  \t" << "Toy Data will be generated - make sure an appropriate config file is existing" << doocore::io::endmsg;
    }
}

void CPTimeFit::saveToyWithToyStudy(bool save_toy_with_toystudy){
    save_toy_with_toystudy_ = save_toy_with_toystudy;
    if(save_toy_with_toystudy_){
        doocore::io::sinfo << "-INFO-  \t" << "Fit Result will be saved using ToyStudy - make sure an appropriate config file is existing" << doocore::io::endmsg;
    }
}

void CPTimeFit::saveGeneratedToyData(bool saveToyData, std::string ToyFile, std::string ToyTree){
    saveGeneratedData_ = saveToyData;
    if(saveGeneratedData_){
        if(ToyFile != "") WriteGeneratedDataFile_ = ToyFile;
        if(ToyTree != "") WriteGeneratedDataTree_ = ToyTree;
        doocore::io::sinfo << "-INFO-  \t" << "Generated toy data will be saved in file " << WriteGeneratedDataFile_ << " and tree " << WriteGeneratedDataTree_ << doocore::io::endmsg;
    }
    else std::cout << "something is wrong  - check" << std::endl;
}

void CPTimeFit::SetMistagPDFOS(bool use_os_mistag_template, std::string os_mistag_file, std::string os_mistag_ws, std::string os_mistag_pdf){
    use_os_mistag_template_ = use_os_mistag_template;
    os_mistag_file_ = os_mistag_file;
    os_mistag_ws_ = os_mistag_ws;
    os_mistag_pdf_ = os_mistag_pdf;
    if(use_os_mistag_template_){
        doocore::io::sinfo << "-INFO-  \t" << "OS mistag template used from file '" << os_mistag_file_ <<
        "' with workspace name '" << os_mistag_ws_ << "' and PDF name '" << os_mistag_pdf_  << "'." << doocore::io::endmsg;
        doocore::io::swarn << "-INFO-  \t" <<
        "Make sure the Mistag observable you use has the same name as the one used in the loaded mistag PDF" << doocore::io::endmsg;
    }
}

void CPTimeFit::SetMistagPDFSS(bool use_ss_mistag_template, std::string ss_mistag_file, std::string ss_mistag_ws, std::string ss_mistag_pdf){
    use_ss_mistag_template_ = use_ss_mistag_template;
    ss_mistag_file_ = ss_mistag_file;
    ss_mistag_ws_ = ss_mistag_ws;
    ss_mistag_pdf_ = ss_mistag_pdf;
    if(use_ss_mistag_template_){
        doocore::io::sinfo << "-INFO-  \t" << "SS mistag template used from file '" << ss_mistag_file_ <<
        "' with workspace name '" << ss_mistag_ws_ << "' and PDF name '" << ss_mistag_pdf_  << "'." << doocore::io::endmsg;
        doocore::io::swarn << "-INFO-  \t" <<
        "Make sure the Mistag observable you use has the same name as the one used in the loaded mistag PDF" << doocore::io::endmsg;
    }
}

void CPTimeFit::SetConstraints(std::vector<std::string> constraints){
    for(auto constraint : constraints){
        constraints_.push_back(constraint);
    }
}

void CPTimeFit::Set4DConstraints(std::vector<std::vector<std::string> > MultiVarconstraints){
    for(auto MultiVarconstraint : MultiVarconstraints){
        MultiVarconstraints_.push_back(MultiVarconstraint);
    }
}

void CPTimeFit::DoContourPlot(bool doContour){
    DoContour_ = doContour;
    if(DoContour_) doocore::io::sinfo << "-INFO-  \t" << "Contour plot for CP observables Sf vs Sfbar will be done" << doocore::io::endmsg;
}

void CPTimeFit::use_acceptance(bool use_acceptance, std::vector<double> acceptance_knots, bool use_B2DXacceptance){
    use_acceptance_ = use_acceptance;
    use_B2DXacceptance_ = use_B2DXacceptance;
    int i = 0;
    if(!use_acceptance_) doocore::io::sinfo << "-INFO-  \t" << "No decay time acceptance is used" << doocore::io::endmsg;
    else doocore::io::sinfo << "-INFO-  \t" << "The following knots will be used for the decay time acceptance: [";
    for(auto acceptance_knot : acceptance_knots){
        acceptance_knots_.push_back(acceptance_knot);
        if(i==0) doocore::io::sinfo << acceptance_knot;
        else doocore::io::sinfo << ", " << acceptance_knot;
        i++;
    }
    doocore::io::sinfo << "]" << doocore::io::endmsg;
}

void CPTimeFit::setData(std::string filename, std::string treename, std::string Cut, std::string Weight, bool truetoytagger){
    this->DataHandler(filename, treename, "", Cut, Weight, false, truetoytagger);
}

void CPTimeFit::setData(std::string filename, std::string wsname, std::string dataname, std::string Cut, std::string Weight, bool truetoytagger){
    this->DataHandler(filename, wsname, dataname, Cut, Weight, true, truetoytagger);
}
