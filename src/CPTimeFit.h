#include "string"
#include "vector"
#include "fstream"
#include "iostream"
#include "math.h"
//Root
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH2.h"
#include "TMatrixTSym.h"
#include "TString.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "RooAbsPdf.h"
#include "RooBDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooMinimizer.h"
#include "RooMultiVarGaussian.h"
#include "RooNumIntConfig.h"
#include "RooPlot.h"
#include "RooPolyVar.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooTreeDataStore.h"
#include "RooTruthModel.h"
#include "RooUnblindUniform.h"
#include "RooWorkspace.h"

//dooSoftware
#include "doofit/config/CommonConfig.h"
#include "doocore/config/EasyConfig.h"
#include "doocore/io/EasyTuple.h"
#include "doocore/io/MsgStream.h"
#include "doocore/lutils/lutils.h"
#include "doofit/fitter/easyfit/EasyFit.h"
#include "doofit/fitter/easyfit/FitResultPrinter.h"
#include "doofit/plotting/Plot/Plot.h"
#include "doofit/plotting/Plot/PlotConfig.h"
#include "doofit/plotting/correlations/CorrelationPlot.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStd.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStdConfig.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStd.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStdConfig.h"
#include "P2VV/RooCubicSplineFun.h"
#include "P2VV/RooGaussEfficiencyModel.h"

#include "DecRateCoeff/DecRateCoeff_Bd.h"

#include "GLMBuilder.hh"

class CPTimeFit{
public:
    CPTimeFit(std::string identifier = "0");
    void printresult_Latex();
    void Fit(int argc, char *argv[], bool do_fit = true, bool compare_two_fits = false);
    void setData(std::string filename, std::string treename, std::string Cut, std::string Weight, bool truetoytagger);
    void setData(std::string filename, std::string wsname, std::string dataname, std::string Cut, std::string Weight, bool truetoytagger);
    void bootstrapData(bool bootstrap, int numBootstrapEntries = -1);
    void useSingleTag(bool use_single_tag);
    void saveToyWithToyStudy(bool save_toy_with_toystudy);
    void saveGeneratedToyData(bool saveToyData, std::string ToyFile, std::string ToyTree);
    void setOutputDir(std::string output_dir);
    void setPlotDir(std::string plot_dir);
    void setConfigDir(std::string config_dir);
    void setTimeVar(std::string time_var, double range_low, double range_high);
    void setFinalstateVar(std::string finalstate_var);
    void setTagOSVar(std::string os_tag_var);
    void setMistagOSVar(std::string os_mistag_var, double range_low = 0.0, double range_high = 0.5);
    void setTagSSVar(std::string ss_tag_var);
    void setMistagSSVar(std::string ss_mistag_var, double range_low = 0.0, double range_high = 0.5);
    void setSingleTagVar(std::string single_tag_var);
    void setSingleMistagVar(std::string single_mistag_var, double range_low = 0.0, double range_high = 0.5);
    void setFitNumCPU(int numcpu);
    void Plot(bool plot_decaytime, bool plot_correlation, bool plot_oscillation, bool both_finalstates = true);
    void BlindResult(bool ResultBlind);
    void build_PDF();
    void use_acceptance(bool use_acceptance, std::vector<double> acceptance_knots, bool use_B2DXacceptance = false);
    void sweightCorrection(bool swCorr);
    void SetConstraints(std::vector<std::string> constraints);
    void Set4DConstraints(std::vector<std::vector<std::string> > MultiVarconstraint);
    void ToyFit(bool toyfit, int seed = 0);
    void UseReducedMinos(bool reduced_minos);
    void UseGLMCalibration(bool use_glm, bool constrained_tagging = false);
    void setXMLOS(std::string xml_file_os);
    void setXMLSS(std::string xml_file_ss);
    void GenerateToy(bool generatetoy);
    void DoContourPlot(bool doContour);
    void SetMistagPDFOS(bool use_os_mistag_template, std::string os_mistag_file, std::string os_mistag_ws, std::string os_mistag_pdf);
    void SetMistagPDFSS(bool use_ss_mistag_template, std::string ss_mistag_file, std::string ss_mistag_ws, std::string ss_mistag_pdf);
protected:
private:
    void DataHandler(std::string filename,
                     std::string treeorwsname,
                     std::string datasetname,
                     std::string Cut,
                     std::string Weight,
                     bool DataSetInWorkspace,
                     bool truetoytagger);
    void AddTrueToyTagger(RooDataSet& data);
    double Round(double number, int digits);
    void build_acceptance(RooRealVar BeautyTime);
    void build_resolution(RooRealVar BeautyTime);
    void build_mistaghistogram(std::string tagger);
    double calculateSweightCorrection(RooDataSet &data, std::string Weight);
    void doBootstrapping(TString weightName);
    void buildConstraint(std::string var);
    void build4DConstraint(std::vector<std::string> vars);
    void build4DConstraint(std::string tagger);
    void Generate_ToySample(int argc, char *argv[]);
    std::string file_;
    std::string tree_;
    std::string output_dir_;
    std::string plot_dir_;
    std::string config_dir_;
    std::string time_var_;
    double range_low_;
    double range_high_;
    std::string os_tag_var_;
    std::string os_mistag_var_;
    double mistag_os_range_low_;
    double mistag_os_range_high_;
    std::string ss_tag_var_;
    std::string ss_mistag_var_;
    double mistag_ss_range_low_;
    double mistag_ss_range_high_;
    std::string finalstate_var_;
    std::string single_tag_var_;
    std::string single_mistag_var_;
    double single_mistag_range_low_;
    double single_mistag_range_high_;
    std::string xml_file_os_;
    std::string xml_file_ss_;
    bool constrained_tagging_;
    Espresso::GLMBuilder *glm_os_;
    Espresso::GLMBuilder *glm_ss_;
    int numcpu_;
    bool make_plot_;
    bool ResultBlind_;
    bool swCorr_;
    bool bootstrap_;
    int numBootstrapEntries_;
    bool toyfit_;
    bool reduced_minos_;
    bool generatetoy_;
    bool saveGeneratedData_;
    std::string WriteGeneratedDataFile_;
    std::string WriteGeneratedDataTree_;
    bool use_acceptance_;
    bool use_B2DXacceptance_;
    bool use_single_tag_;
    bool save_toy_with_toystudy_;
    bool use_GlM_;
    bool second_fit_when_comparing_fits_;
    bool did_Fit_;
    bool DoContour_;
    std::string tag_var_to_cut_;
    std::string os_mistag_file_;
    std::string os_mistag_ws_;
    std::string os_mistag_pdf_;
    bool use_os_mistag_template_;
    std::string ss_mistag_file_;
    std::string ss_mistag_ws_;
    std::string ss_mistag_pdf_;
    bool use_ss_mistag_template_;
    int seed_;
    int toy_study_run_id_;
    std::string toyIdentifier_;
    std::vector<std::string> constraints_;
    std::vector<std::vector<std::string> > MultiVarconstraints_;
    std::vector<double> acceptance_knots_;
    RooWorkspace ws_;
    RooWorkspace ws_for_alt_pdf_;
};
