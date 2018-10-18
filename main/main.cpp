#include "vector"
#include "CPTimeFit.h"
#include "boost/program_options.hpp"

#include "doocore/config/EasyConfig.h"

int main(int argc, char *argv[]) {

    std::string config_file;
    std::string data_file = "";
    std::string identifier = "0";

    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("config,C", boost::program_options::value<std::string>(), "config file for fitter")
        ("toyfile,t", boost::program_options::value<std::string>(), "datafile which - if given - will be used instead of file in config")
        ("identifier,i", boost::program_options::value<std::string>(), "identifier which - if given - will be used for toys (e.g. as seed)");

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
    boost::program_options::notify(vm);

    if(vm.count("help")){
        std::cout << desc << std::endl;
        return 1;
    }

    if(vm.count("config")){
        config_file = vm["config"].as<std::string>();
    }
    else{
        doocore::io::serr << "-ERROR- \t" << "Configfile not set." << doocore::io::endmsg;
        std::cout << desc << std::endl;
        return 1;
    }

    if(vm.count("toyfile")){
        data_file = vm["toyfile"].as<std::string>();
    }

    if(vm.count("identifier")){
        identifier = vm["identifier"].as<std::string>();
    }

    doocore::config::EasyConfig config(config_file);

    CPTimeFit fitter(identifier);
    fitter.setOutputDir(config.getString("config.output_dir"));
    fitter.setPlotDir(config.getString("config.plot_dir"));
    fitter.setConfigDir(config.getString("config.config_dir"));
    fitter.sweightCorrection(config.getBool("fit_settings.sweightCorrection"));
    fitter.BlindResult(config.getBool("fit_settings.BlindResult"));
    fitter.bootstrapData(config.getBool("fit_settings.doBootstrap"), config.getInt("fit_settings.CandidatesForBootstrap"));
    fitter.useSingleTag(config.getBool("observables.single_tag"));
    fitter.GenerateToy(config.getBool("toy.GenerateToy"));
    if(config.getBool("toy.GenerateToy")){
        fitter.SetMistagPDFOS(config.getBool("toy.MistagOS_Template.use"),
                              config.getString("toy.MistagOS_Template.file"),
                              config.getString("toy.MistagOS_Template.workspace"),
                              config.getString("toy.MistagOS_Template.pdf"));
        fitter.SetMistagPDFSS(config.getBool("toy.MistagSS_Template.use"),
                              config.getString("toy.MistagSS_Template.file"),
                              config.getString("toy.MistagSS_Template.workspace"),
                              config.getString("toy.MistagSS_Template.pdf"));
    }
    fitter.saveToyWithToyStudy(config.getBool("toy.SaveWithToyStudy"));
    fitter.saveGeneratedToyData(config.getBool("toy.SaveGeneratedData.action"),
                                config.getString("toy.SaveGeneratedData.file"),
                                config.getString("toy.SaveGeneratedData.tree"));

    std::string::size_type sz; // converting the string from the command line into an integer as seed for the toyfit
    fitter.ToyFit(config.getBool("toy.FitToy"), std::stoi(identifier, &sz));
    fitter.UseReducedMinos(config.getBool("toy.ReducedMinos"));
    fitter.setFitNumCPU(config.getInt("fit_settings.NumCPU"));

    fitter.setTimeVar(config.getString("observables.decaytime"), config.getDouble("observables.dt_low"), config.getDouble("observables.dt_high"));
    fitter.setFinalstateVar(config.getString("observables.finalstate"));
    if(config.getBool("observables.single_tag")){
        fitter.setSingleMistagVar(config.getString("observables.single_mistag"),
                                  config.getDouble("observables.single_mistag_low"),
                                  config.getDouble("observables.single_mistag_high"));
        fitter.setSingleTagVar(config.getString("observables.single_tag_obs"));
    }
    else{
        fitter.setTagOSVar(config.getString("observables.tag_OS"));
        fitter.setMistagOSVar(config.getString("observables.mistag_OS"),
                              config.getDouble("observables.mistag_OS_low"),
                              config.getDouble("observables.mistag_OS_high"));
        fitter.setTagSSVar(config.getString("observables.tag_SS"));
        fitter.setMistagSSVar(config.getString("observables.mistag_SS"),
                              config.getDouble("observables.mistag_SS_low"),
                              config.getDouble("observables.mistag_SS_high"));
    }

    fitter.use_acceptance(config.getBool("PDF.use_acceptance"), config.GetVector<double>("PDF.acceptance_knots"), config.getBool("PDF.B2DX_acceptance"));
    fitter.UseGLMCalibration(config.getBool("fit_settings.UseGLMCalibration"), config.getBool("constraints.constraintsWithGLM"));
    if(config.getBool("fit_settings.UseGLMCalibration")){
        fitter.setXMLOS(config.getString("fit_settings.XML_OS"));
        if(!config.getBool("observables.single_tag")) fitter.setXMLSS(config.getString("fit_settings.XML_SS"));
    }

    fitter.DoContourPlot(config.getBool("fit_settings.DoContourplot"));
    // simple constraints
    std::vector<std::string> constraints = config.getVoStrings("constraints.1Dconstraints");
    fitter.SetConstraints(constraints);

    // Multidimensional constraints
    std::vector<std::vector<std::string> > MultiVarconstraint;

    //OS Tagging
    std::vector<std::string> OS_constraint = config.getVoStrings("constraints.OS_Constraint");
    // MultiVarconstraint.push_back(OS_constraint);

    //SS Tagging
    std::vector<std::string> SS_constraint = config.getVoStrings("constraints.SS_Constraint");
    // MultiVarconstraint.push_back(SS_constraint);

    fitter.Set4DConstraints(MultiVarconstraint);

    if(config.getBool("data.DataSetInWorkspace")){
        if(data_file == "") data_file = config.getString("data.workspace.file");
        fitter.setData(data_file,
                       config.getString("data.workspace.workspace"),
                       config.getString("data.workspace.dataset"),
                       config.getString("data.cut"),
                       config.getString("data.weight"),
                       config.getBool("observables.truetoytagger"));
    }
    else{
        if(data_file == "") data_file = config.getString("data.ttree.file");
        fitter.setData(data_file,
                       config.getString("data.ttree.tree"),
                       config.getString("data.cut"),
                       config.getString("data.weight"),
                       config.getBool("observables.truetoytagger"));
    }
    fitter.build_PDF();
    fitter.Fit(argc, argv, config.getBool("plot.plot_only"), config.getBool("comparing_fits.do_two_fits"));
    fitter.Plot(config.getBool("plot.plot_decaytime"),
                config.getBool("plot.plot_correlation"),
                config.getBool("plot.plot_oscillation"),
                config.getBool("plot.oscillation_both_finalstates"));

    if(config.getBool("comparing_fits.do_two_fits")){
        fitter.setConfigDir("/home/abirnkraut/tank/CPTimeFit_BD2Dpi/config_SSTagging/");
        fitter.setSingleTagVar("obsTagSS");
        fitter.UseGLMCalibration(false);
        fitter.setData(config.getString("data.ttree.file"),
                       config.getString("data.ttree.tree"),
                       config.getString("data.cut"),
                       config.getString("data.weight"),
                       config.getBool("observables.truetoytagger"));
        fitter.build_PDF();
        fitter.Fit(argc, argv, true);
    }
}
