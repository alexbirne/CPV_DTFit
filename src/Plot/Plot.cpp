#include "CPTimeFit.h"

void CPTimeFit::Plot(bool plot_decaytime, bool plot_correlation, bool plot_oscillation, bool both_finalstates){

	// first plotting the decay time distribution
	// therefore getting the time observable, data and pdf from the workspace
	RooRealVar *BeautyTime = (RooRealVar*)ws_.obj(time_var_.c_str());
	RooDataSet *data = dynamic_cast<RooDataSet*>(ws_.data("dataset"));
    if(bootstrap_) data = dynamic_cast<RooDataSet*>(ws_.data("bootstrapped_dataset"));
	RooAbsPdf *pdf = ws_.pdf("pdfSigTime_conditional"); //checked that the pdf in the workspace is already the updated one

	// Timeplot
    if(plot_decaytime){
        doofit::plotting::PlotConfig cfg_plot_time("cfg_plot_time");
        cfg_plot_time.InitializeOptions();
        cfg_plot_time.set_plot_appendix("");
        cfg_plot_time.set_plot_directory(plot_dir_+"time");
        cfg_plot_time.set_num_cpu(10);
        std::vector<std::string> components_time;
        doofit::plotting::Plot Time_plot(cfg_plot_time, *BeautyTime, *data, *pdf, components_time);
        Time_plot.set_scaletype_y(doofit::plotting::kBoth);
        Time_plot.PlotIt();
    }

    if(plot_correlation&&!save_toy_with_toystudy_){
        if(did_Fit_) {
            const RooFitResult *fitresult = (const RooFitResult*)ws_.obj("fitresult_pdfSigTime_conditional_dataset");
            doofit::plotting::correlations::CorrelationPlot corr_plot(*fitresult);
            corr_plot.Plot(plot_dir_ + "correlation/");
        }
        else{
            doocore::io::sinfo << "-INFO-  \t" << "Plotting correlation matrix required but no Fit performed." <<  doocore::io::endmsg;
            doocore::io::sinfo << "-INFO-  \t" << "Therefore no correlation matrix plotted!" <<  doocore::io::endmsg;
        }
    }

    if(plot_oscillation){
        if(use_single_tag_){
            double array[52];
            for(int i=0; i<50; i++){
                array[i] = 0.4 + i*11.6/100.0;
            }
            for(int i=0; i<10; i++){
                array[i+50] = array[49] + (i+1)*0.226;
            }
            array[60] = 10.0;
            array[61] = 12.0;
            RooBinning binning(61,array);

            RooCategory *TagDec = (RooCategory*)ws_.obj(single_tag_var_.c_str());
            doocore::io::sinfo << "-INFO-  \t" << "Plotting asymmetry for finalstate f..." <<  doocore::io::endmsg;
            std::string finalstate_cut = finalstate_var_ + "==1";
            RooDataSet *dataset_f = (RooDataSet*)data->reduce(finalstate_cut.c_str());

            doofit::plotting::PlotConfig cfg_plot_asymmetry("cfg_plot_asymmetry");
            cfg_plot_asymmetry.InitializeOptions();
            cfg_plot_asymmetry.set_plot_appendix("_f");
            cfg_plot_asymmetry.set_plot_directory(plot_dir_+"asymmetry_f");
            cfg_plot_asymmetry.set_num_cpu(1);
            std::vector<std::string> components_asymmetry;
            doofit::plotting::Plot Asymmetry_plot(cfg_plot_asymmetry, *BeautyTime, *dataset_f, *pdf, components_asymmetry);
            Asymmetry_plot.AddPlotArg(RooFit::Asymmetry(*TagDec));
            Asymmetry_plot.AddPlotArg(RooFit::ProjWData(*dataset_f));
            Asymmetry_plot.AddPlotArgData(RooFit::Asymmetry(*TagDec));
            Asymmetry_plot.AddPlotArgData(RooFit::Binning(binning));
            Asymmetry_plot.set_plot_asymmetry(true);
            Asymmetry_plot.PlotIt();

            if(both_finalstates){
                doocore::io::sinfo << "-INFO-  \t" << "Plotting asymmetry for finalstate fbar..." <<  doocore::io::endmsg;
                std::string finalstate_cut = finalstate_var_ + "==-1";
                RooDataSet *dataset_fbar = (RooDataSet*)data->reduce(finalstate_cut.c_str());

                doofit::plotting::PlotConfig cfg_plot_asymmetry("cfg_plot_asymmetry");
                cfg_plot_asymmetry.InitializeOptions();
                cfg_plot_asymmetry.set_plot_appendix("_fbar");
                cfg_plot_asymmetry.set_plot_directory(plot_dir_+"asymmetry_fbar");
                cfg_plot_asymmetry.set_num_cpu(1);
                std::vector<std::string> components_asymmetry;
                doofit::plotting::Plot Asymmetry_plot(cfg_plot_asymmetry, *BeautyTime, *dataset_fbar, *pdf, components_asymmetry);
                Asymmetry_plot.AddPlotArg(RooFit::Asymmetry(*TagDec));
                Asymmetry_plot.AddPlotArg(RooFit::ProjWData(*dataset_fbar));
                Asymmetry_plot.AddPlotArgData(RooFit::Asymmetry(*TagDec));
                Asymmetry_plot.AddPlotArgData(RooFit::Binning(binning));
                Asymmetry_plot.set_plot_asymmetry(true);
                Asymmetry_plot.PlotIt();
            }
        }
        else {
            doocore::io::sinfo << "-INFO-  \t" << "Asymmetry plotting not defined for two taggers." <<  doocore::io::endmsg;
            doocore::io::sinfo << "-INFO-  \t" << "Therefore no asymmetry plot created!" <<  doocore::io::endmsg;
        }
    }
}

