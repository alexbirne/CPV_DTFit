#include "CPTimeFit.h"

void CPTimeFit::build_PDF(){

    RooRealVar *BeautyTime = NULL;
    if(!second_fit_when_comparing_fits_) BeautyTime = (RooRealVar*)ws_.obj(time_var_.c_str());
    else BeautyTime = (RooRealVar*)ws_for_alt_pdf_.obj(time_var_.c_str());
    // BeautyTime->Print();
    if(use_acceptance_) this->build_acceptance(*BeautyTime);

    this->build_resolution(*BeautyTime);

    // now resolution function should be in workspace
    RooResolutionModel* resolution = NULL;
    if(use_acceptance_&&!second_fit_when_comparing_fits_) resolution = (RooGaussEfficiencyModel*)ws_.obj("resolution");
    else if(use_acceptance_&&second_fit_when_comparing_fits_) resolution = (RooGaussEfficiencyModel*)ws_for_alt_pdf_.obj("resolution");
    else if(!use_acceptance_&&!second_fit_when_comparing_fits_) resolution = (RooGaussModel*)ws_.obj("resolution");
    else resolution = (RooGaussModel*)ws_for_alt_pdf_.obj("resolution");

    resolution->Print();

    // creating CP observables
    RooRealVar C_f = RooRealVar("C_f","C_f",0.9995,-5.0,5.0);
    RooFormulaVar C_fbar = RooFormulaVar("C_fbar", "C_fbar", "-1.0*@0", RooArgSet(C_f));
    RooRealVar *S_f = new RooRealVar("S_f","S_{f}",-0.031,-5.0,5.0);
    RooUnblindUniform S_f_blind = RooUnblindUniform("Sf_blind", "S_{f} (blind)", "TD_Dpi_3fb_Sf", 0.5, *S_f);
    RooRealVar *S_fbar = new RooRealVar("S_fbar","S_{#bar f}",0.029,-5.0,5.0);
    RooUnblindUniform S_fbar_blind = RooUnblindUniform("Sfbar_blind", "S_{fbar} (blind)", "TD_Dpi_3fb_Sfbar", 0.5, *S_fbar);
    RooRealVar sinh_f = RooRealVar("sinh_f","D_{f}",0.0,-5.0,5.0);
    RooRealVar sinh_fbar = RooRealVar("sinh_fbar","D_{#bar f}",0.0,-5.0,5.0);
    RooRealVar cosh_f = RooRealVar("cosh_f","cosh_f",1.0,-5.0,5.0);
    RooRealVar cosh_fbar = RooRealVar("cosh_fbar","cosh_fbar",1.0,-5.0,5.0);

    // creating calibration parameters for the Flavour Tagging
    RooRealVar p0 = RooRealVar("p0","p_{0}^{OS}",0.37);
    RooRealVar p1 = RooRealVar("p1","p_{1}^{OS}",1.0);
    RooRealVar Deltap0 = RooRealVar("Deltap0","#Delta p_{0}^{OS}",0.0);
    RooRealVar Deltap1 = RooRealVar("Deltap1","#Delta p_{1}^{OS}",0.0);
    RooRealVar MeanEta = RooRealVar("MeanEta","MeanEta",0.37);

    RooRealVar p0SS = RooRealVar("p0SS","p_{0}^{SS}",0.46);
    RooRealVar p1SS = RooRealVar("p1SS","p_{1}^{SS}",1.0);
    RooRealVar Deltap0SS = RooRealVar("Deltap0SS","#Delta p_{0}^{SS}",0.0);
    RooRealVar Deltap1SS = RooRealVar("Deltap1SS","#Delta p_{1}^{SS}",0.0);
    RooRealVar MeanEtaSS = RooRealVar("MeanEtaSS","MeanEtaSS",0.46);

    // creating all further asymmetries which are implemented in the DecRateCoeff_Bd class at the moment: production asymmetry, detection asymmetry
    RooRealVar production_asym = RooRealVar("production_asym","a_{prod}",0);
    RooRealVar detection_asym = RooRealVar("detection_asym","a_{det}",0);
    // creating tagging efficiencies and their asymmetries
    RooRealVar tageff_os = RooRealVar("tageff_os", "#varepsilon_{tag}^{OS}", 0.37, 0.0, 1.0);
    RooRealVar tageff_asym_os = RooRealVar("tageff_asym_os","tageff_asym_os",0.0);
    RooRealVar tageff_ss = RooRealVar("tageff_ss", "#varepsilon_{tag}^{SS}", 0.71, 0.0, 1.0);
    RooRealVar tageff_asym_ss = RooRealVar("tageff_asym_ss","tageff_asym_ss",0.0);

    RooRealVar *MistagOS = NULL;
    RooCategory *TagDecOS = NULL;
    RooCategory *TagDecSS = NULL;
    RooRealVar *MistagSS = NULL;
    RooCategory *BacCharge = NULL;

    if(!second_fit_when_comparing_fits_){
      if(use_single_tag_) {
        TagDecOS = (RooCategory*)ws_.obj(single_tag_var_.c_str());
        MistagOS = (RooRealVar*)ws_.obj(single_mistag_var_.c_str());
      }
      else{
        TagDecOS = (RooCategory*)ws_.obj(os_tag_var_.c_str());
        MistagOS = (RooRealVar*)ws_.obj(os_mistag_var_.c_str());
        TagDecSS = (RooCategory*)ws_.obj(ss_tag_var_.c_str());
        MistagSS = (RooRealVar*)ws_.obj(ss_mistag_var_.c_str());
      }
      BacCharge = (RooCategory*)ws_.obj(finalstate_var_.c_str());
    }
    else{
      if(use_single_tag_){
        TagDecOS = (RooCategory*)ws_for_alt_pdf_.obj(single_tag_var_.c_str());
        MistagOS = (RooRealVar*)ws_for_alt_pdf_.obj(single_mistag_var_.c_str());
      }
      else{
        TagDecOS = (RooCategory*)ws_for_alt_pdf_.obj(os_tag_var_.c_str());
        MistagOS = (RooRealVar*)ws_for_alt_pdf_.obj(os_mistag_var_.c_str());
        TagDecSS = (RooCategory*)ws_for_alt_pdf_.obj(ss_tag_var_.c_str());
        MistagSS = (RooRealVar*)ws_for_alt_pdf_.obj(ss_mistag_var_.c_str());
      }
      BacCharge = (RooCategory*)ws_for_alt_pdf_.obj(finalstate_var_.c_str());
    }
    // here the coefficients are cerated
    CPfitter::DecRateCoeff_Bd *coeff_c_fit = NULL;
    CPfitter::DecRateCoeff_Bd *coeff_s_fit = NULL;
    CPfitter::DecRateCoeff_Bd *coeff_sh_fit = NULL;
    CPfitter::DecRateCoeff_Bd *coeff_ch_fit = NULL;

    if(use_single_tag_){
        if(use_GlM_){
          glm_os_ = new Espresso::GLMBuilder("Calibration_for_the_OS",
                                                          "OS_Calibration",
                                                          *MistagOS,
                                                          "Calibration_not_used_name",
                                                          xml_file_os_);

          RooRealVar* calib_coeff = NULL;
          RooRealVar* calib_deltacoeff = NULL;
          for(int c = 0; c< glm_os_->coefficients().getSize(); c++){
            calib_coeff = (RooRealVar*)glm_os_->coefficients().at(c);
            calib_coeff->setConstant(true);
            calib_deltacoeff = (RooRealVar*)glm_os_->delta_coefficients().at(c);
            calib_deltacoeff->setConstant(true);
            if(constrained_tagging_){
              calib_coeff->setConstant(false);
              calib_deltacoeff->setConstant(false);
            }
          }

          // in case of constrains adding a vector of strings to the class variable MultiVarconstraints_
          if(constrained_tagging_){
            RooArgList params = glm_os_->coefficients();
            params.add(glm_os_->delta_coefficients());
            std::vector<std::string> MultiVarconstraint;
            for(int i = 0; i < params.getSize(); i++){
              MultiVarconstraint.push_back(params.at(i)->GetName());
            }
            MultiVarconstraints_.push_back(MultiVarconstraint);
          }

          Espresso::RooGLMFunction omega_b = glm_os_->b_mistag();
          Espresso::RooGLMFunction omega_bbar = glm_os_->bbar_mistag();

          coeff_c_fit = new CPfitter::DecRateCoeff_Bd("coef_cos_fit",
                                                      "coef_cos_fit",
                                                      CPfitter::DecRateCoeff_Bd::kCos,
                                                      *BacCharge,
                                                      C_f,
                                                      C_fbar,
                                                      *TagDecOS,
                                                      omega_b,
                                                      omega_bbar,
                                                      tageff_os,
                                                      tageff_asym_os,
                                                      production_asym,
                                                      detection_asym);
          if(ResultBlind_) coeff_s_fit = new CPfitter::DecRateCoeff_Bd("coef_sin_fit",
                                                                       "coef_sin_fit",
                                                                       CPfitter::DecRateCoeff_Bd::kSin,
                                                                       *BacCharge,
                                                                       S_f_blind,
                                                                       S_fbar_blind,
                                                                       *TagDecOS,
                                                                       omega_b,
                                                                       omega_bbar,
                                                                       tageff_os,
                                                                       tageff_asym_os,
                                                                       production_asym,
                                                                       detection_asym);
            else coeff_s_fit = new CPfitter::DecRateCoeff_Bd("coef_sin_fit",
                                                             "coef_sin_fit",
                                                             CPfitter::DecRateCoeff_Bd::kSin,
                                                             *BacCharge,
                                                             *S_f,
                                                             *S_fbar,
                                                             *TagDecOS,
                                                             omega_b,
                                                             omega_bbar,
                                                             tageff_os,
                                                             tageff_asym_os,
                                                             production_asym,
                                                             detection_asym);
          coeff_sh_fit = new CPfitter::DecRateCoeff_Bd("coef_sinh_fit",
                                                       "coef_sinh_fit",
                                                       CPfitter::DecRateCoeff_Bd::kSinh,
                                                       *BacCharge,
                                                       sinh_f,
                                                       sinh_fbar,
                                                       *TagDecOS,
                                                       omega_b,
                                                       omega_bbar,
                                                       tageff_os,
                                                       tageff_asym_os,
                                                       production_asym,
                                                       detection_asym);
          coeff_ch_fit = new CPfitter::DecRateCoeff_Bd("coef_cosh_fit",
                                                       "coef_cosh_fit",
                                                       CPfitter::DecRateCoeff_Bd::kCosh,
                                                       *BacCharge,
                                                       cosh_f,
                                                       cosh_fbar,
                                                       *TagDecOS,
                                                       omega_b,
                                                       omega_bbar,
                                                       tageff_os,
                                                       tageff_asym_os,
                                                       production_asym,
                                                       detection_asym);
          if(second_fit_when_comparing_fits_) {
            ws_for_alt_pdf_.import(*coeff_c_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
            ws_for_alt_pdf_.import(*coeff_s_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
            ws_for_alt_pdf_.import(*coeff_sh_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
            ws_for_alt_pdf_.import(*coeff_ch_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
          }
          else{
            ws_.import(*coeff_c_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
            ws_.import(*coeff_s_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
            ws_.import(*coeff_sh_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
            ws_.import(*coeff_ch_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
          }
        }
        else{
          coeff_c_fit = new CPfitter::DecRateCoeff_Bd("coef_cos_fit",
                                                      "coef_cos_fit",
                                                      CPfitter::DecRateCoeff_Bd::kCos,
                                                      *BacCharge,
                                                      C_f,
                                                      C_fbar,
                                                      *TagDecOS,
                                                      *MistagOS,
                                                      p0,
                                                      p1,
                                                      Deltap0,
                                                      Deltap1,
                                                      MeanEta,
                                                      tageff_os,
                                                      tageff_asym_os,
                                                      production_asym,
                                                      detection_asym);
          if(ResultBlind_) coeff_s_fit = new CPfitter::DecRateCoeff_Bd("coef_sin_fit",
                                                                       "coef_sin_fit",
                                                                       CPfitter::DecRateCoeff_Bd::kSin,
                                                                       *BacCharge,
                                                                       S_f_blind,
                                                                       S_fbar_blind,
                                                                       *TagDecOS,
                                                                       *MistagOS,
                                                                       p0,
                                                                       p1,
                                                                       Deltap0,
                                                                       Deltap1,
                                                                       MeanEta,
                                                                       tageff_os,
                                                                       tageff_asym_os,
                                                                       production_asym,
                                                                       detection_asym);
            else coeff_s_fit = new CPfitter::DecRateCoeff_Bd("coef_sin_fit",
                                                             "coef_sin_fit",
                                                             CPfitter::DecRateCoeff_Bd::kSin,
                                                             *BacCharge,
                                                             *S_f,
                                                             *S_fbar,
                                                             *TagDecOS,
                                                             *MistagOS,
                                                             p0,
                                                             p1,
                                                             Deltap0,
                                                             Deltap1,
                                                             MeanEta,
                                                             tageff_os,
                                                             tageff_asym_os,
                                                             production_asym,
                                                             detection_asym);
          coeff_sh_fit = new CPfitter::DecRateCoeff_Bd("coef_sinh_fit",
                                                       "coef_sinh_fit",
                                                       CPfitter::DecRateCoeff_Bd::kSinh,
                                                       *BacCharge,
                                                       sinh_f,
                                                       sinh_fbar,
                                                       *TagDecOS,
                                                       *MistagOS,
                                                       p0,
                                                       p1,
                                                       Deltap0,
                                                       Deltap1,
                                                       MeanEta,
                                                       tageff_os,
                                                       tageff_asym_os,
                                                       production_asym,
                                                       detection_asym);
          coeff_ch_fit = new CPfitter::DecRateCoeff_Bd("coef_cosh_fit",
                                                       "coef_cosh_fit",
                                                       CPfitter::DecRateCoeff_Bd::kCosh,
                                                       *BacCharge,
                                                       cosh_f,
                                                       cosh_fbar,
                                                       *TagDecOS,
                                                       *MistagOS,
                                                       p0,
                                                       p1,
                                                       Deltap0,
                                                       Deltap1,
                                                       MeanEta,
                                                       tageff_os,
                                                       tageff_asym_os,
                                                       production_asym,
                                                       detection_asym);
          if(second_fit_when_comparing_fits_) {
            ws_for_alt_pdf_.import(*coeff_c_fit, RooFit::RecycleConflictNodes());
            ws_for_alt_pdf_.import(*coeff_s_fit, RooFit::RecycleConflictNodes());
            ws_for_alt_pdf_.import(*coeff_sh_fit, RooFit::RecycleConflictNodes());
            ws_for_alt_pdf_.import(*coeff_ch_fit, RooFit::RecycleConflictNodes());
          }
          else{
            ws_.import(*coeff_c_fit, RooFit::RecycleConflictNodes());
            ws_.import(*coeff_s_fit, RooFit::RecycleConflictNodes());
            ws_.import(*coeff_sh_fit, RooFit::RecycleConflictNodes());
            ws_.import(*coeff_ch_fit, RooFit::RecycleConflictNodes());
          }
        }
    }
    else{
      if(use_GlM_){
        // OS tagging
        glm_os_ = new Espresso::GLMBuilder("Calibration_for_the_OS",
                                           "OS_Calibration",
                                           *MistagOS,
                                           "Calibration_not_used_name_OS",
                                           xml_file_os_);

        RooRealVar* calib_coeff = NULL;
        RooRealVar* calib_deltacoeff = NULL;
        for(int c = 0; c< glm_os_->coefficients().getSize(); c++){
          calib_coeff = (RooRealVar*)glm_os_->coefficients().at(c);
          calib_coeff->setConstant(true);
          calib_deltacoeff = (RooRealVar*)glm_os_->delta_coefficients().at(c);
          calib_deltacoeff->setConstant(true);
          if(constrained_tagging_){
            calib_coeff->setConstant(false);
            calib_deltacoeff->setConstant(false);
          }
        }

        // in case of constrains adding a vector of strings to the class variable MultiVarconstraints_
        if(constrained_tagging_){
          RooArgList params = glm_os_->coefficients();
          params.add(glm_os_->delta_coefficients());
          std::vector<std::string> MultiVarconstraint;
          for(int i = 0; i < params.getSize(); i++){
            MultiVarconstraint.push_back(params.at(i)->GetName());
          }
          MultiVarconstraints_.push_back(MultiVarconstraint);
        }

        Espresso::RooGLMFunction omega_b_os = glm_os_->b_mistag();
        Espresso::RooGLMFunction omega_bbar_os = glm_os_->bbar_mistag();
        ws_.import(omega_b_os, RooFit::RecycleConflictNodes(), RooFit::Silence());
        ws_.import(omega_bbar_os, RooFit::RecycleConflictNodes(), RooFit::Silence());
        // SS tagging
        glm_ss_ = new Espresso::GLMBuilder("Calibration_for_the_SS",
                                           "SS_Calibration",
                                           *MistagSS,
                                           "Calibration_not_used_name_SS",
                                           xml_file_ss_);

        for(int c = 0; c< glm_ss_->coefficients().getSize(); c++){
          calib_coeff = (RooRealVar*)glm_ss_->coefficients().at(c);
          calib_coeff->setConstant(true);
          calib_deltacoeff = (RooRealVar*)glm_ss_->delta_coefficients().at(c);
          calib_deltacoeff->setConstant(true);
          if(constrained_tagging_){
            calib_coeff->setConstant(false);
            calib_deltacoeff->setConstant(false);
          }
        }

        // in case of constrains adding a vector of strings to the class variable MultiVarconstraints_
        if(constrained_tagging_){
          RooArgList params = glm_ss_->coefficients();
          params.add(glm_ss_->delta_coefficients());
          std::vector<std::string> MultiVarconstraint;
          for(int i = 0; i < params.getSize(); i++){
            MultiVarconstraint.push_back(params.at(i)->GetName());
          }
          MultiVarconstraints_.push_back(MultiVarconstraint);
        }

        Espresso::RooGLMFunction omega_b_ss = glm_ss_->b_mistag();
        Espresso::RooGLMFunction omega_bbar_ss = glm_ss_->bbar_mistag();
        ws_.import(omega_b_ss, RooFit::RecycleConflictNodes(), RooFit::Silence());
        ws_.import(omega_bbar_ss, RooFit::RecycleConflictNodes(), RooFit::Silence());

        // building coefficients
        coeff_c_fit = new CPfitter::DecRateCoeff_Bd("coef_cos_fit",
                                                    "coef_cos_fit",
                                                    CPfitter::DecRateCoeff_Bd::kCos,
                                                    *BacCharge,
                                                    C_f,
                                                    C_fbar,
                                                    *TagDecOS,
                                                    *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_OS"))),
                                                    *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_OS"))),
                                                    tageff_os,
                                                    tageff_asym_os,
                                                    *TagDecSS,
                                                    *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_SS"))),
                                                    *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_SS"))),
                                                    tageff_ss,
                                                    tageff_asym_ss,
                                                    production_asym,
                                                    detection_asym);
        if(ResultBlind_) coeff_s_fit = new CPfitter::DecRateCoeff_Bd("coef_sin_fit",
                                                                     "coef_sin_fit",
                                                                     CPfitter::DecRateCoeff_Bd::kSin,
                                                                     *BacCharge,
                                                                     S_f_blind,
                                                                     S_fbar_blind,
                                                                     *TagDecOS,
                                                                     *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_OS"))),
                                                                     *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_OS"))),
                                                                     tageff_os,
                                                                     tageff_asym_os,
                                                                     *TagDecSS,
                                                                     *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_SS"))),
                                                                     *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_SS"))),
                                                                     tageff_ss,
                                                                     tageff_asym_ss,
                                                                     production_asym,
                                                                     detection_asym);
          else coeff_s_fit = new CPfitter::DecRateCoeff_Bd("coef_sin_fit",
                                                           "coef_sin_fit",
                                                           CPfitter::DecRateCoeff_Bd::kSin,
                                                           *BacCharge,
                                                           *S_f,
                                                           *S_fbar,
                                                           *TagDecOS,
                                                           *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_OS"))),
                                                           *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_OS"))),
                                                           tageff_os,
                                                           tageff_asym_os,
                                                           *TagDecSS,
                                                           *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_SS"))),
                                                           *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_SS"))),
                                                           tageff_ss,
                                                           tageff_asym_ss,
                                                           production_asym,
                                                           detection_asym);
        coeff_sh_fit = new CPfitter::DecRateCoeff_Bd("coef_sinh_fit",
                                                     "coef_sinh_fit",
                                                     CPfitter::DecRateCoeff_Bd::kSinh,
                                                     *BacCharge,
                                                     sinh_f,
                                                     sinh_fbar,
                                                     *TagDecOS,
                                                     *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_OS"))),
                                                     *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_OS"))),
                                                     tageff_os,
                                                     tageff_asym_os,
                                                     *TagDecSS,
                                                     *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_SS"))),
                                                     *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_SS"))),
                                                     tageff_ss,
                                                     tageff_asym_ss,
                                                     production_asym,
                                                     detection_asym);
        coeff_ch_fit = new CPfitter::DecRateCoeff_Bd("coef_cosh_fit",
                                                     "coef_cosh_fit",
                                                     CPfitter::DecRateCoeff_Bd::kCosh,
                                                     *BacCharge,
                                                     cosh_f,
                                                     cosh_fbar,
                                                     *TagDecOS,
                                                     *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_OS"))),
                                                     *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_OS"))),
                                                     tageff_os,
                                                     tageff_asym_os,
                                                     *TagDecSS,
                                                     *((RooAbsReal*)(ws_.obj("b_Calibration_for_the_SS"))),
                                                     *((RooAbsReal*)(ws_.obj("bbar_Calibration_for_the_SS"))),
                                                     tageff_ss,
                                                     tageff_asym_ss,
                                                     production_asym,
                                                     detection_asym);
      }
      else{
        coeff_c_fit = new CPfitter::DecRateCoeff_Bd("coef_cos_fit",
                                                    "coef_cos_fit",
                                                    CPfitter::DecRateCoeff_Bd::kCos,
                                                    *BacCharge,
                                                    C_f,
                                                    C_fbar,
                                                    *TagDecOS,
                                                    *MistagOS,
                                                    p0,
                                                    p1,
                                                    Deltap0,
                                                    Deltap1,
                                                    MeanEta,
                                                    tageff_os,
                                                    tageff_asym_os,
                                                    *TagDecSS,
                                                    *MistagSS,
                                                    p0SS,
                                                    p1SS,
                                                    Deltap0SS,
                                                    Deltap1SS,
                                                    MeanEtaSS,
                                                    tageff_ss,
                                                    tageff_asym_ss,
                                                    production_asym,
                                                    detection_asym);
        if(ResultBlind_) {
          coeff_s_fit = new CPfitter::DecRateCoeff_Bd("coef_sin_fit",
                                                      "coef_sin_fit",
                                                      CPfitter::DecRateCoeff_Bd::kSin,
                                                      *BacCharge,
                                                      S_f_blind,
                                                      S_fbar_blind,
                                                      *TagDecOS,
                                                      *MistagOS,
                                                      p0,
                                                      p1,
                                                      Deltap0,
                                                      Deltap1,
                                                      MeanEta,
                                                      tageff_os,
                                                      tageff_asym_os,
                                                      *TagDecSS,
                                                      *MistagSS,
                                                      p0SS,
                                                      p1SS,
                                                      Deltap0SS,
                                                      Deltap1SS,
                                                      MeanEtaSS,
                                                      tageff_ss,
                                                      tageff_asym_ss,
                                                      production_asym,
                                                      detection_asym);
        }
        else{
          coeff_s_fit =  new CPfitter::DecRateCoeff_Bd("coef_sin_fit",
                                                       "coef_sin_fit",
                                                       CPfitter::DecRateCoeff_Bd::kSin,
                                                       *BacCharge,
                                                       *S_f,
                                                       *S_fbar, // this line only changed for investigating the bias in deltaM
                                                       *TagDecOS,
                                                       *MistagOS,
                                                       p0,
                                                       p1,
                                                       Deltap0,
                                                       Deltap1,
                                                       MeanEta,
                                                       tageff_os,
                                                       tageff_asym_os,
                                                       *TagDecSS,
                                                       *MistagSS,
                                                       p0SS,
                                                       p1SS,
                                                       Deltap0SS,
                                                       Deltap1SS,
                                                       MeanEtaSS,
                                                       tageff_ss,
                                                       tageff_asym_ss,
                                                       production_asym,
                                                       detection_asym);
        }
        coeff_sh_fit = new CPfitter::DecRateCoeff_Bd("coef_sinh_fit",
                                                     "coef_sinh_fit",
                                                     CPfitter::DecRateCoeff_Bd::kSinh,
                                                     *BacCharge,
                                                     sinh_f,
                                                     sinh_fbar,
                                                     *TagDecOS,
                                                     *MistagOS,
                                                     p0,
                                                     p1,
                                                     Deltap0,
                                                     Deltap1,
                                                     MeanEta,
                                                     tageff_os,
                                                     tageff_asym_os,
                                                     *TagDecSS,
                                                     *MistagSS,
                                                     p0SS,
                                                     p1SS,
                                                     Deltap0SS,
                                                     Deltap1SS,
                                                     MeanEtaSS,
                                                     tageff_ss,
                                                     tageff_asym_ss,
                                                     production_asym,
                                                     detection_asym);
        coeff_ch_fit = new CPfitter::DecRateCoeff_Bd("coef_cosh_fit",
                                                     "coef_cosh_fit",
                                                     CPfitter::DecRateCoeff_Bd::kCosh,
                                                     *BacCharge,
                                                     cosh_f,
                                                     cosh_fbar,
                                                     *TagDecOS,
                                                     *MistagOS,
                                                     p0,
                                                     p1,
                                                     Deltap0,
                                                     Deltap1,
                                                     MeanEta,
                                                     tageff_os,
                                                     tageff_asym_os,
                                                     *TagDecSS,
                                                     *MistagSS,
                                                     p0SS,
                                                     p1SS,
                                                     Deltap0SS,
                                                     Deltap1SS,
                                                     MeanEtaSS,
                                                     tageff_ss,
                                                     tageff_asym_ss,
                                                     production_asym,
                                                     detection_asym);
      }
      if(second_fit_when_comparing_fits_) {
        ws_for_alt_pdf_.import(*coeff_c_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
        ws_for_alt_pdf_.import(*coeff_s_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
        ws_for_alt_pdf_.import(*coeff_sh_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
        ws_for_alt_pdf_.import(*coeff_ch_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
      }
      else{
        ws_.import(*coeff_c_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
        ws_.import(*coeff_s_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
        ws_.import(*coeff_sh_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
        ws_.import(*coeff_ch_fit, RooFit::RecycleConflictNodes(), RooFit::Silence());
      }
    }

    // RooFormulaVar* coeff_ch_fit = new RooFormulaVar("coef_cosh_fit",
    //                                                 "coef_cosh_fit",
    //                                                 "(@0==+1)*(@1*(1+@3)*(1-@4*@5))+(@0==-1)*(@2*(1-@3)*(1-@4*@5))",
    //                                                 RooArgList(*BacCharge,cosh_f,cosh_fbar,detection_asym,*TagDecOS,production_asym));
    // RooFormulaVar* coeff_sh_fit = new RooFormulaVar("coef_sinh_fit",
    //                                                 "coef_sinh_fit",
    //                                                 "(@0==+1)*(@1*(1+@3)*(1-@4*@5))+(@0==-1)*(@2*(1-@3)*(1-@4*@5))",
    //                                                 RooArgList(*BacCharge,sinh_f,sinh_fbar,detection_asym,*TagDecOS,production_asym));
    // RooFormulaVar* coeff_s_fit = new RooFormulaVar("coef_sin_fit",
    //                                                "coef_sin_fit",
    //                                                "(@0==+1)*-1*@1*(1+@4)*(@3-@5)+(@0==-1)*-1*@2*(1-@4)*(@3-@5)",
    //                                                RooArgList(*BacCharge,*S_f,*S_fbar,*TagDecOS,detection_asym,production_asym));
    // RooFormulaVar* coeff_c_fit = new RooFormulaVar("coef_cos_fit",
    //                                                "coef_cos_fit",
    //                                                "(@0==+1)*@1*@3+(@0==-1)*@2*@3",
    //                                                RooArgList(*BacCharge,C_f,C_fbar,*TagDecOS,detection_asym,production_asym));

    // Bringing all together in RooBDecay (plus creating RooRealVars for the lifetime tau, delta Gamma and delta M)
    RooRealVar tau = RooRealVar("tau","#tau",1.519,1.0,2.0);
    tau.setConstant(true);
    RooRealVar dgamma = RooRealVar("dgamma","dgamma",0.0);
    dgamma.setConstant(true);
    RooRealVar deltaM = RooRealVar("deltaM","#Delta m",0.51);
    deltaM.setConstant(true);
    RooBDecay *pdfSigTime_noExtension_fit = NULL;
    if(second_fit_when_comparing_fits_){
      pdfSigTime_noExtension_fit = new RooBDecay("pdfSigTime_noExtension_fit",
                                                 "pdfSigTime_noExtension_fit",
                                                 *BeautyTime,
                                                 tau,
                                                 dgamma,
                                                 *((RooAbsReal*)(ws_for_alt_pdf_.obj("coef_cosh_fit"))),
                                                 *((RooAbsReal*)(ws_for_alt_pdf_.obj("coef_sinh_fit"))),
                                                 *((RooAbsReal*)(ws_for_alt_pdf_.obj("coef_cos_fit"))),
                                                 *((RooAbsReal*)(ws_for_alt_pdf_.obj("coef_sin_fit"))),
                                                 deltaM,
                                                 *resolution,
                                                 RooBDecay::SingleSided);
    }
    else{
      pdfSigTime_noExtension_fit = new RooBDecay("pdfSigTime_noExtension_fit",
                                                 "pdfSigTime_noExtension_fit",
                                                 *BeautyTime,
                                                 tau,
                                                 dgamma,
                                                 *((RooAbsReal*)(ws_.obj("coef_cosh_fit"))),
                                                 *((RooAbsReal*)(ws_.obj("coef_sinh_fit"))),
                                                 *((RooAbsReal*)(ws_.obj("coef_cos_fit"))),
                                                 *((RooAbsReal*)(ws_.obj("coef_sin_fit"))),
                                                 deltaM,
                                                 *resolution,
                                                 RooBDecay::SingleSided);
    }

    this->build_mistaghistogram("OS");
    if(!use_single_tag_) this->build_mistaghistogram("SS");

    RooAbsPdf *mistagPDF_OS = ws_.pdf("MistagPDF_OS");
    if(second_fit_when_comparing_fits_) mistagPDF_OS = ws_for_alt_pdf_.pdf("MistagPDF_OS");
    RooAbsPdf *mistagPDF_SS = NULL;
    if(!use_single_tag_) mistagPDF_SS = ws_.pdf("MistagPDF_SS");

    // std::cout << "Print mistag pdfs" << std::endl;
    // mistagPDF_OS->Print("V");
    // if(!use_single_tag_) mistagPDF_SS->Print("V");

    // building conditional PDF
    RooArgSet comp_dimensions = RooArgSet(*pdfSigTime_noExtension_fit);
    RooArgSet* comp_conditional_dimensions = NULL;
    RooArgSet* comp_conditional_observables = NULL;
    if(use_single_tag_){
        comp_conditional_dimensions = new RooArgSet(*mistagPDF_OS);
        comp_conditional_observables = new RooArgSet(*BeautyTime,*TagDecOS,*BacCharge);
    }
    else{
        comp_conditional_dimensions = new RooArgSet(*mistagPDF_OS,*mistagPDF_SS);
        comp_conditional_observables = new RooArgSet(*BeautyTime,*TagDecOS,*TagDecSS,*BacCharge);
    }
    RooProdPdf pdfSigTime_conditional =  RooProdPdf("pdfSigTime_conditional",
                                                    "pdfSigTime_conditional",
                                                    *comp_conditional_dimensions,
                                                    RooFit::Conditional(comp_dimensions, *comp_conditional_observables));

    // RooRealVar SigYield = RooRealVar("SigYield","SigYield",154010,0,1000000);
    // SigYield.setConstant(false);
    // RooExtendPdf pdfSigTime_fit = RooExtendPdf("pdfSigTime_fit","pdfSigTime_fit",*pdfSigTime_conditional,SigYield);

    // RooCategory *MixState = (RooCategory*)ws_.obj("catMixingTrue");
    // RooRealVar omega = RooRealVar("omega", "omega", 0);
    // RooRealVar deltaOmega = RooRealVar("deltaOmega", "deltaOmega", 0);
    // RooBMixDecay pdfSigTime_conditional = RooBMixDecay("pdfSigTime_conditional",8
    //                                                    "pdfSigTime_conditional",
    //                                                    *BeautyTime,
    //                                                    *MixState,
    //                                                    *TagDecOS,
    //                                                    tau,
    //                                                    deltaM,
    //                                                    omega,
    //                                                    deltaOmega,
    //                                                    *resolution,
    //                                                    RooBMixDecay::SingleSided);
    if(second_fit_when_comparing_fits_) ws_for_alt_pdf_.import(pdfSigTime_conditional, RooFit::RecycleConflictNodes(), RooFit::Silence());
    else ws_.import(pdfSigTime_conditional, RooFit::RecycleConflictNodes(), RooFit::Silence());
}

