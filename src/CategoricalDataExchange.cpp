#include "CategoricalDataExchange.h"
#include "coclust/src/Models/CategoricalLBModel.h"

void CategoricalDataExchange::Output(Rcpp::S4& obj,ICoClustModel* model,bool successful){
  if(!successful)
  {
    obj.slot("successful") = false;
    obj.slot("message") = model->GetErrormsg();
  }
  else
  {
    obj.slot("successful") = true;
    obj.slot("message") = "Co-Clustering successfully terminated!";
    CategoricalLBModel* ptrLBM;
    MatrixReal dispersion;

    ptrLBM = dynamic_cast<CategoricalLBModel*>(model);
    const std::vector<MatrixReal> mean = ptrLBM->Getmean();
    std::vector<std::vector<std::vector<double> > > tempmean(mean.size(),std::vector<std::vector<double> >(Mparam_.nbrowclust_,(std::vector<double>(Mparam_.nbcolclust_))));
    //std::vector<Rcpp::NumericMatrix> tempmean(mean.size());
    for (int h = 0; h < mean.size(); ++h) {
      for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
        for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
         tempmean[h][k][l] = mean[h](k,l);
        }
      }
    }

    obj.slot("classmean") = Rcpp::wrap(tempmean);
    obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBM->GetArrangedDataClusters());
    obj.slot("rowclass") = convertvector<Rcpp::IntegerVector,VectorInteger>(model->GetRowClassificationVector());
    obj.slot("colclass") = convertvector<Rcpp::IntegerVector,VectorInteger>(model->GetColumnClassificationVector());
    obj.slot("rowproportions") = convertvector<Rcpp::NumericVector,VectorReal>(model->GetRowProportions());
    obj.slot("columnproportions") = convertvector<Rcpp::NumericVector,VectorReal>(model->GetColProportions());
    obj.slot("rowposteriorprob") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(model->GetRowPosteriorprob());
    obj.slot("colposteriorprob") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(model->GetColPosteriorprob());
    obj.slot("likelihood") = model->GetLikelihood();
    obj.slot("ICLvalue") = model->ICLCriteriaValue();
  }
}

void CategoricalDataExchange::DataInput(Rcpp::S4& obj){
  Rcpp::NumericMatrix data(SEXP(obj.slot("data")));
  convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(data,m_Dataij_);
  Mparam_.nbrowdata_ = m_Dataij_.rows();
  Mparam_.nbcoldata_ = m_Dataij_.cols();

  //Get Strategy
  Rcpp::S4 strategy(obj.slot("strategy"));
  //get hyper-parameters
  if(Rcpp::as<bool>(strategy.slot("bayesianform")))
  {
    Rcpp::IntegerVector hyperparam(SEXP(strategy.slot("hyperparam")));
    a_ = hyperparam(0);
    b_ = hyperparam(1);
  }
}

void CategoricalDataExchange::instantiateModel(ICoClustModel*& model){
  if(!strategy_.Bayesianform_){
    if(!strategy_.SemiSupervised){
      switch (strategy_.Model_)
      {
        case pi_rho_multi:
          Mparam_.fixedproportions_ = true;
          model = new CategoricalLBModel(m_Dataij_,Mparam_);
          break;
        case pik_rhol_multi:
          Mparam_.fixedproportions_ = false;
          model = new CategoricalLBModel(m_Dataij_,Mparam_);
          break;
        default:
          break;
      }
    }else{
      switch (strategy_.Model_){
        case pi_rho_multi:
          Mparam_.fixedproportions_ = true;
          model = new CategoricalLBModel(m_Dataij_,v_rowlabels_,v_collabels_,Mparam_);
          break;
        case pik_rhol_multi:
          Mparam_.fixedproportions_ = false;
          model = new CategoricalLBModel(m_Dataij_,v_rowlabels_,v_collabels_,Mparam_);
          break;
        default:
          break;
      }
    }
  }else{
    if(!strategy_.SemiSupervised){
      switch (strategy_.Model_)
      {
        case pi_rho_multi:
          Mparam_.fixedproportions_ = true;
          model = new CategoricalLBModel(m_Dataij_,Mparam_,a_,b_);
          break;
        case pik_rhol_multi:
          Mparam_.fixedproportions_ = false;
          model = new CategoricalLBModel(m_Dataij_,Mparam_,a_,b_);
          break;
        default:
          break;
      }
    }else{
      switch (strategy_.Model_){
        case pi_rho_multi:
          Mparam_.fixedproportions_ = true;
          model = new CategoricalLBModel(m_Dataij_,v_rowlabels_,v_collabels_,Mparam_,a_,b_);
          break;
        case pik_rhol_multi:
          Mparam_.fixedproportions_ = false;
          model = new CategoricalLBModel(m_Dataij_,v_rowlabels_,v_collabels_,Mparam_,a_,b_);
          break;
        default:
          break;
      }
    }
  }
}
