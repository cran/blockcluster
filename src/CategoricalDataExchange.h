#ifndef CATEGORICALDATAEXCHANGE_H_
#define CATEGORICALDATAEXCHANGE_H_
/**@file CategoricalDataExchange.h
 * @brief 
 */
#include "IDataExchange.h"
#include "coclust/src/Models/CategoricalLBModel.h"
class CategoricalDataExchange: public IDataExchange
{
  public:
    CategoricalDataExchange(){};
    virtual void Output(Rcpp::S4& obj,ICoClustModel*,bool);
    virtual void DataInput(Rcpp::S4 & obj);
    virtual void instantiateModel(ICoClustModel*& model);
    inline const MatrixInteger& GetData() const {return m_Dataij_;}
    virtual ~CategoricalDataExchange(){};
  protected:
    MatrixInteger m_Dataij_;
    int a_,b_;
};

#endif /* CATEGORICALDATAEXCHANGE_H_ */
