/************************************************************************************/
/* Largely inspired by the similar class in RAT, written by P.G.Jones and M.Mottram */
/************************************************************************************/

// Look at refactoring this 
#ifndef __OXSX_MINUIT__
#define __OXSX_MINUIT__
#include <Optimiser.h>
#include <string>
#include <vector>
#include <MinuitFCN.h>
#include <Minuit2/MnApplication.h>
#include <FitResult.h>

class TestStatistic;

class Minuit : public Optimiser{
 public:
    Minuit() :  fMethod("Migrad"),
                fMinimiser(NULL), fMaxCalls(0), 
                fTolerance(0.1), fMaximising(false) {}
    ~Minuit();

    virtual const FitResult& Optimise(TestStatistic*);

    void Fix(size_t index_);        // these three methods only work after a call to Initialise
    void Release(size_t index_);
    void RemoveLimits(size_t index_);

    void SetMethod(const std::string&);
    std::string GetMethod() const;

    void SetInitialValues(const std::vector<double>&);
    void SetInitialErrors(const std::vector<double>&);

    void   SetUpperContourEdge(double);
    double GetUpperContourEdge() const;

    void SetMinima(const std::vector<double>& minima_);
    std::vector<double> GetMinima() const;

    void SetMaxima(const std::vector<double>& maxima_);
    std::vector<double> GetMaxima() const;
   
    void     SetMaxCalls(unsigned);
    unsigned GetMaxCalls() const;
    
    void   SetTolerance(double);
    double GetTolerance() const;

    void SetMaximising(bool b_) {fMaximising = b_;}
    bool GetMaximising() const  {return fMaximising;}

    FitResult GetFitResult() const;
    
 private:
	void Initialise();
    MinuitFCN   fMinuitFCN; // wrapper on evaluator so migrad can call it
    std::vector<double> fInitialValues;
    std::vector<double> fInitialErrors;

    std::vector<double> fMinima;
    std::vector<double> fMaxima;

    unsigned fMaxCalls;
    double   fTolerance;

    std::string fMethod;
    ROOT::Minuit2::MnApplication* fMinimiser;

    FitResult fFitResult;
    bool fMaximising;
};
#endif
