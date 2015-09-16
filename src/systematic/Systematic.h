/***************************************************************************************/
/* Responsible for constructing a response matrix and matching systematic/pdf indicies */
/***************************************************************************************/

#ifndef __SYSTEMATIC__
#define __SYSTEMATIC__
#include <PdfMapping.h>
#include <BinnedPdf.h>
#include <DataRepresentation.h>
#include <vector>

class Systematic{
 public:
    Systematic()   {}
    virtual ~Systematic()  {}

    BinnedPdf 
    operator()(const BinnedPdf& pdf_) const;
        
    void 
    SetResponse(const std::vector<std::vector<double> >& responseMatrix_);
    const PdfMapping& GetResponse() const;
        
    void SetDataRep(const DataRepresentation&);
    DataRepresentation GetDataRep() const;

    void SetPdfDataRep(const DataRepresentation&);
    DataRepresentation GetPdfDataRep() const;

    void SetParameters(const std::vector<double>&);
    const std::vector<double>& GetParameters() const;
    size_t GetParamCount() const;

 protected:
    PdfMapping fPdfMapping;
    std::vector<double> fParams;
    DataRepresentation fDataRep;     // the data indicies that this systematic acts on
    DataRepresentation fPdfDataRep;  
    // the data indices  of the pdfs it will act on, needs to be at least the lenth of the 
    // systematics representation
    
    bool BinsCompatible(size_t bin1 , size_t bin2) const;
    bool VectorContains(const std::vector<size_t>&, size_t) const;
};
#endif
