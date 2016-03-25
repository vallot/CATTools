#ifndef sampleHelpers_h
#define sampleHelpers_h

#include <vector>
#include <string>

#include <TString.h>







/// Namespace to treat different analysis eras as enum types
namespace Era{
    
    /// All analysis eras as needed
    enum Era{run1_8tev, run2_13tev_50ns, run2_13tev_25ns, undefined};
    
    /// Convert an era from string to enum
    Era convert(const TString& era);
    
    /// Convert an era from enum to string
    TString convert(const Era& era);
    
    /// Return energy for given era in TeV
    double energyInTev(const Era era);
}







/// Namespace to define enums needed for b-tag related stuff
/// Also needed for b-tag scale factors, since ROOT dictionary does not want definitions done in ZTopUtils/
namespace Btag{
    
    /// Enum for the implemented b-tagging algorithms
    enum Algorithm{
        csv,
        csvv2,
        csvv2_50ns,
        undefinedAlgorithm
    };
    
    /// Enum for the working points
    enum WorkingPoint{
        L, M, T, undefinedWP
    };
    
    /// Enum for the implemented modes of btag corrections
    enum CorrectionMode{
        noCorrection,               // Do not apply any corrections, i.e. scale factors event SF=1
        greaterEqualOneTagReweight, // Correct selection efficiency for given working point via event SF for >=1 b-tag
        randomNumberRetag,          // Random-number based tag flipping for b-/c-/l-jets to correct for selection efficiency
        discriminatorReweight,      // Reweight with event-wise SF to describe b-tag discriminator distribution
        undefinedCorrectionMode     // Undefined
    };
    
    
    
    /// Convert an Algorithm from string to enum
    Algorithm convertAlgorithm(const std::string& algo);
    
    /// Convert an Algorithm from enum to string
    std::string convertAlgorithm(const Algorithm& algo);
    
    /// Convert a WorkingPoint from string to enum
    WorkingPoint convertWorkingPoint(const std::string& wp);
    
    /// Convert a WorkingPoint from enum to string
    std::string convertWorkingPoint(const WorkingPoint& wp);
    
    /// Convert a CorrectionMode from string to enum
    CorrectionMode convertCorrectionMode(const std::string& mode);
    
    /// Convert a CorrectionMode from enum to string
    std::string convertCorrectionMode(const CorrectionMode& mode);
}







/// Namespace to treat systematics as enum types
namespace Systematic{
    
    /// All systematic types as needed in any part of the framework
    enum Type{
        nominal,            // nominal, i.e. no systematic variation applied
        mH110,              // Higgs mass of 110 GeV
        mH115,              // Higgs mass of 115 GeV
        mH120,              // Higgs mass of 120 GeV
        mH1225,             // Higgs mass of 122.5 GeV
        mH1275,             // Higgs mass of 127.5 GeV
        mH130,              // Higgs mass of 130 GeV
        mH135,              // Higgs mass of 135 GeV
        mH140,              // Higgs mass of 140 GeV
        lept,               // scale lepton ID/ISO data-to-MC scale factors
        trig,               // scale trigger data-to-MC scale factors
        pu,                 // scale pileup data-to-MC scale factors
        dy,                 // uncertainty on the Drell-Yan same-flavour background
        bg,                 // general background uncertainty
        kin,                // scale kinematic reconstruction scale factors
        btag,               // scale b-tagging data-to-MC scale factors of the b-/c-jets
        btagPt,             // median method: scale b-tagging data-to-MC scale factors of the b-/c-jets below/above median pt down/up or up/down
        btagEta,            // median method: scale b-tagging data-to-MC scale factors of the b-/c-jets below/above median eta down/up or up/down
        btagLjet,           // scale b-tagging data-to-MC scale factors of the l-jets
        btagLjetPt,         // median method: scale b-tagging data-to-MC scale factors of the l-jets below/above median pt down/up or up/down
        btagLjetEta,        // median method: scale b-tagging data-to-MC scale factors of the l-jets below/above median eta down/up or up/down
        btagBeff,           // scale the b-tagging efficiencies as estimated from MC for b-jets for stat. uncertainty (not applied anywhere, should it be removed?)
        btagCeff,           // scale the b-tagging efficiencies as estimated from MC for c-jets for stat. uncertainty (not applied anywhere, should it be removed?)
        btagLeff,           // scale the b-tagging efficiencies as estimated from MC for l-jets for stat. uncertainty (not applied anywhere, should it be removed?)
        btagDiscrBpurity,   // for b-tag discriminator reweighting: purity of the HF sample used for the LF SF determination
        btagDiscrLpurity,   // for b-tag discriminator reweighting: purity of the LF sample used for the HF SF determination
        btagDiscrBstat1,    // for b-tag discriminator reweighting: scale part 1 of the statistical uncertainty for b-jets
        btagDiscrBstat2,    // for b-tag discriminator reweighting: scale part 2 of the statistical uncertainty for b-jets
        btagDiscrLstat1,    // for b-tag discriminator reweighting: scale part 1 of the statistical uncertainty for l-jets
        btagDiscrLstat2,    // for b-tag discriminator reweighting: scale part 2 of the statistical uncertainty for l-jets
        btagDiscrCerr1,     // for b-tag discriminator reweighting: scale part 1 of the total uncertainty for c-jets
        btagDiscrCerr2,     // for b-tag discriminator reweighting: scale part 2 of the total uncertainty for c-jets
        jer,                // scale jet energy resolution scale factors
        jes,                // scale jet energy scale scale factors
        frac_tthf,          // correction factor for the fraction of tt+HF events from the template fit
        frac_ttother,       // correction factor for the fraction of tt+Other events from the template fit
        lumi,               // luminosity uncertainty
        xsec_tt2b,          // cross-section uncertainty of tt2b process
        xsec_ttcc,          // cross-section uncertainty of ttcc process
        xsec_ttother,       // cross-section uncertainty of tt+light jets process
        xsec_ttZ,           // cross-section uncertainty of ttZ process
        xsec_ttH,           // cross-section uncertainty of ttH process
        topPt,              // scale top pt as estimated in ttbar differential cross-section measurements
        mass,               // variations of masses used in process generation (here top quark mass)
        match,              // matching uncertainty in process generation
        scale,              // scale uncertainty in process generation
        powheg,             // POWHEG event generator matched to PYTHIA shower
        powhegHerwig,       // POWHEG event generator matched to HERWIG shower
        mcatnlo,            // MC@NLO event generator
        perugia11,          // Perugia11 parton shower tune
        perugia11NoCR,      // Perugia11 parton shower tune, no colour-reconnection
        pdf,                // PDF variations
        closure,            // Closure test
        allAvailable,       // All systematics which are available
        all,                // All allowed systematics
        undefinedType       // No systematic defined (also not nominal)
    };
    
    
    
    /// Convert a type from string to enum
    Type convertType(const TString& type);
    
    /// Convert a type from enum to string
    TString convertType(const Type& type);
    
    /// Convert a vector of types from string to enum
    std::vector<Type> convertType(const std::vector<TString>& types);
    
    /// Convert a vector of types from string to enum
    std::vector<Type> convertType(const std::vector<std::string>& types);
    
    /// Convert a vector of types from string to enum
    std::vector<TString> convertType(const std::vector<Type>& types);
    
    
    
    
    /// All variations as needed in any part of the framework
    enum Variation{up, down, central, undefinedVariation};
    
    
    
    /// Convert a variation from string to enum
    Variation convertVariation(const TString& variation);
    
    /// Convert a variation from enum to string
    TString convertVariation(const Variation& variation);
    
    /// Convert a vector of variations from string to enum
    std::vector<Variation> convertVariation(const std::vector<TString>& variations);
    
    /// Convert a vector of variations from string to enum
    std::vector<Variation> convertVariation(const std::vector<std::string>& variations);
    
    /// Convert a vector of variations from string to enum
    std::vector<TString> convertVariation(const std::vector<Variation>& variations);
    
    
    
    
    
    
    
    /// Define for which systematics up/down variations are allowed
    const std::vector<Type> upDownTypes{
        lept, trig, pu,
        dy, bg, kin,
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
        btagBeff, btagCeff, btagLeff,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrBpurity, btagDiscrLpurity,
        btagDiscrCerr1, btagDiscrCerr2,
        jer, jes,
        frac_tthf, frac_ttother,
        lumi,
        xsec_tt2b, xsec_ttcc, xsec_ttother, 
        xsec_ttH, xsec_ttZ,
        topPt,
        mass, match, scale,
        pdf
    };
    
    /// Define for which systematics central variations are allowed
    /// This is also used to identify for which systematics variation numbers can be assigned
    const std::vector<Type> centralTypes{
        pdf
    };
    
    
    
    /// Check the validity of a variation for a given type
    void isValid(const Type& type, const Variation& variation, const int variationNumber =-1);
    
    
    
    
    
    /// Define b-tag systematics, valid for all b-tag corrections
    const std::vector<Type> btagTypes{
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
        btagBeff, btagCeff, btagLeff,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrBpurity, btagDiscrLpurity,
        btagDiscrCerr1, btagDiscrCerr2,
    };
    
    /// Define b-tag systematics, valid for b-tag corrections concerning discriminator reweighting
    const std::vector<Type> btagDiscriminatorReweightTypes{
        btag, btagLjet,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrBpurity, btagDiscrLpurity,
        btagDiscrCerr1, btagDiscrCerr2,
    };
    
    /// Define b-tag systematics, valid for b-tag corrections concerning efficiency
    const std::vector<Type> btagEfficiencyCorrectionTypes{
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
        btagBeff, btagCeff, btagLeff,
    };
    
    /// Define ttbar systematics, i.e. variations of the ttbar sample (e.g. mass or scale variations)
    const std::vector<Type> ttbarTypes{
        topPt,
        mass, match, scale,
        powheg, powhegHerwig, mcatnlo, perugia11, perugia11NoCR,
        pdf,
        closure,
    };
    
    /// Define cross-section uncertainty systematics, which use nominal samples, and change only the scaling
    const std::vector<Type> crossSectionTypes{
        xsec_tt2b, xsec_ttcc, xsec_ttother, 
        xsec_ttH, xsec_ttZ,
    };
    
    /// Define uncertainties due to tt+HF fraction scale factor from the fit, which use nominal samples, and change only the scaling
    const std::vector<Type> tthfFractionTypes{
        frac_tthf, frac_ttother,
    };
    
    /// Define systematics that do not require dedicated root files
    const std::vector<Type> fileIndependentTypes{
        xsec_tt2b, xsec_ttcc, xsec_ttother, 
        xsec_ttH, xsec_ttZ,
        frac_tthf, frac_ttother,
        lumi
    };
    
    
    
    
    
    /// Class for proper handling of systematic
    class Systematic{
        
    public:
        
        Systematic();
        
        Systematic(const Type& type, const Variation& variation, const int variationNumber =-1);
        
        Systematic(const TString& systematicName);
        
        ~Systematic(){}
        
        bool operator<(const Systematic& rhs)const{return this->name() < rhs.name();}
        
        TString name()const;
        
        Type type()const{return type_;}
        
        Variation variation()const{return variation_;}
        
        int variationNumber()const{return variationNumber_;}
        
        
        
    private:
        
        Type type_;
        
        Variation variation_;
        
        int variationNumber_;
    };
    
    
    
    /// Set all systematics from a list of allowed types, using the defined variations
    std::vector<Systematic> allowedSystematicsAnalysis(const std::vector<Type>& allowedTypes);
    
    /// Set up systematics from vector of systematicNames
    std::vector<Systematic> setSystematics(const std::vector<std::string>& systematicNames);
    
    /// Set up systematic for nominal (i.e. no systematic variation)
    Systematic nominalSystematic();
    
    /// Set up undefined systematic
    Systematic undefinedSystematic();
}









/// Namespace to treat decay channels as enum types
namespace Channel{
    
    /// All dileptonic decay channels as needed in any part of the framework
    enum Channel{ee, emu, mumu, combined, tautau, undefined};
    
    
    
    /// All dileptonic decay channels allowed for analysis step
    /// (allow undefined to select all channels if no option is set, i.e. option is empty)
    const std::vector<Channel> allowedChannelsAnalysis
        {ee, emu, mumu, undefined};
    
    /// All dileptonic decay channels allowed for plotting step
    const std::vector<Channel> allowedChannelsPlotting
        {ee, emu, mumu, combined};
    
    /// Real analysis channels, i.e. all channels which describe a real final state
    const std::vector<Channel> realChannels
        {ee, emu, mumu};
    
    /// Possible Drell-Yan decay channels
    const std::vector<Channel> dyDecayChannels
        {ee, mumu, tautau};
    
    
    
    /// Convert a channel from string to enum
    Channel convert(const TString& channel);
    
    /// Convert a channel from enum to string
    TString convert(const Channel& channel);
    
    /// Return the label of a channel as used for drawing
    TString label(const Channel& channel);
    
    /// Convert a vector of channels from string to enum
    std::vector<Channel> convert(const std::vector<TString>& channels);
    
    /// Convert a vector of channels from string to enum
    std::vector<Channel> convert(const std::vector<std::string>& channels);
    
    /// Convert a vector of channels from string to enum
    std::vector<TString> convert(const std::vector<Channel>& channels);
}







namespace common{
    
    /// Create and assign an output folder depending on the channel and systematic
    TString assignFolder(const char* baseDir, const Channel::Channel& channel, const Systematic::Systematic& systematic, const char* subDir ="");
    
    /// Access an already existing input folder
    TString accessFolder(const char* baseDir, const Channel::Channel& channel,
                         const Systematic::Systematic& systematic, const bool allowNonexisting =false);
    
    /// Access the real final state from a filename, ie. only "ee", "emu", "mumu", but not "combined"
    Channel::Channel finalState(const TString& filename);
    
    /// Find file list for a given channel and systematic, and return its name
    /// In case it does not exist, return empty string
    TString findFilelist(const TString& filelistDirectory,
                         const Channel::Channel& channel,
                         const Systematic::Systematic& systematic);
    
    /// Find from vector of given systematics those for which a file list exists for all given channels
    std::vector<Systematic::Systematic> findSystematicsFromFilelists(const TString& filelistDirectory,
                                                                     const std::vector<Channel::Channel>& v_channel,
                                                                     const std::vector<Systematic::Systematic>& v_systematic);
    
    /// Read the file list for given channel and systematic, and return the input file names
    /// In case a vector of patterns is specified, only files containing this pattern in the full path name will be read
    std::vector<TString> readFilelist(const TString& filelistDirectory,
                                      const Channel::Channel& channel,
                                      const Systematic::Systematic& systematic,
                                      const std::vector<TString>& v_pattern =std::vector<TString>());
    
    /// Read a file for given file name, and return the lines each as element in vector
    /// In case a vector of patterns is specified, only lines containing this pattern will be read
    std::vector<TString> readFile(const TString& filename,
                                  const std::vector<TString>& v_pattern =std::vector<TString>());
}







#endif






