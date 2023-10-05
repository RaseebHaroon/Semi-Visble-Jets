// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "fastjet/tools/Filter.hh"
#include "Rivet/Projections/Thrust.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
namespace Rivet {

class MC_Darkmatter : public Analysis
  {
  public:
  RIVET_DEFAULT_ANALYSIS_CTOR(MC_Darkmatter);

    void init()
    {
        // Projection of final state particles
        const FinalState fs(Cuts::abseta < 4.9);

        // Using anti-kt algorithm for both, ignoring muons and neutrinos
        // Fat Jet (Radius = 1)
        FastJets fj(fs, FastJets::ANTIKT, 1.0, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(fj, "FJets");

        // Small Jet (Radius = 0.4)
        FastJets sj(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(sj, "Jets");

        // TRIMMER
        // Predefined criteria for the trimmed jets
        _trimmer = fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));

        // MET
        // Only used for MET construction using particles
        VetoedFinalState vfs;
        vfs.addVetoPairId(PID::MUON);
        vfs.vetoNeutrinos();
        declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");
        declare(VisibleFinalState(Cuts::abseta < 4.9), "vfs");

        // Histogram booking
        book(_h["dj0_lha"]     , "lha_0"    , 35,   0    ,    1   );
        book(_h["dj0_c2"]      , "c2_0"     , 35,   0    ,    0.6 );
        book(_h["dj0_t21"]     , "tau21_0"  , 35,  -0.05 ,    1.05);
        book(_h["dj0_t32"]     , "tau32_0"  , 35,  -0.05 ,    1.05);
        book(_h["pTFJ"]        , "FpT"      , 35, 300    , 1600   );
        book(_h["dj0_tm"]      , "tmass0"   , 35,   0    ,  500   );
        book(_h["transThrust"] , "thrust"   , 25,  -9    ,    0   );
        book(_h["eTmiss_p"]    , "MET_Pinv" , 40,   0    , 1000   );
        book(_h["eTmiss_J"]    , "MET_J"    , 40,   0    , 1500   );
        book(_h["deltaPhi_F0M"], "F0-M_Dphi", 25,   0    ,    4   );
        book(_h["deltaPhi_F1M"], "F1-M_Dphi", 25,   0    ,    4   );
    }

    void analyze(const Event &event)
    {
        // Using small jets for Thrust and MET calculation only
        const Jets &jet = apply<FastJets>(event, "Jets").jetsByPt(Cuts::abseta < 4.9 && Cuts::pT > 30 * GeV);
        
        // Thrust Calculation
        std::vector<const Jet *> goodJets;
        goodJets.clear();
        std::vector<Vector3> momenta;
        momenta.clear();

        for (const Jet &j : jet) {
            if (j.pt() > 30 * GeV) {
                // using pt critera for goodJets 
                goodJets.push_back(&j);
                Vector3 jet3 = j.p3();
                //Constructing momenta 3 vector for thrust calculations
                momenta.push_back(jet3);
            }
        }

        Thrust thrust;
        thrust.calc(momenta);
        const double transThrust = log(1.0 - thrust.thrust());

        // Using goodJets as a sieve
        if (goodJets.size() > 0) {
            _h["transThrust"]->fill(transThrust);
        }

        // Setting initial values for most variables
        double lha(-99), c2(-99), tau32(-99), tau21(-99);
        lha = 0.0;
        double Dphi0 = -99;
        double Dphi1 = -99;

        const double beta = 1, Rcut = 1;

        // MET calculations using particles
        FourMomentum pTmiss;
        for (const Particle &p : apply<VisibleFinalState>(event, "vfs").particles()) {
            pTmiss -= p.momentum();
        }

        double eTmiss = pTmiss.pT();
        _h["eTmiss_p"]->fill(eTmiss);

        // MET calculations using small jets
        FourMomentum p4missj;
        for (const Jet &j : jet)
            p4missj -= j;
        _h["eTmiss_J"]->fill(p4missj.pT());

        // Introducing Fat Jet
        const Jets &fjets = apply<FastJets>(event, "FJets").jetsByPt(Cuts::abseta < 4.9 && Cuts::pT > 250 * GeV);

        // Criteria for successfull Events
        if (fjets.size() < 1)// atlest one Fat Jets
            vetoEvent;
        if (p4missj.size() < 1)// atlest one MET vector
            vetoEvent;
        if (p4missj.pT() < 200 * GeV)// atleast 200 MeV MET
            vetoEvent;
        if (p4missj.absrap() > 4.9)// |Eta| > 4.9
            vetoEvent;

        // Trimming the Fat Jet
        PseudoJets tr_ljets;
        for (size_t i = 0; i < fjets.size(); ++i) {
            tr_ljets += _trimmer(fjets[i]);
            tr_ljets[tr_ljets.size() - 1].set_user_index(i);
        }

        // Selection criteria for the Trimmed Jets
        ifilter_select(tr_ljets, [](const PseudoJet &j) {
            return j.perp() > 250 * GeV && j.eta() < 4.9 && j.eta() > -4.9;
        });

        // Atlest 1 trimmed jet is needed every event
        if (tr_ljets.size() < 1) 
            vetoEvent;

        // sorting Trimmed jets by pt in a vector
        if (tr_ljets.size() > 0)
            tr_ljets = sorted_by_pt(tr_ljets);

        // Leading Trimmed Fat Jet
        const Jet &fjet0 = tr_ljets[0]; 

        // Seperation b/w Trimmed Leading Jet and MET
        Dphi0 = mapAngle0ToPi(p4missj.phi() - fjet0.phi());
        _h["deltaPhi_F0M"]->fill(Dphi0);

        // First Loop
        ////////////////////////////////////////////////////////////////////////////
        //// Mono-Jet events: Events with only 1 Fat Jet overlapping with MET /////
        //////////////////////////////////////////////////////////////////////////
        if (tr_ljets.size() == 1 && Dphi0 < 1.0) {

            // Trimmed Leading Fat Jet (from mono-jet event)
            PseudoJet LJet = fjet0;

            //  Les Houches Angularity(LHA)
            for (const PseudoJet &p : LJet.constituents()) {
                double trpt = p.pt();
                double trtheta = p.squared_distance(LJet);
                lha += pow(trpt, 1.0) * pow(trtheta, 0.25);
            }
            double lterm = pow(LJet.pt(), 1.0) * pow(1.0, 0.5);
            if (lterm != 0)
                lha /= lterm;
            else
                lha = -99;

            // Energy Correlation Function (C2)
            fastjet::contrib::EnergyCorrelator ECF3(3, beta, fastjet::contrib::EnergyCorrelator::pt_R);
            fastjet::contrib::EnergyCorrelator ECF2(2, beta, fastjet::contrib::EnergyCorrelator::pt_R);
            fastjet::contrib::EnergyCorrelator ECF1(1, beta, fastjet::contrib::EnergyCorrelator::pt_R);
            double recf3 = ECF3(LJet);
            double recf2 = ECF2(LJet);
            double recf1 = ECF1(LJet);
            c2 = (recf2 != 0 ? recf3 * recf1 / (recf2 * recf2) : -1);

            // NSubjettiness
            fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
            fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
            fastjet::contrib::Nsubjettiness nSub3(3, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
            double tau1 = nSub1.result(LJet);
            double tau2 = nSub2.result(LJet);
            double tau3 = nSub3.result(LJet);
            if (tau1 != 0)
                tau21 = tau2 / tau1;
            else
                tau21 = -99;
            if (tau2 != 0)
                tau32 = tau3 / tau2;
            else
                tau32 = -99;

            // Trimmed Jet pT
            const double pt = LJet.pt();

            // Trimmed Jet Mass
            const double mass = LJet.m();

            _h["pTFJ"]->fill(pt);
            _h["dj0_tm"]->fill(mass);
            _h["dj0_c2"]->fill(c2);
            _h["dj0_lha"]->fill(lha);
            _h["dj0_t21"]->fill(tau21);
            _h["dj0_t32"]->fill(tau32);

        ///////////////////////////////////////////////////////////////////////////////////
        // Multi-Jet events: Events with one or both of the 2 Fat Jets overlap with MET///
        /////////////////////////////////////////////////////////////////////////////////
        } else if (tr_ljets.size() > 1) {

            // Sub-Leading Trimmed Fat Jet
            const Jet &fjet1 = tr_ljets[1]; 

            // Seperation b/w Trimmed Sub Jet and MET
            Dphi1 = mapAngle0ToPi(p4missj.phi() - fjet1.phi());
            _h["deltaPhi_F1M"]->fill(Dphi1);

            // First Condition
            // If the Leading trimmed Jet is closer to MET than Sub-Leading Trimmed Jet and both are within 1.0 
            if (Dphi0 < Dphi1 && Dphi0 < 1.0) {
                // Trimmed Leading fatjet (from multijet event)
                PseudoJet LJet = fjet0;

                // Les Houches Angularity(LHA)
                for (const PseudoJet &p : LJet.constituents()) {
                    double trpt = p.pt();
                    double trtheta = p.squared_distance(LJet);
                    lha += pow(trpt, 1.0) * pow(trtheta, 0.25);
                }
                double lterm = pow(LJet.pt(), 1.0) * pow(1.0, 0.5);
                if (lterm != 0)
                    lha /= lterm;
                else
                    lha = -99;

                // Energy Correlation Function (C2)
                fastjet::contrib::EnergyCorrelator ECF3(3, beta, fastjet::contrib::EnergyCorrelator::pt_R);
                fastjet::contrib::EnergyCorrelator ECF2(2, beta, fastjet::contrib::EnergyCorrelator::pt_R);
                fastjet::contrib::EnergyCorrelator ECF1(1, beta, fastjet::contrib::EnergyCorrelator::pt_R);
                double recf3 = ECF3(LJet);
                double recf2 = ECF2(LJet);
                double recf1 = ECF1(LJet);
                c2 = (recf2 != 0 ? recf3 * recf1 / (recf2 * recf2) : -1);

                // NSubjettiness
                fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
                fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
                fastjet::contrib::Nsubjettiness nSub3(3, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
                double tau1 = nSub1.result(LJet);
                double tau2 = nSub2.result(LJet);
                double tau3 = nSub3.result(LJet);
                if (tau1 != 0)
                    tau21 = tau2 / tau1;
                else
                    tau21 = -99;
                if (tau2 != 0)
                    tau32 = tau3 / tau2;
                else
                    tau32 = -99;

                // Trimmed Jet pT
                const double pt = LJet.pt();

                // Trimmed Jet Mass
                const double mass = LJet.m();

                _h["pTFJ"]->fill(pt);
                _h["dj0_tm"]->fill(mass);
                _h["dj0_c2"]->fill(c2);
                _h["dj0_lha"]->fill(lha);
                _h["dj0_t21"]->fill(tau21);
                _h["dj0_t32"]->fill(tau32);


                // Second Condition
                // If the Subleading trimmed Jet is closer to MET than Leading Trimmed Jet and both are within 1.0 
            } else if (Dphi1 < Dphi0 && Dphi1 < 1.0) {
                // Trimmed Sub-Leading fatjet (from multijet event)
                PseudoJet LJet = fjet1;

                // LHA
                for (const PseudoJet &p : LJet.constituents()) {
                    double trpt = p.pt();
                    double trtheta = p.squared_distance(LJet);
                    lha += pow(trpt, 1.0) * pow(trtheta, 0.25);
                }
                double lterm = pow(LJet.pt(), 1.0) * pow(1.0, 0.5);
                if (lterm != 0)
                    lha /= lterm;
                else
                    lha = -99;

                // Energy Correlation Function (C2)
                fastjet::contrib::EnergyCorrelator ECF3(3, beta, fastjet::contrib::EnergyCorrelator::pt_R);
                fastjet::contrib::EnergyCorrelator ECF2(2, beta, fastjet::contrib::EnergyCorrelator::pt_R);
                fastjet::contrib::EnergyCorrelator ECF1(1, beta, fastjet::contrib::EnergyCorrelator::pt_R);
                double recf3 = ECF3(LJet);
                double recf2 = ECF2(LJet);
                double recf1 = ECF1(LJet);
                c2 = (recf2 != 0 ? recf3 * recf1 / (recf2 * recf2) : -1);

                // NSubjettiness
                fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
                fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
                fastjet::contrib::Nsubjettiness nSub3(3, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta, Rcut));
                double tau1 = nSub1.result(LJet);
                double tau2 = nSub2.result(LJet);
                double tau3 = nSub3.result(LJet);
                if (tau1 != 0)
                    tau21 = tau2 / tau1;
                else
                    tau21 = -99;
                if (tau2 != 0)
                    tau32 = tau3 / tau2;
                else
                    tau32 = -99;

                // Trimmed Jet pT
                const double pt = LJet.pt();

                // Trimmed Jet Mass
                const double mass = LJet.m();

                _h["pTFJ"]->fill(pt);
                _h["dj0_tm"]->fill(mass);
                _h["dj0_c2"]->fill(c2);
                _h["dj0_lha"]->fill(lha);
                _h["dj0_t21"]->fill(tau21);
                _h["dj0_t32"]->fill(tau32);
            }
        }
    }

    void finalize()
    {
        // Area under the curve normalize to 1
        normalize(_h, 1.0);
    }

private:
    fastjet::Filter _trimmer;
    map<string, Histo1DPtr> _h;
};

RIVET_DECLARE_PLUGIN(MC_Darkmatter);

} // namespace Rivet
