#include "othercorrel.hh"
#include <cmath>
#include <stdexcept>
#include <fastjet/Selector.hh>
#include <iostream>
#include <string.h>
using namespace std;

namespace OtherCorrelators
{
    CorrelatorsContainer::CorrelatorsContainer()
    : fr()
    , fw()
    , frxw()
    , findx1()
    , findx2()
    {
        ;
    }
    
    CorrelatorsContainer::~CorrelatorsContainer()
    {
        ;
    }

    void CorrelatorsContainer::clear()
    {
        fw.clear();
        fr.clear();
        frxw.clear();
    }

    void CorrelatorsContainer::addwr(const double &w, const double &r)
    {
	cout << "weight pushed = " << w <<endl;
	cout << "distance pushed = " << r <<endl;
        fw.push_back(w);
        fr.push_back(r);
    }

    void CorrelatorsContainer::addwr(const double &w, const double &r, const int &indx1, const int &indx2)
    {

	fr.push_back(r);
        fw.push_back(w);
	findx1.push_back(indx1);
        findx2.push_back(indx2);
    }

    std::vector<double> *CorrelatorsContainer::weights()
    {
        return &fw;
    }

    std::vector<double> *CorrelatorsContainer::rs()
    {
        return &fr;
    }

    std::vector<int> *CorrelatorsContainer::indices1()
    {
        return &findx1;
    }

    std::vector<int> *CorrelatorsContainer::indices2()
    {
        return &findx2;
    }

    const double *CorrelatorsContainer::wa()
    {
        return &fw[0];
    }

    const double *CorrelatorsContainer::ra()
    {
        return &fr[0];
    }

    std::vector<double> *CorrelatorsContainer::rxw()
    {
        frxw.clear();
        for (size_t i = 0; i < fr.size(); i++)
        {
            frxw.push_back(fr[i] * fw[i]);
        }
        return &frxw;
    }

    std::vector<fastjet::PseudoJet> constituents_as_vector(const fastjet::PseudoJet &jet)
    {
        std::vector<fastjet::PseudoJet> _v;
        for (auto &c : jet.constituents())
        {
            _v.push_back(c);
        }
        return _v;
    }

    OtherCorrelatorBuilder::OtherCorrelatorBuilder()
    : fec()
    , fncmax(0)
    {
        ;
    }
    
    OtherCorrelatorBuilder::OtherCorrelatorBuilder(const std::vector<fastjet::PseudoJet> &parts,  const double &scale, const int &nmax, const int &power, const double dphi_cut = -9999, const double deta_cut = -9999, const char* correltype = "")
    : fec()
    , fncmax(nmax)
    {
        // std::cout << "Initializing n point correlator with power " << power << " for " << parts.size() << " paritlces" << std::endl;
        if (fncmax < 2)
        {
            throw std::overflow_error("asking for n-point correlator with n < 2?");
        }
        if (fncmax > 3)
        {
            throw std::overflow_error("max n for n-point correlator is currently 2");
        }
        for (int i = 0; i < fncmax - 2 + 1; i++)
        {
            fec.push_back(new CorrelatorsContainer());
        }
        for (size_t i = 0; i < parts.size(); i++)
        {
            for (size_t j = 0; j < parts.size(); j++)
            {
                double _phi12 = fabs(parts[i].delta_phi_to(parts[j])); // expecting delta_phi_to() to return values in [-pi, pi]
                double _eta12 = parts[i].eta() - parts[j].eta();
                if (dphi_cut > -1)
                { // if dphi_cut is on, apply it to pairs
                    double _pt1 = parts[i].pt();
                    double _pt2 = parts[j].pt();
                    int _q1 = 1; // FIX ME: just dummy (no charge info available yet in data and full sim)
                    int _q2 = 1;
                    if ( !ApplyDeltaPhiRejection(dphi_cut, _q1, _q2, _pt1, _pt2, _phi12) ) continue;
                }
                if (deta_cut > -1)
                { // if deta_cut is on, apply it to pairs
                    if ( !ApplyDeltaEtaRejection(deta_cut, _eta12) ) continue;
                }
                double _d12 = 0;
                
                cout << "correltype " << correltype << "; " << strcmp(correltype, "deltap") << endl;

                if (strcmp(correltype, "deltap") == 0) {
                    // _d12 = fabs(parts[i].p() - parts[j].p()); //parts[i].delta_R(parts[j]);
                    double p1 = sqrt(parts[i].px()*parts[i].px() + parts[i].py()*parts[i].py() + parts[i].pz()*parts[i].pz()); //parts[i].charge();
                    double p2 = sqrt(parts[j].px()*parts[j].px() + parts[j].py()*parts[j].py() + parts[j].pz()*parts[j].pz());
                    _d12 = fabs(p1 - p2);
                    cout << "delta p: " << p1 << " - " << p2 << " = " << _d12 << endl;
                } else if (strcmp(correltype, "deltapt") == 0) {
                    _d12 = fabs(parts[i].pt() - parts[j].pt());
                    cout << "delta pT: " << parts[i].pt() << " - " << parts[j].pt() << " = " << _d12 << endl;
                } // else if (correltype == "charge") {
                //     double q1 = parts[i].python_info().charge; //parts[i].charge();
                //     double q2 = parts[j].python_info().charge; //parts[j].charge();
                //     _d12 = q1*q2;
                // }
		
		
                double _w2 = 0;
            
                // _w2 = parts[i].perp() * parts[j].perp() / std::pow(scale, 2);
                // _w2 = pow(_w2, power);

                fec[2 - 2]->addwr(_w2, _d12, i, j); // save weight, distance and indices of the pair
                
                if (fncmax < 3)
                    continue;
                
            }
        }
    }

    



    CorrelatorsContainer* OtherCorrelatorBuilder::correlator(int n)
    {
        if (n > fncmax)
        {
            throw std::overflow_error("requesting n-point correlator with too large n");
        }
        if (n < 2)
        {
            throw std::overflow_error("requesting n-point correlator with n < 2?");
        }
        return fec[n - 2];
    }

    OtherCorrelatorBuilder::~OtherCorrelatorBuilder()
    {
        for (auto p : fec)
        {
            delete p;
        }
        fec.clear();
    }

    bool OtherCorrelatorBuilder::ApplyDeltaPhiRejection(const double dphi_cut, const double q1, const double q2, const double pt1, const double pt2, const double phi12)
    {
        double R = 1.1; // reference radius for TPC
        double Bz = 0.5;
        double phi_star = phi12 + q1*asin(-0.015*Bz*R/pt1) - q2*asin(-0.015*Bz*R/pt2);
        if ( fabs(phi_star)<dphi_cut ) return false;  
        return true;
    }

    bool OtherCorrelatorBuilder::ApplyDeltaEtaRejection(const double deta_cut, const double eta12)
    {
        if ( fabs(eta12)<deta_cut ) return false;
        return true;
    }

	std::vector<fastjet::PseudoJet> merge_signal_background_pjvectors(const std::vector<fastjet::PseudoJet> &signal, 
																	  const std::vector<fastjet::PseudoJet> &background,
																      const double pTcut,
																	  const int bg_index_start)
    {
        std::vector<fastjet::PseudoJet> _vreturn;
        auto _selector = fastjet::SelectorPtMin(pTcut);
        auto _signal = _selector(signal);
        for (auto &_p : _signal)
        {
            _p.set_user_index(_p.user_index());
            _vreturn.push_back(_p);
        }
        auto _background = _selector(background);
        for (auto &_p : _background)
        {
            if (bg_index_start > 0)
            {
                int _index = &_p - &_background[0];
                _p.set_user_index(bg_index_start + _index);
            }
            else
            {
                _p.set_user_index(_p.user_index());
            }
            _vreturn.push_back(_p);
        }
        return _vreturn;
    }

}
