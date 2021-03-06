#include "PhysicsTools/PatAlgos/plugins/PATSingleVertexSelector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include <DataFormats/BeamSpot/interface/BeamSpot.h>

#include <algorithm>

using pat::PATSingleVertexSelector;


PATSingleVertexSelector::Mode
PATSingleVertexSelector::parseMode(const std::string &mode) {
    if (mode == "firstVertex") {
        return First;
    } else if (mode == "nearestToCandidate") {
        return NearestToCand;
    } else if (mode == "fromCandidate") {
        return FromCand;
    } else if (mode == "beamSpot") {
        return FromBeamSpot;
    } else {
        throw cms::Exception("Configuration") << "PATSingleVertexSelector: Mode '" << mode << "' not recognized or not supported.\n";
    }
}


PATSingleVertexSelector::PATSingleVertexSelector(const edm::ParameterSet & iConfig)
  : doFilterEvents_(false)
{
   using namespace std;

   modes_.push_back( parseMode(iConfig.getParameter<std::string>("mode")) );
   if (iConfig.exists("fallbacks")) {
      vector<string> modes = iConfig.getParameter<vector<string> >("fallbacks");
      for (vector<string>::const_iterator it = modes.begin(), ed = modes.end(); it != ed; ++it) {
        modes_.push_back( parseMode(*it) );
      }
   }
   if (hasMode_(First) || hasMode_(NearestToCand)) {
        vertices_ = iConfig.getParameter<edm::InputTag>("vertices");
        if (iConfig.existsAs<string>("vertexPreselection")) {
            string presel = iConfig.getParameter<string>("vertexPreselection");
            if (!presel.empty()) vtxPreselection_ = auto_ptr<VtxSel>(new VtxSel(presel));
        }
   }
   if (hasMode_(NearestToCand) || hasMode_(FromCand)) {
        candidates_ = iConfig.getParameter<vector<edm::InputTag> >("candidates");
        if (iConfig.existsAs<string>("candidatePreselection")) {
            string presel = iConfig.getParameter<string>("candidatePreselection");
            if (!presel.empty()) candPreselection_ = auto_ptr<CandSel>(new CandSel(presel));
        }
   }
   if (hasMode_(FromBeamSpot)) {
        beamSpot_ = iConfig.getParameter<edm::InputTag>("beamSpot");
   }

   if ( iConfig.exists("filter") ) doFilterEvents_ = iConfig.getParameter<bool>("filter");

   produces<vector<reco::Vertex> >();
}


PATSingleVertexSelector::~PATSingleVertexSelector()
{
}

bool PATSingleVertexSelector::hasMode_(Mode mode) const {
    return (std::find(modes_.begin(), modes_.end(), mode) != modes_.end());
}

bool
PATSingleVertexSelector::filter(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    // Clear
    selVtxs_.clear(); bestCand_ = 0;

    // Gather data from the Event
    // -- vertex data --
    if (hasMode_(First) || hasMode_(NearestToCand)) {
        Handle<vector<reco::Vertex> > vertices;
        iEvent.getByLabel(vertices_, vertices);
        for (vector<reco::Vertex>::const_iterator itv = vertices->begin(), edv = vertices->end(); itv != edv; ++itv) {
            if ((vtxPreselection_.get() != 0) && !((*vtxPreselection_)(*itv)) ) continue;
            selVtxs_.push_back( &*itv );
        }
    }
    // -- candidate data --
    if (hasMode_(NearestToCand) || hasMode_(FromCand)) {
       vector<pair<double, const reco::Candidate *> > cands;
       for (vector<edm::InputTag>::const_iterator itt = candidates_.begin(), edt = candidates_.end(); itt != edt; ++itt) {
          Handle<View<reco::Candidate> > theseCands;
          iEvent.getByLabel(*itt, theseCands);
          for (View<reco::Candidate>::const_iterator itc = theseCands->begin(), edc = theseCands->end(); itc != edc; ++itc) {
            if ((candPreselection_.get() != 0) && !((*candPreselection_)(*itc))) continue;
            cands.push_back( pair<double, const reco::Candidate *>(-itc->pt(), &*itc) );
          }
       }
       if (!cands.empty()) bestCand_ = cands.front().second;
    }

    bool passes = false;
    auto_ptr<vector<reco::Vertex> > result;
    // Run main mode + possible fallback modes
    for (std::vector<Mode>::const_iterator itm = modes_.begin(), endm = modes_.end(); itm != endm; ++itm) {
        result = filter_(*itm, iEvent, iSetup);
        // Check if we got any vertices.  If so, take them.
        if (result->size()) {
          passes = true;
          break;
        }
    }
    iEvent.put(result);
    // Check if we want to apply the EDFilter
    if (doFilterEvents_)
      return passes;
    else return true;
}

std::auto_ptr<std::vector<reco::Vertex> >
PATSingleVertexSelector::filter_(Mode mode, const edm::Event &iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;
    std::auto_ptr<std::vector<reco::Vertex> > result(
        new std::vector<reco::Vertex>());
    switch(mode) {
        case First: {
            if (selVtxs_.empty()) return result;
            result->push_back(*selVtxs_.front());
            return result;
            }
        case FromCand: {
            if (bestCand_ == 0) return result;
            reco::Vertex vtx;
            if (typeid(*bestCand_) == typeid(reco::VertexCompositeCandidate)) {
                vtx = reco::Vertex(bestCand_->vertex(), bestCand_->vertexCovariance(),
                        bestCand_->vertexChi2(), bestCand_->vertexNdof(), bestCand_->numberOfDaughters() );
            } else {
                vtx = reco::Vertex(bestCand_->vertex(), reco::Vertex::Error(), 0, 0, 0);
            }
            result->push_back(vtx);
            return result;
            }
        case NearestToCand: {
            if (selVtxs_.empty() || (bestCand_ == 0)) return result;
            const reco::Vertex * which = 0;
            float dzmin = 9999.0;
            for (vector<const reco::Vertex *>::const_iterator itv = selVtxs_.begin(), edv = selVtxs_.end(); itv != edv; ++itv) {
                float dz = std::abs((*itv)->z() - bestCand_->vz());
                if (dz < dzmin) { dzmin = dz; which = *itv; }
            }
            if (which != 0) // actually it should not happen, but better safe than sorry
              result->push_back(*which);
            return result;
            }
        case FromBeamSpot: {
            Handle<reco::BeamSpot> beamSpot;
            iEvent.getByLabel(beamSpot_, beamSpot);
            reco::Vertex bs(beamSpot->position(), beamSpot->covariance3D(), 0, 0, 0);
            result->push_back(bs);
            return result;
            }
        default:
            // Return an empty vector signifying no vertices found.
            return result;
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATSingleVertexSelector);
