#ifndef EffAnalyzer_TreeProducer_ProtonKinematicsUtils_h
#define EffAnalyzer_TreeProducer_ProtonKinematicsUtils_h

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "EffAnalyzer/TreeProducer/interface/RPAlignmentConstants.h"

#include <utility> // std::pair

namespace ProtonUtils
{
  typedef std::pair<unsigned int, const TotemRPLocalTrack&> localtrack_t;

  float tracksDistance( const CTPPSAlCa::RPAlignmentConstants&, const localtrack_t& near_p, const localtrack_t& far_p );
}

#endif
