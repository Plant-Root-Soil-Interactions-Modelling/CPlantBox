#pragma once

#include <vector>
#include <memory>
#include <unordered_map>

#include "Perirhizal.h"
#include "MappedOrganism.h"

namespace CPlantBox
{

class New_MultiPerirhizalManager {
public:
    struct Entry {
        std::shared_ptr<Perirhizal>     peri;
        std::shared_ptr<MappedSegments> ms;
        std::vector<double>             hsr;   // interface heads proposed by this plant
    };

    New_MultiPerirhizalManager() = default;

    void addPlant(const std::shared_ptr<Perirhizal>& peri,
                  const std::shared_ptr<MappedSegments>& ms)
    {
        Entry e;
        e.peri = peri;
        e.ms   = ms;
        entries_.push_back(e);

        // keep outer-radii storage in sync
        shared_outer_.emplace_back();  // will be filled in recomputeSharedOuterRadii
    }

    void setProposal(std::size_t plantIdx, const std::vector<double>& hsr)
    {
        entries_.at(plantIdx).hsr = hsr;
    }

    // existing one
    void mergeSharedCells();

    // NEW: compute per-segment outer radii, but across ALL plants
    // type = 0 volume, 1 surface, 2 length (same as Perirhizal::segOuterRadii)
    void recomputeSharedOuterRadii(int type = 2);

    const std::vector<double>& getMerged(std::size_t plantIdx) const
    {
        return entries_.at(plantIdx).hsr;
    }

    // NEW: Python can fetch the shared outer radii for plant i
    const std::vector<double>& getSharedOuter(std::size_t plantIdx) const
    {
        return shared_outer_.at(plantIdx);
    }

    std::size_t size() const { return entries_.size(); }

private:
    std::vector<Entry> entries_;
    // per plant: per-segment shared outer radius
    std::vector<std::vector<double>> shared_outer_;
};

} // namespace CPlantBox
