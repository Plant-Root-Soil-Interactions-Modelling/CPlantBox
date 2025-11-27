#include "MultiPerirhizalManager.h"
#include <cmath> // std::sqrt

namespace CPlantBox
{

void MultiPerirhizalManager::mergeSharedCells()
{
    // cell_id -> list of (plant_index, seg_index)
    std::unordered_map<int, std::vector<std::pair<std::size_t,int>>> cell2list;

    // 1) build the global map using the already-grouped per-plant cell2seg
    for (std::size_t p = 0; p < entries_.size(); ++p) {
        const auto& ms  = entries_[p].ms;
        if (!ms) continue;

        // ms->cell2seg is kept up to date by MappedSegments::mapSegments(...) etc.
        for (const auto& kv : ms->cell2seg) {
            int cellId = kv.first;
            const auto& segs = kv.second;
            for (int sidx : segs) {
                cell2list[cellId].push_back({p, sidx});
            }
        }
    }

    // 2) merge cells that are touched by more than one segment (possibly from different plants)
    for (auto& kv : cell2list) {
        auto& entries = kv.second;
        if (entries.size() < 2)
            continue;

        double wsum = 0.0;
        double hsum = 0.0;

        // weighted average of interface heads, weight = segment length
        for (const auto& ps : entries) {
            std::size_t p = ps.first;
            int sidx      = ps.second;
            auto ms       = entries_[p].ms;

            double w = ms->segLength()[sidx];
            double h = entries_[p].hsr[sidx];

            wsum += w;
            hsum += w * h;
        }

        if (wsum > 0.0) {
            double merged = hsum / wsum;
            for (const auto& ps : entries) {
                std::size_t p = ps.first;
                int sidx      = ps.second;
                entries_[p].hsr[sidx] = merged;
            }
        }
    }
}


// -----------------------------------------------------------------------------
// NEW: true-ish parity — compute outer radii over ALL plants together
// -----------------------------------------------------------------------------
void MultiPerirhizalManager::recomputeSharedOuterRadii(int type)
{
    // global: cell -> list of (plant, seg)
    std::unordered_map<int, std::vector<std::pair<std::size_t,int>>> cell2list;

    // prepare storage per plant
    shared_outer_.clear();
    shared_outer_.resize(entries_.size());

    // build global map and init per-plant arrays
    for (std::size_t p = 0; p < entries_.size(); ++p) {
        const auto& ms = entries_[p].ms;
        if (!ms) continue;

        // allocate result array for this plant
        shared_outer_[p].assign(ms->segments.size(), 0.0);

        // we use cell2seg here because it's already grouped per cell
        for (const auto& kv : ms->cell2seg) {
            int cellId = kv.first;
            const auto& segs = kv.second;
            for (int segIdx : segs) {
                cell2list[cellId].push_back({p, segIdx});
            }
        }
    }

    if (entries_.empty())
        return;

    // assume rectangular grid like in Perirhizal::segOuterRadii
    auto ref_ms = entries_[0].ms;
    auto width  = ref_ms->maxBound.minus(ref_ms->minBound);
    double cellVolume = width.x * width.y * width.z /
                        (ref_ms->resolution.x * ref_ms->resolution.y * ref_ms->resolution.z);

    // now do the real sharing
    for (auto& kv : cell2list) {
        int cellId = kv.first;
        auto& seglist = kv.second;

        // collect weights across all plants for this cell
        double sumW = 0.0;
        std::vector<double> weights;
        weights.reserve(seglist.size());

        // first pass: compute total weight
        for (auto& ps : seglist) {
            std::size_t p = ps.first;
            int sidx      = ps.second;
            auto ms       = entries_[p].ms;

            // we only apply to roots, like Perirhizal does
            if (ms->organTypes[sidx] != Organism::ot_root) {
                weights.push_back(0.0);
                continue;
            }

            double l = ms->segLength()[sidx];
            double a = ms->getEffectiveRadii()[sidx];

            double w = 0.0;
            if (type == 0) {            // volume
                w = M_PI * a * a * l;
            } else if (type == 1) {     // surface
                w = 2.0 * M_PI * a * l;
            } else {                    // length
                w = l;
            }
            weights.push_back(w);
            sumW += w;
        }

        if (sumW <= 0.0) {
            continue;
        }

        // second pass: assign outer radii like in Perirhizal::segOuterRadii
        for (std::size_t k = 0; k < seglist.size(); ++k) {
            auto [p, sidx] = seglist[k];
            auto ms        = entries_[p].ms;

            if (ms->organTypes[sidx] != Organism::ot_root) {
                // keep 0.0 in shared_outer_ for non-roots
                continue;
            }

            double l = ms->segLength()[sidx];
            double a = ms->getEffectiveRadii()[sidx];

            double t = weights[k] / sumW;       // proportionality factor
            double targetV = t * cellVolume;    // target soil vol. for this segment

            double outer_r = std::sqrt(targetV / (M_PI * l) + a * a);

            shared_outer_[p][sidx] = outer_r;
        }
    }
}

} // namespace CPlantBox
