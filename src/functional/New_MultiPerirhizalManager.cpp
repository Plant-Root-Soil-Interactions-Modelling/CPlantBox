#include "New_MultiPerirhizalManager.h"
#include <cmath> // std::sqrt

namespace CPlantBox
{

void New_MultiPerirhizalManager::mergeSharedCells()
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
void New_MultiPerirhizalManager::recomputeSharedOuterRadii(int type)
{
    // prepare storage per plant
    shared_outer_.clear();
    shared_outer_.resize(entries_.size());

    // initialise per-plant result arrays
    for (std::size_t p = 0; p < entries_.size(); ++p) {
        const auto& ms = entries_[p].ms;
        if (!ms) continue;
        shared_outer_[p].assign(ms->segments.size(), 0.0);
    }

    if (entries_.empty()) {
        return;
    }

    // find a reference MappedSegments to get grid info
    std::shared_ptr<MappedSegments> ref_ms = nullptr;
    for (const auto& e : entries_) {
        if (e.ms) {
            ref_ms = e.ms;
            break;
        }
    }
    if (!ref_ms) {
        return; // no valid ms at all
    }

    // assume rectangular, equidistant grid like Perirhizal::segOuterRadii
    auto width = ref_ms->maxBound.minus(ref_ms->minBound);
    double cellVolume =
        width.x * width.y * width.z /
        (ref_ms->resolution.x * ref_ms->resolution.y * ref_ms->resolution.z);

    // -------------------------------------------------------------------------
    // FIRST PASS: accumulate total weight per cell over ALL plants
    // -------------------------------------------------------------------------
    std::unordered_map<int, double> cellSumW;

    // optional: reserve roughly the number of distinct cells
    std::size_t estimatedCells = 0;
    for (const auto& e : entries_) {
        if (e.ms) {
            estimatedCells += e.ms->cell2seg.size();
        }
    }
    cellSumW.reserve(estimatedCells);

    for (std::size_t p = 0; p < entries_.size(); ++p) {
        const auto& ms = entries_[p].ms;
        if (!ms) continue;

        const auto& segLength = ms->segLength();
        const auto& effRadii = ms->getEffectiveRadii();
        const auto& organTypes = ms->organTypes;

        for (const auto& kv : ms->cell2seg) {
            int cellId = kv.first;
            const auto& segs = kv.second;

            // negative cell ids: treated separately (see second pass)
            if (cellId < 0) {
                continue;
            }

            double& sumW = cellSumW[cellId]; // default-initialises to 0 if new

            for (int sidx : segs) {
                if (organTypes[sidx] != Organism::ot_root) {
                    continue;
                }

                double l = segLength[sidx];
                double a = effRadii[sidx];

                double w = 0.0;
                if (type == 0) {            // volume
                    w = M_PI * a * a * l;
                }
                else if (type == 1) {     // surface
                    w = 2.0 * M_PI * a * l;
                }
                else {                    // length
                    w = l;
                }

                sumW += w;
            }
        }
    }

    // -------------------------------------------------------------------------
    // SECOND PASS: compute outer radii per segment for each plant
    // -------------------------------------------------------------------------
    for (std::size_t p = 0; p < entries_.size(); ++p) {
        const auto& ms = entries_[p].ms;
        if (!ms) continue;

        const auto& segLength = ms->segLength();
        const auto& effRadii = ms->getEffectiveRadii();
        const auto& organTypes = ms->organTypes;

        for (const auto& kv : ms->cell2seg) {
            int cellId = kv.first;
            const auto& segs = kv.second;

            // behaviour for cellId < 0: match Perirhizal::segOuterRadii
            if (cellId < 0) {
                for (int sidx : segs) {
                    if (organTypes[sidx] == Organism::ot_root) {
                        shared_outer_[p][sidx] = 1.0; // arbitrary fallback value
                    }
                }
                continue;
            }

            auto it = cellSumW.find(cellId);
            if (it == cellSumW.end()) {
                // no weight recorded for this cell → leave as 0.0 for roots
                continue;
            }

            double sumW = it->second;
            if (sumW <= 0.0) {
                continue;
            }

            for (int sidx : segs) {
                if (organTypes[sidx] != Organism::ot_root) {
                    // keep 0.0 for non-root segments
                    continue;
                }

                double l = segLength[sidx];
                double a = effRadii[sidx];

                double w = 0.0;
                if (type == 0) {            // volume
                    w = M_PI * a * a * l;
                }
                else if (type == 1) {     // surface
                    w = 2.0 * M_PI * a * l;
                }
                else {                    // length
                    w = l;
                }

                double t = w / sumW;          // proportionality factor
                double targetV = t * cellVolume;    // target soil volume
                double outer_r = std::sqrt(targetV / (M_PI * l) + a * a);

                shared_outer_[p][sidx] = outer_r;
            }
        }
    }
}
} // namespace CPlantBox
