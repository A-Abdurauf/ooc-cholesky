#pragma once
#include <unordered_map>
namespace plasma::dist {
using RankMapping = std::unordered_map<int, int>;

int PLASMA_Init_GPU_Dist(const RankMapping *rankMapping);

int PLASMA_Init_Affinity_GPU_Dist(const RankMapping *rankMapping,
                                  int *coresbind);

int PLASMA_Finalize_GPU_Dist();

RankMapping getVerifyRankMapping();
}  // namespace plasma::dist
