import ROOT

ORTHO_LUMI_CPP = r"""
#include <cstdint>

namespace ortholumi {

    static uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    static double uniform01(unsigned int run,
                            unsigned int lumi,
                            unsigned long long event,
                            unsigned long long fullEventId,
                            unsigned int seed = 12345) {
        uint64_t x = uint64_t(seed);

        x ^= splitmix64(uint64_t(run));
        x ^= splitmix64(uint64_t(lumi) << 16);
        x ^= splitmix64(uint64_t(event));
        x ^= splitmix64(uint64_t(fullEventId));

        uint64_t h = splitmix64(x);

        return (h >> 11) * (1.0 / 9007199254740992.0);
    }

    static int era_tag_2024_2025_2026(unsigned int run,
                                      unsigned int lumi,
                                      unsigned long long event,
                                      unsigned long long fullEventId,
                                      unsigned int seed = 12345) {
        const double lumi_2024 = 110.0;
        const double lumi_2025 = 111.0;
        const double lumi_2026 = 26.0;

        const double total = lumi_2024 + lumi_2025 + lumi_2026;

        const double edge_2024 = lumi_2024 / total;
        const double edge_2025 = (lumi_2024 + lumi_2025) / total;

        const double u = uniform01(run, lumi, event, fullEventId, seed);

        if (u < edge_2024) return 0;
        if (u < edge_2025) return 1;
        return 2;
    }
}
"""

_declared = False


def _declare_once():
    global _declared
    if not _declared:
        ROOT.gInterpreter.Declare(ORTHO_LUMI_CPP)
        _declared = True


def era_to_orthogonal_tag(era):
    era = str(era)
    if "2024" in era:
        return 0
    if "2025" in era:
        return 1
    if "2026" in era:
        return 2

    return None


def ApplyOrthogonalLumiFilter(df, era, seed=12345, keep_tag_column=True):
    _declare_once()

    target_tag = era_to_orthogonal_tag(era)

    if target_tag is None:
        return df, []

    df = df.Define(
        "OrthogonalEraTag",
        (
            "ortholumi::era_tag_2024_2025_2026("
            f"run, luminosityBlock, event, FullEventId, {int(seed)})"
        )
    )

    df = df.Filter(
        f"OrthogonalEraTag == {target_tag}",
        f"Orthogonal lumi split for {era}"
    )

    cols = ["OrthogonalEraTag"] if keep_tag_column else []

    return df, cols