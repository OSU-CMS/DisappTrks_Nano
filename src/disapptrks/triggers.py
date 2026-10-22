"""Run-dependent trigger definitions used by the tau background estimate."""

from __future__ import annotations


ISO_MUON_REFERENCE_TRIGGER = "IsoMu24"

# HPS30 (IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1) is not
# prescaled (brilcalc: HLT and L1-seed prescale 1) but fires ~1/330 as often
# as HPS27 in 2022-2024 data, so it is not usable as the nominal path.  HPS27
# (the dissertation/legacy choice) is used for 2022-2023; PNet26 exists from
# 2024 onward and is used there.
_TAU_CROSS_TRIGGERS = {
    "2022": "IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
    "2023": "IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
    "2024": "IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1",
    "2025": "IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1",
    "2026": "IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1",
}


def tau_cross_trigger_for_year(year: str) -> str:
    """Return the nominal unprescaled muon-tau cross-trigger for a Run-3 year."""

    calendar_year = str(year).split("_", 1)[0]
    try:
        return _TAU_CROSS_TRIGGERS[calendar_year]
    except KeyError as error:
        raise ValueError(
            f"no nominal tau cross-trigger is configured for year {year!r}"
        ) from error
