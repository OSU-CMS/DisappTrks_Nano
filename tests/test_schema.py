from disapptrks.schema import (
    EVENT_ID_BRANCHES,
    ISOTRACK_BRANCHES,
    JET_BRANCHES,
    STORED_NO_MU_MET_BRANCHES,
    TAG_AND_PROBE_BRANCHES,
    audit_branches,
)


def test_complete_background_schema_is_ready():
    branches = (
        EVENT_ID_BRANCHES
        | ISOTRACK_BRANCHES
        | JET_BRANCHES
        | STORED_NO_MU_MET_BRANCHES
        | TAG_AND_PROBE_BRANCHES
    )
    report = audit_branches(branches)
    assert report.ready_for("backgrounds")
    assert report.has_no_mu_met_inputs


def test_incomplete_isotrack_schema_is_rejected():
    report = audit_branches({"run", "event", "nIsoTrack", "IsoTrack_pt"})
    assert not report.ready_for("search")
    assert not report.has_no_mu_met_inputs
