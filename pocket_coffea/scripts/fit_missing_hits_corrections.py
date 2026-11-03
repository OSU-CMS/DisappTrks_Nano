"""Fit the missing-hits (hit-drop) correction probabilities (SIG-02/SIG-03).

Port of the legacy DisappTrks/StandardAnalysis/scripts/
findMissingMiddleHitsCorrection.py and findMissingOuterHitsCorrection.py
(Run 2 method, AN Appendix B) to the missing_hits_control outputs:

- P0 (extra missing middle hit): data missing-middle distribution in the
  muon control region without the missing inner/middle requirements vs MC
  (missing middle x strip layers), fit over bins 0-5.
- P_TOB, P_pre-TOB: full muon control region, data missing-outer vs MC
  (missing outer x TOB layers x strip layers), simultaneous fit over bins 1-6.
- P_post-TOB: same histograms, fit over bins 6-9. As in the legacy script this
  second fit also floats P_TOB, but that value is discarded; P_TOB is the one
  from the pre-TOB fit (the discarded value is reported for reference).

Objective, ranges, starting values, limits and the correction algorithms are
the legacy ones (chi2 per bin used, data and MC errors in quadrature, both
unit-normalized; Minuit2 Migrad). Uncertainties (absent in the legacy
scripts) are HESSE errors of the summed chi2 (errordef 1), which has the same
minimum. ``--closure`` refits pseudo-data built from MC with known
probabilities to validate the port.

Run from pocket_coffea/ inside ./shell:
    python scripts/fit_missing_hits_corrections.py [--closure]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from coffea.util import load
from iminuit import Minuit

MIDDLE = "missingHitsNoInnerMiddle_missingMiddle_stripLayers"
OUTER = "missingHitsMuonCtrl_missingOuter_TOBLayers_stripLayers"

# label: (data period, data tag, MC period, MC tag)
FITS = {
    "2022 CD": ("2022CD", "data_v1", "2022CD", "bkgmc_v1"),
    "2022 EFG": ("2022EFG", "data_v1", "2022EFG", "bkgmc_v1"),
    "2023 C": ("2023C", "data_v1", "2023C", "bkgmc_v1"),
    "2023 D": ("2023D", "data_v1", "2023D", "bkgmc_v1"),
    "2024": ("2024", "data_v1", "2023D", "bkgmc_v1_pu2024"),
    "2025": ("2025", "data_v1", "2023D", "bkgmc_v1_pu2025"),
    "2026": ("2026", "data_v1", "2023D", "bkgmc_v1_pu2026"),
}

MIDDLE_BINS = range(0, 6)  # legacy BIN_RANGE
NMIDDLE = 6  # corrected middle histogram: 6 bins, -0.5..5.5
NOUTER = 16  # corrected outer histogram: 16 bins, -0.5..15.5
PRE_TOB_BINS = range(1, 7)
POST_TOB_BINS = range(6, 10)
STRIP_LAYERS = range(1, 20)  # legacy "for num_hits in range(1, 20)"


# ----------------------------------------------------------------- inputs

def _leaves(node):
    if isinstance(node, dict):
        for child in node.values():
            yield from _leaves(child)
    else:
        yield node


def _summed(output, name):
    """Sum over samples/datasets; nominal variation; values and variances."""
    values = variances = None
    for hist in _leaves(output["variables"][name]):
        h = hist[{"cat": hist.axes["cat"][0]}]
        if "variation" in [a.name for a in h.axes]:
            h = h[{"variation": "nominal"}]
        v, w2 = h.values(), h.variances()
        if w2 is None:
            w2 = v.copy()
        values = v if values is None else values + v
        variances = w2 if variances is None else variances + w2
    return values, variances


def load_inputs(data_period, data_tag, mc_period, mc_tag):
    base = "analysis_output/{}/missing_hits_control/{}/output_all.coffea"
    data = load(base.format(data_period, data_tag))
    mc = load(base.format(mc_period, mc_tag))
    data_middle, _ = _summed(data, MIDDLE)
    data_outer, _ = _summed(data, OUTER)
    return {
        "data_middle": data_middle.sum(axis=1),  # project to missing middle
        "data_outer": data_outer.sum(axis=(1, 2)),  # project to missing outer
        "mc_middle": _summed(mc, MIDDLE),  # [missing middle, strip layers]
        "mc_outer": _summed(mc, OUTER),  # [missing outer, TOB layers, strip layers]
    }


def data_pdf(counts):
    """Unit-normalized data distribution with Poisson errors (TH1::Scale)."""
    total = counts.sum()
    return counts / total, np.sqrt(counts) / total


# ------------------------------------------------------- legacy algorithms

def middle_extra_pdf(n, p):
    """Legacy get_missing_hits_pdf (middle): P(k extra missing middle hits)."""
    pdf = np.zeros(n + 1)
    for k in range(n + 1):
        prob = 0.0
        for i in range(1, n + 1):
            if k > n - i:
                continue
            prob += (p ** (i - 1)) * (1 - p) * math.comb(n - i, k) * (p**k) * ((1 - p) ** (n - i - k))
        pdf[k] = prob
    return pdf / pdf.sum()


def corrected_middle(p, values, variances, normalize=True):
    """Legacy get_corrected_pdf (middle) on a [missing middle, strip layers] histogram."""
    out = np.zeros(NMIDDLE)
    var = np.zeros(NMIDDLE)
    for n in STRIP_LAYERS:
        if n >= values.shape[1]:
            continue
        pdf = middle_extra_pdf(n, p)
        extra = np.zeros(NMIDDLE)
        extra[: min(len(pdf), NMIDDLE)] = pdf[:NMIDDLE]
        for orig in MIDDLE_BINS:
            tracks, tracks_var = values[orig, n], variances[orig, n]
            if tracks == 0 and tracks_var == 0:
                continue
            for e in MIDDLE_BINS:
                target = orig + e
                if target < NMIDDLE:  # beyond 5.5 is overflow in the legacy histogram
                    out[target] += tracks * extra[e]
                    var[target] += tracks_var * extra[e] ** 2
    if normalize and out.sum() > 0:
        scale = 1.0 / out.sum()
        return out * scale, np.sqrt(var) * scale
    return out, np.sqrt(var)


def outer_extra_pdf(n, p):
    """Legacy get_missing_hits_pdf (outer): 16 bins, normalized over 0..15."""
    pdf = np.zeros(NOUTER)
    for k in range(min(n, NOUTER - 1) + 1):
        pdf[k] = (p**k) * (1 - p)
    if n < NOUTER:
        pdf[n] = p**n  # all hits missed
    return pdf / pdf.sum()


_O, _T, _E = np.meshgrid(np.arange(NOUTER), np.arange(NOUTER), np.arange(NOUTER), indexing="ij")


def corrected_outer(p, p_tob, values, variances, tob_status, normalize=True):
    """Legacy get_corrected_pdf (outer) on a [missing outer, TOB, strip] histogram."""
    if tob_status == "pre":
        idx_kept, idx_dropped = _O + _E, _O + _T
    else:
        idx_kept, idx_dropped = _O, _O + _E + _T
    out = np.zeros(NOUTER)
    var = np.zeros(NOUTER)
    for n in STRIP_LAYERS:
        if n >= values.shape[2]:
            continue
        w = values[:NOUTER, :NOUTER, n]
        w2 = variances[:NOUTER, :NOUTER, n]
        if not w.any() and not w2.any():
            continue
        pdf = outer_extra_pdf(n, p)
        for frac, idx in ((1.0 - p_tob, idx_kept), (p_tob, idx_dropped)):
            weight = w[:, :, None] * pdf[None, None, :] * frac
            weight2 = w2[:, :, None] * (pdf[None, None, :] * frac) ** 2
            inside = idx < NOUTER  # beyond 15.5 is overflow in the legacy histogram
            out += np.bincount(idx[inside], weights=weight[inside], minlength=NOUTER)
            var += np.bincount(idx[inside], weights=weight2[inside], minlength=NOUTER)
    if normalize and out.sum() > 0:
        scale = 1.0 / out.sum()
        return out * scale, np.sqrt(var) * scale
    return out, np.sqrt(var)


def chi2(data, data_err, mc, mc_err, bin_range):
    """Legacy get_chi2 summed (returns sum and number of bins used)."""
    total, used = 0.0, 0
    for i in bin_range:
        err = math.hypot(data_err[i], mc_err[i])
        if err > 1e-5:
            total += (data[i] - mc[i]) ** 2 / err**2
            used += 1
    return total, used


# ------------------------------------------------------------------ fits

def _minimize(fcn, starts, steps, names):
    m = Minuit(fcn, *starts, name=names)
    for name, step in zip(names, steps):
        m.errors[name] = step
        m.limits[name] = (0.0, 1.0)
    m.errordef = Minuit.LEAST_SQUARES
    m.tol = 1e-3
    m.migrad(ncall=20000)
    m.hesse()
    return m


def fit_middle(d_pdf, d_err, mc_values, mc_variances):
    def fcn(p0):
        mc, mc_err = corrected_middle(p0, mc_values, mc_variances)
        return chi2(d_pdf, d_err, mc, mc_err, MIDDLE_BINS)[0]

    m = _minimize(fcn, [0.005], [0.001], ["P0"])
    raw = mc_values[:, :].sum(axis=1)
    raw_err = np.sqrt(mc_variances.sum(axis=1))
    before, used = chi2(d_pdf, d_err, raw / raw.sum(), raw_err / raw.sum(), MIDDLE_BINS)
    after, _ = chi2(d_pdf, d_err, *corrected_middle(m.values["P0"], mc_values, mc_variances), MIDDLE_BINS)
    return {
        "P0": (m.values["P0"], m.errors["P0"]),
        "chi2_per_bin_before": before / used,
        "chi2_per_bin_after": after / used,
        "valid": bool(m.valid),
    }


def fit_outer(d_pdf, d_err, mc_values, mc_variances):
    def pre(p, p_tob):
        mc, mc_err = corrected_outer(p, p_tob, mc_values, mc_variances, "pre")
        return chi2(d_pdf, d_err, mc, mc_err, PRE_TOB_BINS)[0]

    def post(p, p_tob):
        mc, mc_err = corrected_outer(p, p_tob, mc_values, mc_variances, "post")
        return chi2(d_pdf, d_err, mc, mc_err, POST_TOB_BINS)[0]

    m_pre = _minimize(pre, [5e-3, 5e-4], [1e-3, 1e-4], ["P_preTOB", "P_TOB"])
    m_post = _minimize(post, [5e-3, 5e-4], [1e-3, 1e-4], ["P_postTOB", "P_TOB_postfit"])
    raw = mc_values.sum(axis=(1, 2))[:NOUTER]
    raw_err = np.sqrt(mc_variances.sum(axis=(1, 2)))[:NOUTER]
    raw_pdf, raw_pdf_err = raw / raw.sum(), raw_err / raw.sum()
    results = {
        "P_TOB": (m_pre.values["P_TOB"], m_pre.errors["P_TOB"]),
        "P_preTOB": (m_pre.values["P_preTOB"], m_pre.errors["P_preTOB"]),
        "P_postTOB": (m_post.values["P_postTOB"], m_post.errors["P_postTOB"]),
        "P_TOB_postfit_discarded": (m_post.values["P_TOB_postfit"], m_post.errors["P_TOB_postfit"]),
        "valid": bool(m_pre.valid and m_post.valid),
    }
    for label, bins, mode, minuit, names in (
        ("pre", PRE_TOB_BINS, "pre", m_pre, ("P_preTOB", "P_TOB")),
        ("post", POST_TOB_BINS, "post", m_post, ("P_postTOB", "P_TOB_postfit")),
    ):
        before, used = chi2(d_pdf, d_err, raw_pdf, raw_pdf_err, bins)
        corrected = corrected_outer(minuit.values[names[0]], minuit.values[names[1]],
                                    mc_values, mc_variances, mode)
        after, _ = chi2(d_pdf, d_err, *corrected, bins)
        results[f"chi2_per_bin_{label}_before"] = before / used
        results[f"chi2_per_bin_{label}_after"] = after / used
    return results


def fit_era(inputs):
    d_mid, d_mid_err = data_pdf(inputs["data_middle"])
    d_out, d_out_err = data_pdf(inputs["data_outer"][:NOUTER])
    return {
        "middle": fit_middle(d_mid, d_mid_err, *inputs["mc_middle"]),
        "outer": fit_outer(d_out, d_out_err, *inputs["mc_outer"]),
        "data_tracks": {"middle": float(inputs["data_middle"].sum()), "outer": float(inputs["data_outer"].sum())},
    }


# --------------------------------------------------------------- closure

def closure(inputs, truth, scale_to):
    """Refit pseudo-data = MC corrected with known probabilities.

    Pseudo-data have the data's track count (Poisson errors of that size). The
    legacy outer model differs between the pre-TOB and post-TOB fits, so each
    fit gets pseudo-data from its own model.
    """
    mid_values, mid_vars = inputs["mc_middle"]
    out_values, out_vars = inputs["mc_outer"]
    n_mid, n_out = scale_to

    pseudo_mid, _ = corrected_middle(truth["P0"], mid_values, mid_vars)
    d_mid = np.concatenate([pseudo_mid, np.zeros(10)]) * n_mid
    middle = fit_middle(*data_pdf(d_mid), mid_values, mid_vars)

    results = {"P0": middle["P0"]}
    for mode, prob in (("pre", "P_preTOB"), ("post", "P_postTOB")):
        pseudo, _ = corrected_outer(truth[prob], truth["P_TOB"], out_values, out_vars, mode)
        outer = fit_outer(*data_pdf(pseudo * n_out), out_values, out_vars)
        results[prob] = outer[prob]
        if mode == "pre":
            results["P_TOB"] = outer["P_TOB"]
    return results


# ------------------------------------------------------------------ main

def _fmt(value_err):
    value, err = value_err
    return f"{value:.3e} +- {err:.1e}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--eras", nargs="+", default=list(FITS))
    parser.add_argument("--closure", action="store_true")
    parser.add_argument("--output", default="analysis_output/missing_hits_corrections/fit_results_v1.json")
    args = parser.parse_args()

    results = {}
    for era in args.eras:
        data_period, data_tag, mc_period, mc_tag = FITS[era]
        try:
            inputs = load_inputs(data_period, data_tag, mc_period, mc_tag)
        except FileNotFoundError as exc:
            print(f"--- {era}: skipped, input missing ({exc.filename})")
            continue
        result = fit_era(inputs)
        result["inputs"] = {"data": f"{data_period}/{data_tag}", "mc": f"{mc_period}/{mc_tag}"}
        results[era] = result
        mid, out = result["middle"], result["outer"]
        print(f"--- {era} (data {data_period}/{data_tag}, MC {mc_period}/{mc_tag})")
        print(f"   P0        {_fmt(mid['P0'])}   chi2/bin {mid['chi2_per_bin_before']:.1f} -> {mid['chi2_per_bin_after']:.2f}")
        print(f"   P_TOB     {_fmt(out['P_TOB'])}   (post fit, discarded: {_fmt(out['P_TOB_postfit_discarded'])})")
        print(f"   P_preTOB  {_fmt(out['P_preTOB'])}   chi2/bin {out['chi2_per_bin_pre_before']:.1f} -> {out['chi2_per_bin_pre_after']:.2f}")
        print(f"   P_postTOB {_fmt(out['P_postTOB'])}   chi2/bin {out['chi2_per_bin_post_before']:.1f} -> {out['chi2_per_bin_post_after']:.2f}")
        print(f"   fits valid: middle {mid['valid']}, outer {out['valid']}")

        if args.closure:
            truth = {"P0": 0.005, "P_TOB": 8e-4, "P_preTOB": 4e-3, "P_postTOB": 0.6}
            test = closure(inputs, truth, (inputs["data_middle"].sum(), inputs["data_outer"].sum()))
            print(f"   closure (truth {truth}):")
            print("      " + "  ".join(f"{k} {_fmt(v)}" for k, v in test.items()))
            result["closure"] = {"truth": truth, "fitted": test}

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2) + "\n")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
