"""The one pass over the held-out years.

Runs the surviving signal, the asset-level volatility-state surface traded
slowly, walk-forward through the full panel including the years from
``HOLDOUT_START`` that no development run touched, and reports the held-out
period on its own next to the development sub-periods, for the wide panel
and for the currencies alone.

Run once: python examples/final_holdout.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel import HOLDOUT_START, load_wide_panel  # noqa: E402
from _trading import evaluate, ic, print_table, stacked_rows, strip  # noqa: E402
from factor_curves import DROP  # noqa: E402

from badlm import BADLM  # noqa: E402
from badlm.backtest import deflated_sharpe, subperiod_sharpe, walk_forward  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
resdir = os.path.join(here, "..", "paper", "results")

if __name__ == "__main__":
    panel = load_wide_panel()
    close, kind = panel["close"], panel["kind"]
    close = close.drop(columns=[c for c in DROP if c in close])
    start = int(close.index.get_indexer([pd.Timestamp("1996-01-01")], method="bfill")[0])
    hold = int(close.index.get_indexer([HOLDOUT_START], method="bfill")[0])
    breaks = [int(close.index.get_indexer([pd.Timestamp(y)], method="bfill")[0]) for y in ("2004-01-01", "2012-01-01")] + [hold]
    pri = Priors(tau_level_scale=0.001, temper=10)
    out = {}
    for label, use in (("wide panel", ("equity_us", "equity_uk", "index", "fx", "commodity")), ("currencies", ("fx",))):
        df = stacked_rows(close, kinds=kind, use=use)
        med = df.groupby("asset")["vol"].transform(lambda s: s.expanding(min_periods=60).median())
        df["volstate"] = np.log(df["vol"] / med).fillna(0.0).clip(-2, 2)
        wf = walk_forward(BADLM(nlag=250, k=25, modifier_k=5, priors=pri), df["x"].to_numpy(), df["y"].to_numpy(),
                          df["time"].to_numpy(), groups=df["asset"].to_numpy(), start=start, step=63,
                          modifier=df["volstate"].to_numpy())
        rows = df.iloc[wf.row].reset_index(drop=True)
        rows["pred"], rows["pvar"] = wf.mean, wf.var
        res = evaluate("volatility surface", (rows["pred"] / rows["pvar"]).to_numpy(), rows, rules=(("gp", {"speed": 0.05}),))
        r = res[0]
        pnl = r["pnl"] - 0.0005 * r["traded"]
        sub = subperiod_sharpe(pnl, r["time"], breaks)
        held = rows["time"] >= hold
        dev = ~held
        print(f"\n{label}: {df['asset'].nunique()} assets; IC development {ic(rows.loc[dev, 'pred'], rows.loc[dev, 'y']):.4f}, "
              f"held out {ic(rows.loc[held, 'pred'], rows.loc[held, 'y']):.4f}")
        print("  net Sharpe (5 bp, slow) by period:", ", ".join(f"{close.index[int(lo)].year}-{close.index[min(int(hi), len(close) - 1)].year}: {s:.2f} ({n} days)"
                                                           for lo, hi, s, n in sub))
        hp = pnl[r["time"] >= hold]
        dsr, sr0 = deflated_sharpe(hp, 12)
        print(f"  held-out years {HOLDOUT_START.year}-{close.index[int(rows['time'].max())].year}: gross "
              f"{np.mean(r['pnl'][r['time'] >= hold]) / np.std(r['pnl'][r['time'] >= hold], ddof=1) * np.sqrt(252):.2f}, "
              f"net {hp.mean() / hp.std(ddof=1) * np.sqrt(252):.2f}; deflated Sharpe probability {dsr:.3f} against twelve trials")
        print_table(res, f"{label}, whole period 1996-{close.index[int(rows['time'].max())].year}")
        out[label] = {"results": strip(res), "subperiods": [(int(lo), int(hi), float(s), int(n)) for lo, hi, s, n in sub],
                      "ic_dev": ic(rows.loc[dev, "pred"], rows.loc[dev, "y"]), "ic_held": ic(rows.loc[held, "pred"], rows.loc[held, "y"]),
                      "held_net_sharpe": float(hp.mean() / hp.std(ddof=1) * np.sqrt(252)), "deflated": dsr}
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir, "final_holdout.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("written")
