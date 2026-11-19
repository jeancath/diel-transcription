"""Figures. Three panels carry the result:

fig1  the reproducibility gap: ortholog-specific phase agreement, within one
      strain across two studies vs between lineages, against the null that
      preserves each dataset's own phase distribution
fig2  the distribution of |phase difference| behind that statistic
fig3  each lineage keeps its own rush hour - the marginal phase distributions
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT, phase_diff, load
from robustness import rank_table
from datasets import BY_ID

FIG = Path(__file__).resolve().parents[1] / "figures"
SURFACE = "#fcfcfb"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, VIOLET = ("#2a78d6", "#eb6834", "#1baf7a",
                                               "#eda100", "#e87ba4", "#6250d6")
INK, INK2, MUTED = "#1a1a19", "#4a4a47", "#8a8a85"
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE, "font.size": 9,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False})
LAB = {d.did: d.lineage for d in BY_ID.values()}
LAB["SE7942_vijayan"] = "S. elongatus (study 2)"


def fig1():
    """Effect size, not p value: how far below the marginal-preserving null the
    mean |phase difference| falls, for each pair and each selection rule."""
    t = pd.read_csv(OUT / "n2_threshold_sensitivity.csv")
    order = ["all orthologs", "either q<0.05", "both q<0.05", "top 200 combined"]
    t = t[t.rule.isin(order)]
    pairs = (t[t.comparison == "cross-lineage"].groupby(["a", "b"])["effect_pct"]
             .median().sort_values())
    win = ("SE7942_ito_WT", "SE7942_vijayan")
    rows = list(pairs.index) + [win]
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ax.axvline(0, color=MUTED, lw=1, zorder=1)
    marks = dict(zip(order, ["o", "s", "^", "D"]))
    for i, (a, b) in enumerate(rows):
        sub = t[(t.a == a) & (t.b == b)]
        c = BLUE if (a, b) == win else ORANGE
        lo, hi = sub.effect_pct.min(), sub.effect_pct.max()
        ax.plot([lo, hi], [i, i], color=c, lw=2, alpha=.35, zorder=2,
                solid_capstyle="round")
        for r in sub.itertuples():
            ax.plot(r.effect_pct, i, marks[r.rule], ms=6, color=c, mec=SURFACE,
                    mew=1.2, zorder=3)
        ax.annotate(f"{sub.effect_pct.median() + 0.0:+.0f}%".replace("+0%", "0%"), (hi, i), xytext=(9, 0),
                    textcoords="offset points", va="center", fontsize=8,
                    color=INK if (a, b) == win else INK2,
                    fontweight="bold" if (a, b) == win else "normal")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{LAB.get(a,a)}  vs  {LAB.get(b,b)}" for a, b in rows], fontsize=8)
    for lab, (a, b) in zip(ax.get_yticklabels(), rows):
        w = (a, b) == win
        lab.set_color(INK if w else INK2)
        if w:
            lab.set_fontweight("bold")
    ax.set_xlim(-12, 100)
    ax.set_xlabel("ortholog-specific phase conservation\n"
                  "(% reduction in mean |phase difference| vs the marginal-preserving null)")
    ax.set_title("The conserved component between lineages is a tenth of\n"
                 "what two studies of one strain recover",
                 color=INK, fontsize=11.5, loc="left", pad=12)
    handles = [plt.Line2D([], [], marker=m, ls="none", color=MUTED, ms=6, label=r)
               for r, m in marks.items()]
    ax.legend(handles=handles, frameon=False, fontsize=7.5, loc="lower right",
              title="ortholog selection rule", title_fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_effect_size.png", dpi=220)
    plt.close(fig)


def fig2():
    """Distribution of |phase difference| behind the effect sizes, using every
    one-to-one ortholog (no rhythm threshold) with the best global shift applied."""
    from null_comparison import prepare, phases, best_shift
    rng = np.random.default_rng(1)
    cross = [("SE7942_ito_WT", "SYN6803"), ("SE7942_ito_WT", "TERY"),
             ("SE7942_ito_WT", "MED4"), ("SE7942_ito_WT", "NPUN"),
             ("SYN6803", "TERY"), ("SYN6803", "MED4"), ("SYN6803", "NPUN"),
             ("TERY", "MED4"), ("TERY", "NPUN"), ("MED4", "NPUN")]

    def collect(pairs):
        d, nul = [], []
        for a, b in pairs:
            ea, ma = prepare(a); eb, mb = prepare(b)
            sh = ea.index.intersection(eb.index)
            if len(sh) < 50:
                continue
            pa, _ = phases(ea.loc[sh], ma)
            pb, _ = phases(eb.loc[sh], mb)
            x, y = pa.loc[sh].values, pb.loc[sh].values
            _, s = best_shift(x, y)
            d.append(np.abs(phase_diff(x, y + s)))
            yp = rng.permutation(y)
            _, s2 = best_shift(x, yp)
            nul.append(np.abs(phase_diff(x, yp + s2)))
        return np.concatenate(d), np.concatenate(nul)

    c, cn = collect(cross)
    w, _ = collect([("SE7942_ito_WT", "SE7942_vijayan")])
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    grid = np.linspace(0, 12, 300)

    def ecdf(v, **kw):
        ax.plot(grid, [(v <= g).mean() for g in grid], **kw)

    ecdf(w, color=BLUE, lw=2.2, label=f"within strain, 7942 x2 (n={len(w)})")
    ecdf(c, color=ORANGE, lw=2.2, label=f"between lineages, 10 pairs (n={len(c)})")
    ecdf(cn, color=MUTED, lw=1.6, ls=(0, (4, 3)),
         label="null: ortholog identity shuffled")
    ax.axvline(3, color=MUTED, lw=.8, alpha=.6)
    ax.text(3.15, .06, "3 h", color=MUTED, fontsize=7.5)
    ax.set_xlabel("|difference in peak time| after the best global shift  (h)")
    ax.set_ylabel("fraction of one-to-one orthologs")
    ax.set_xlim(0, 12); ax.set_ylim(0, 1)
    ax.set_title("Every one-to-one ortholog, no rhythm threshold applied",
                 color=INK, fontsize=11, loc="left", pad=12)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_phase_difference_ecdf.png", dpi=220)
    plt.close(fig)


def fig3():
    df = load()
    order = ["SE7942_ito_WT", "SE7942_vijayan", "SYN6803", "TERY", "MED4", "NPUN"]
    cols = [BLUE, BLUE, ORANGE, AQUA, VIOLET, YELLOW]
    fig, axes = plt.subplots(2, 3, figsize=(8.4, 4.6), sharex=True, sharey=True)
    bins = np.arange(0, 25, 2)
    for ax, d, c in zip(axes.ravel(), order, cols):
        g = df[(df.dataset == d) & (df.q < 0.05)]
        h, _ = np.histogram(g.phase_h, bins=bins)
        frac = h / h.sum()
        ax.bar(bins[:-1] + 1, frac, width=1.6, color=c, linewidth=0)
        ax.set_title(f"{LAB[d]}   n={len(g)}", fontsize=8.5, color=INK, loc="left")
        ax.set_xticks([0, 6, 12, 18, 24])
        ax.axvspan(12, 24, color="#000000", alpha=.045, lw=0,
                   zorder=0) if BY_ID.get(d) and BY_ID[d].light == "LD" else None
        ax.tick_params(labelsize=7.5)
    for ax in axes[1]:
        ax.set_xlabel("peak time (h)", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("fraction of rhythmic\ntranscripts", fontsize=8)
    fig.suptitle("Each lineage concentrates its rhythmic transcripts at a different time of day",
                 color=INK, fontsize=11, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG / "fig3_phase_marginals.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    fig1(); fig2(); fig3()
    print("wrote", *[p.name for p in sorted(FIG.glob("*.png"))])
