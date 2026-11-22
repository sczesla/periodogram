"""Draws the figures in doc/images from tool/plot_data.json.

Run from the package root:

    dart run tool/plot_data.dart && python3 tool/plots.py
"""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

with open("tool/plot_data.json") as f:
    d = json.load(f)

t, y, err = np.array(d["t"]), np.array(d["y"]), np.array(d["err"])
res = np.array(d["residuals"])
freq = np.array(d["frequencies"])
style = dict(fmt="o", ms=3.5, color="#1f4e79", ecolor="#7f9fbf", elinewidth=1,
             capsize=0)
model_color = "#d1495b"
true_color = "#2a9d8f"
dpi = 160

# Figure 1: data with the best-fit sine (top), residuals after subtracting
# it (bottom); against time (left), and folded on the best period (top
# right) or on the period found in the residuals (bottom right)
def fold(x, period, t0):
    return ((x - t0) / period) % 1


def folded_curve(period, t0, model):
    ph = fold(tf, period, t0)
    order = np.argsort(ph)
    return ph[order], model[order]


period, t0 = d["bestPeriod"], d["bestT0"]
period2, t02 = d["secondPeriod"], d["secondT0"]
tf, mf = np.array(d["tFine"]), np.array(d["modelFine"])
mf2 = np.array(d["secondModelFine"])

fig = plt.figure(figsize=(9, 4.6), facecolor="white")
outer = fig.add_gridspec(1, 2, width_ratios=[1.6, 1], wspace=0.06)
left = outer[0].subgridspec(2, 1, height_ratios=[3, 1], hspace=0)
right = outer[1].subgridspec(2, 1, height_ratios=[3, 1], hspace=0)
ax1 = fig.add_subplot(left[0])
ax3 = fig.add_subplot(left[1], sharex=ax1)
ax2 = fig.add_subplot(right[0], sharey=ax1)
ax4 = fig.add_subplot(right[1], sharex=ax2, sharey=ax3)

ax1.errorbar(t, y, yerr=err, **style)
ax1.plot(tf, mf, color=model_color, lw=1.2)
ax1.set_ylabel("Value")
ax1.set_title("Data and best-fit sine", fontsize=10)
ax1.tick_params(labelbottom=False)

ax3.errorbar(t, res, yerr=err, **style)
ax3.axhline(0, color=model_color, lw=1)
ax3.set_ylabel("Residuals")
ax3.set_xlabel("Time [days]")
lim = 1.15 * np.max(np.abs(res) + err)
ax3.set_ylim(-lim, lim)

ax2.errorbar(fold(t, period, t0), y, yerr=err, **style)
ax2.plot(*folded_curve(period, t0, mf), color=model_color, lw=1.2)
ax2.set_title("Folded", fontsize=10)
ax2.text(0.97, 0.95, f"on {period:.2f} days", transform=ax2.transAxes,
         ha="right", va="top", fontsize=9)
ax2.tick_params(labelleft=False, labelbottom=False)

ax4.errorbar(fold(t, period2, t02), res, yerr=err, **style)
ax4.plot(*folded_curve(period2, t02, mf2), color=model_color, lw=1.2)
ax4.text(0.97, 0.92, f"on {period2:.2f} days", transform=ax4.transAxes,
         ha="right", va="top", fontsize=9)
ax4.set_xlabel("Phase")
ax4.tick_params(labelleft=False)
ax4.set_xlim(0, 1)
fig.savefig("doc/images/data.png", dpi=dpi, facecolor="white",
            bbox_inches="tight")


# Figure 2: periodogram, spectral window, and periodogram after subtracting
# the best-fit sine
def periodogram_panel(ax, power, levels, label_true):
    for k, f in enumerate(d["trueFrequencies"]):
        ax.axvline(f, color=true_color, lw=1, ls="--", alpha=0.8,
                   label="true frequencies" if label_true and k == 0 else None)
    ax.plot(freq, power, color="#1f4e79", lw=1.2, label="GLS")
    ax.axhline(levels[0], color="0.35", ls=":", lw=1, label="FAP = 1 %")
    ax.axhline(levels[1], color="0.35", ls="-.", lw=1, label="FAP = 0.1 %")
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Power")


def peak_near(power, f, width=0.05):
    """Frequency and power of the highest point within width of f."""
    sel = np.abs(freq - f) < width
    k = np.argmax(np.where(sel, power, -1))
    return freq[k], power[k]


line = dict(arrowstyle="-", color="0.4", lw=0.8)
fig, axs = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True,
                        facecolor="white")

ax = axs[0]
gls = np.array(d["gls"])
periodogram_panel(ax, gls, d["glsLevels"], True)
f1 = 1 / period
ax.annotate(f"{period:.2f} days", peak_near(gls, f1), xytext=(f1 + 0.1, 1.0),
            arrowprops=line, fontsize=9)
fa, pa = peak_near(gls, 1 - f1)
ax.annotate("alias at 1/day − f", (fa, pa), xytext=(fa - 0.11, 0.95),
            ha="right", arrowprops=line, fontsize=9)
ax.set_title("Periodogram", fontsize=10)
ax.legend(loc="upper right", fontsize=8, frameon=False)

ax = axs[1]
ax.plot(freq, d["window"], color="#6a4c93", lw=1.2)
ax.set_ylim(0, 1.05)
ax.set_ylabel("Window")
ax.set_title("Spectral window: peak at 1/day from nightly sampling",
             fontsize=10)

ax = axs[2]
rgls = np.array(d["residualGls"])
periodogram_panel(ax, rgls, d["residualLevels"], False)
p2 = d["secondPeriod"]
ax.annotate(f"{p2:.2f} days", peak_near(rgls, 1 / p2),
            xytext=(1 / p2 + 0.1, 0.75), arrowprops=line, fontsize=9)
ax.set_xlabel("Frequency [1/day]")
ax.set_title(f"After subtracting the {period:.2f}-day sine", fontsize=10)
ax.set_xlim(0, freq.max())
fig.tight_layout()
fig.savefig("doc/images/periodogram.png", dpi=dpi, facecolor="white")
print("Written doc/images/data.png and doc/images/periodogram.png")
