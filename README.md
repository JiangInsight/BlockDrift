# BlockDrift

*Exact Top-k Changes Between Historical Windows from Paired Prefix Bounds*

The two histories below have identical endpoint count distributions. In the first, both keys have zero change between windows. In the second, A loses ten and B gains ten. BlockDrift keeps differences between prefixes of the **same key** to bound its signed change.

<img src="figures/window-changes.png" width="590" alt="Two histories with identical endpoint marginals yield changes of (0,0) and (-10,+10).">

#### Which windows favor the index?

![Fastest independent control for each window family and log. Each cell gives the control's service time divided by full-word BlockDrift's time.](figures/window-families.png)

In Figure 3, values below 1 favor the named control. On nine snapshots, full-word BlockDrift won seven ordinary-window comparisons against the fastest independent control. Fixed pairing retained all seven wins.

**Try the ranking.** Run `build.cmd` followed by `python run_example.py` on x64 Windows with Visual Studio C++ Build Tools and Python 3.10+. The example checks 1,323 signed top-k queries against direct counts.
