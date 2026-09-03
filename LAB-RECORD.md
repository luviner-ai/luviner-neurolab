# Lab record

A curated export of the laboratory notebook kept while the work in
`paper/` was done. Every section below was written in a private
repository *before* the run it describes, or immediately after it, and has
not been edited since — corrections are appended, never applied in place,
which is why several sections end by contradicting their own opening.

Each section carries the short hash of the commit that introduced it in
the private repository. Those commits are the timestamp; this file cites
them, it does not replace them. Commit hashes refer to the private
working record in which each entry was committed before the corresponding
run; that record is available to reviewers on request.

Registered predictions and their falsifiers appear before results, in the
same words used to score them. Negative results are kept. Work outside the
biophysics line — engine, deployment, product — is not part of this
release and is not included here; nor is the project's task board, its
costs, or the identifiers of the sessions that did the work.

---


<!-- private working record: commit 50b3891, 2026-08-30 -->

### The neutral point of a plasticity rule is a ratio, and the network sets it

The registered falsifier for N6 was: *the learning must vanish when the
window integral is set to zero*. That experiment cannot be run, and the
reason is worth more than the experiment would have been.

`window_integral` — `A_plus*tau_plus - A_minus*tau_minus` — is the net
change per pair **when the lags are spread evenly over the window**. A
stimulated network does not spread them evenly. It concentrates them near
zero, where the causal side of the window is the taller of the two, so a
rule with no area still potentiates. Measured: at integral 0 the expected
change is **+0.00009 per pair**, not zero, and it has the same sign under
every stimulation protocol tried.

The correct statement is scale-free and is a ratio. The net change is
`A_plus*P - A_minus*M`, where P and M are the means, over the lags the
network actually produces, of the two exponentials of the window. The
rule is neutral when

    A_minus / A_plus  =  P / M

The window integral is the special case `P/M = tau_plus/tau_minus`, which
is what "lags uniform over the window" means.

**The test.** Predict P/M from a frozen dish's lags — frozen, so the
prediction cannot depend on the rule it is about to be tested against.
Then run a ladder in `A_minus/A_plus` on plastic dishes and find where
the realised mean weight change crosses zero. Five seeds, paired, nothing
within 0.2% of a bound and weights moving 11-19% of their initial value,
so the ladder is reading drift and not clipping.

| A_minus/A_plus | 0.30 | 0.42 | **0.50** | 0.55 | 0.62 | 0.75 | 0.95 |
|---|---|---|---|---|---|---|---|
| mean dw, x10^-3 | +1.09 | +0.66 | **+0.35** | +0.17 | -0.18 | -0.86 | -1.90 |

The bold column is the zero-area rule: at 0.50 the window integral is
-0.0002, which is zero to the precision of the ladder.

                      predicted from the lags   measured crossing
    seed 0                         0.5655             0.5785
    seed 1                         0.6316             0.5699
    seed 2                         0.5975             0.5724
    seed 3                         0.5613             0.5791
    seed 4                         0.6694             0.6908
    ------------------------------------------------------------
    mean                    0.6051 +/- 0.0458   0.5981 +/- 0.0520

The paired difference is **-0.0069 +/- 0.0359, t = -0.43**: the
prediction has no bias. It tracks the individual network, **r = +0.738**
across seeds — the seed with the highest P/M is the seed with the highest
crossing. And the zero-area ratio is rejected: the measured neutral point
sits **+0.0996 +/- 0.0520 above it, t = +4.29, 5/5 seeds**. At the ratio
that makes the window integral vanish the network potentiates, at
+0.00035 per synapse — a number that would have been read as "the
plasticity was switched off".

**The cutoff is not a matter of taste, and this nearly went in wrong.**
P and M are normalised by the number of lags collected, so a collection
window that truncates the wider (depression) side inflates P/M. Across
cutoffs of 40, 60, 100, 150, 250, 400 and 600 ms the ratio reads 0.804,
0.661, 0.608, 0.592, 0.586, 0.586, 0.586. It is flat past 250 and wrong
by a factor of 3.5 in the neutral integral at 40. The first version of
this measurement used 60 ms, chosen because it looked generous next to a
33.7 ms time constant, and would have published -0.066 for a quantity
that is -0.035. Convergence curve in `experiments/dish/lag_cutoff.py`;
the default cutoff is now 250 ms and a test holds it there.

**What is and is not new here.** The mathematics is classical: that a
pair rule's mean effect is its window integrated against the pre/post
correlation function is Kempter-Gerstner-van Hemmen 1999, and
Song-Miller-Abbott chose `A_minus*tau_minus > A_plus*tau_plus` in 2000
precisely to be net-depressing for uncorrelated input. Three things are
not in that literature:

1. **the ablation consequence.** Every claim of the form "the effect is
   STDP, because it disappears when the window integral is zeroed" is
   using a control that does not control. The mechanism ablation has to
   be defined on the observed lag distribution, and
   `rule_with_expected_change(0, lags)` is the instrument that does it.
2. **the measured value** of P/M for a stimulated conductance-based E/I
   network, with the convergence requirement that makes it measurable.
3. **its stability.** This is where a claim died. Read at a 60 ms
   cutoff, P/M looked strongly protocol-dependent — 0.586 under
   multi-site stimulation against 0.816 under the feedback trains, which
   would have meant that "is this rule potentiating here?" has no answer
   belonging to the network. Read past the knee, the same four protocols
   give 0.598, 0.590, 0.606 and 0.613: a spread of **4%**. The
   protocol-dependence was the truncation, not the network. What survives
   is the better result — P/M is a *constant of the network*, measurable
   once and reusable, sitting 20% above the value that zeroes the area.

That is what the registered prediction should have said, and it is why
the falsifier as written could not have decided anything.

---


<!-- private working record: commit 52a24ed, 2026-09-01 -->

## Run: stg-P/M — the bursts do read the weight order, and not by the mechanism I registered

*Recorded 2026-09-01. `experiments/dish/stg_pm.py`, registered
before any number existed. Ten wirings, three batches, ~37 minutes.
Spikes, weights, classes and burst windows persisted; the follow-up
analyses below cost no compute.*

`## Run: N8` measured `P/M` by weight class on 18–46 ms events and found
bursts to be "approximately a no-op": `+0.0041 ± 0.0432`, a tenth of the
between-seed spread, both classes on the same side in 3 of 5 seeds.
`## Run: dish-slow2c` built bursts of 22–25 spikes over 444–492 ms. This
is the same instrument on those.

| seed | bursts | dur | ρ | P/M strong | P/M weak | diff |
|---|---|---|---|---|---|---|
| 0 | 17 | 313 | +0.068 | 0.7180 | 0.5795 | **+0.1386** |
| 1 | 22 | 448 | +0.122 | 0.8879 | 0.8601 | +0.0278 |
| 2 | 22 | 445 | −0.059 | 0.8743 | 0.8519 | +0.0225 |
| 3 | 22 | 449 | −0.029 | 0.8615 | 0.8451 | +0.0163 |
| 4 | 22 | 449 | −0.183 | 0.8709 | 0.8490 | +0.0219 |
| 5 | 22 | 448 | −0.031 | 0.8698 | 0.8567 | +0.0131 |
| 6 | 22 | 446 | +0.142 | 0.8775 | 0.8422 | +0.0354 |
| 7 | 22 | 445 | +0.184 | 0.8607 | 0.8396 | +0.0212 |
| 8 | 17 | 388 | +0.050 | 0.6936 | 0.5778 | **+0.1158** |
| 9 | 22 | 445 | −0.406 | 0.8812 | 0.8732 | +0.0080 |

### (a) The mechanism I registered: REFUTED

> *Spearman ρ < 0 between a cell's median within-burst first-spike
> latency and its total incoming weight, in ≥ 8/10 wirings.*

**Negative in 5 of 10. Mean −0.014 ± 0.176, t = −0.26.** Centred on
zero, not on anything. **Within-burst firing order, as first-spike
latency, does not track synaptic weight at all.**

### (b) The outcome: one clause met, one not

> *strong-minus-weak `P/M` positive in ≥ 8/10 wirings, **and** the effect
> exceeding the between-wiring spread.*

    positive in 10/10          t = +2.90        SIGN CLAUSE MET
    mean +0.0421, sd 0.0458    effect/spread 0.92   SPREAD CLAUSE NOT MET

**I wrote both clauses and one of them fails, so (b) is not met as
registered.** The direction is unambiguous — ten wirings out of ten, a
sign test at p ≈ 0.001, and every wiring's effect is larger than N8's
mean. The magnitude varies between wirings by slightly more than its own
mean, which is exactly what the second clause was built to catch, and it
catches it. Dropping seed 0 does not rescue it (0.96), so this is not an
outlier problem.

**An observation about the spread, which cannot and does not rescue the
clause.** The ten wirings are two regimes, separable by **burst count**
— an independent variable, not the outcome:

| | n | diff | effect/spread | classes vs neutral 0.8299 |
|---|---|---|---|---|
| 17-burst | 2 | +0.1272 ± 0.0161 | 7.9 | both **below** (0.706/0.579) |
| 22-burst | 8 | +0.0208 ± 0.0085 | 2.4 | both **above** (0.873/0.852) |

Each regime passes the spread clause comfortably; pooling two regimes
inflates the spread and fails it. **The split was noticed after seeing
the result and was not registered, so the registered clause stands as
failed.** It is recorded because it is the honest description of the
data, not because it changes the score.

### What carries the sorting — the pre-authorised re-analyses, read only after scoring

|  | strong | weak | difference | same sign |
|---|---|---|---|---|
| spikes per burst | 20.33 | 20.50 | −0.17 | 2/10 |
| burst participation | 0.894 | 0.901 | −0.007 | 2/10 |
| **causal pair fraction** | **0.5147** | **0.5125** | **+0.0022** | **10/10** |

Not spike count, not participation — both are coins. The carrier is a
**0.22-percentage-point excess of causal over anticausal pairs** on
strong synapses relative to weak ones, and it is the same sign in every
one of the ten wirings.

**That is why (a) failed while (b)'s direction held.** (a) asked about
firing order at the scale of a burst — which cell starts first — and at
that scale there is nothing. The signal lives in the *pairwise ordering
statistics inside* the burst: across roughly 300,000 pre-post pairs per [CORRECTION APPENDED 2026-09-03: this figure is the per-WIRING total (~3×10⁶ over ~340 synapses; recomputed from persisted PM-4 spikes, median 8,960 pairs per synapse at |lag| ≤ 250 ms inside burst windows). The "per synapse" reading propagated into the preprint abstract and §5 and was caught by the figure-6 redraw; corrected everywhere, original text left as written.] per
synapse, strong synapses see the presynaptic spike first slightly more
often. The exponential STDP kernels integrate that sub-percent bias into
a `P/M` difference of a few percent.

**A tiny, utterly consistent asymmetry, amplified by the rule.** It
rhymes with `## Run: phase-link`, where a per-step drift six orders below
the weight scale moved the operating point by 2% through a `tau/dt` of
4000. This project keeps finding that the interesting quantity is a small
bias with a large amplifier, not a large effect.

### Verdict

- **(a) refuted.** The mechanism I proposed is not the mechanism. I gave
  (b) a 45% prior and (a) no explicit prior; had (a) been the only
  measurement, the correct conclusion would have been "no order, so no
  sorting", and it would have been wrong.
- **(b) not met as registered** — sign yes, spread no. The consolidation
  claim is therefore **not** established at the bar I set. What is
  established is weaker and still substantial: the direction is
  consistent in 10/10 wirings and an order of magnitude above N8's, and
  the mechanism carrying it is identified.
- **N8's null is explained rather than merely contradicted.** N8's events
  had one or two spikes; the carrier here is a bias in pairwise ordering
  that needs thousands of pairs per synapse to express. With 18–46 ms
  events there were not enough pairs for a 0.2% asymmetry to become
  anything. The bursts were not too brief to *have* an order — they were
  too brief to *accumulate* one.
- **Where this leaves consolidation.** In a preparation built to give it
  every advantage, spontaneous activity sorts the weight order in the
  right direction every time, by an amount that varies more than its own
  mean across wirings. That is neither the clean positive nor the clean
  null, and saying so is the result.

### What transfers

**The double-barrelled registration earned its keep in the direction it
was designed to.** (a) could never rescue (b) — that was the binding —
and in the event (a) *failed* while (b)'s direction held, which is the
configuration where a single-barrelled registration would have been
actively misleading. Measuring the mechanism alongside the outcome told
us the outcome is real and our story about it was wrong. That is worth
more than either number alone.

**A criterion with two clauses will eventually split.** I wrote "positive
in ≥ 8/10 **and** exceeding the spread" precisely so a consistent-but-
variable effect could not be called a win, and then produced exactly
that. The discipline is to report the split as a split rather than to
lead with the clause that passed.

**On persisting the wiring.** All three follow-ups needed `pre`/`post`,
which I did not persist — recoverable from the seed, and verified so
because `w_in` was computed inside the run and reproduced to `atol=0`.
That witness was luck. Persist the wiring, not a way to rebuild it.

---


<!-- private working record: commit ea95d3d, 2026-09-01 -->

## Run: PM-4

Direction-only replication on ten wirings that did not exist when the
predictions were written. Registered and **pushed** as `a6b18f9` before
the first simulation step; results below.

    seed  bursts   P/M diff  causal diff
      10      22    +0.0277     +0.00128
      11      22    +0.0389     +0.00202
      12      22    +0.0244     +0.00113
      13      22    +0.0317     +0.00148
      14      22    +0.0101     +0.00050
      15      22    +0.0257     +0.00114
      16      22    +0.0658     +0.00264
      17      30    +0.2296     +0.01522
      18      10    +0.1354     +0.00744
      19      22    +0.0396     +0.00176

    (1) P/M diff positive          10/10   MET  (>= 9)
    (2) causal excess positive     10/10   MET  (>= 9)
    (3) Spearman rho               +0.988  MET  (>= 0.60)

All three met. **But (3) is worth less than its number, and the reason
has to be on the record rather than in the pleased silence after a
pass.** `P/M` and the causal fraction are computed from the same spike
pairs — one exponentially weighted, one counted. Their ratio across the
ten wirings is 20.8 +/- 2.7, a coefficient of variation of 0.13: the
`P/M` difference is very nearly the pair-count difference times a
constant. So rho = 0.988 is close to arithmetic. What it genuinely
establishes is narrower and still useful: the effect is carried by
**how many** causal pairs there are, not by their fine timing within the
STDP kernel. The exponential weighting adds nothing. Predictions (1) and
(2) are not independent of each other either, for the same reason.

So PM-4 is honestly **one** out-of-sample result, not three: the sign of
the strong-minus-weak difference replicated 10/10 on fresh wirings
(p ~ 0.001 against a fair-coin null), and it is a counting effect. That
is exactly the accumulation reading, and it now survives a test that
could have killed it. It is not three independent confirmations, and
`## Run: stg-P/M`'s claim should be read at this reduced strength too.

Cost: 15.3 min for four (estimated 14.8), then two workers x 11.1 min
(estimated 11.1) — 3.7 s/s holds to 3%.

### Descriptive, and deliberately not scored

The magnitude spread is still the size of the mean: +0.0629 +/- 0.0683,
effect/spread 0.92. PM-3 said this would happen and PM-4 was registered
knowing it would; it decides nothing here.

One thing in the descriptive column is more interesting than the scored
column, so it needs the strongest possible label: **this is post-hoc,
found after unblinding, and predicts nothing yet.** Eight of the ten
wirings burst 22 times; the two that did not (seed 17 at 30, seed 18 at
10) have the two largest `P/M` differences, and every off-regime wiring
exceeds every on-regime one (0.135-0.230 vs 0.010-0.066). Both a faster
and a slower network give a bigger effect, so it is not monotone in
rate — rho(magnitude, burst count) is only +0.24, while
rho(magnitude, burst duration) is -0.47. The candidate reading is that
shorter bursts sort the order harder, which would make PM-3's
"intrinsic, unstratifiable" variance partly explainable after all, by a
variable measured *during* the run rather than assignable before it.

**With n = 2 off-regime wirings this is an anecdote with a mechanism
attached.** It is written down as the seed of the next registration, not
as a finding, and the next run must fix burst duration as a predictor
*before* seeing which way it goes.

## Run: SM-2 (scoping only — no simulation run, no registration yet)

Does contingency appear on a substrate that demonstrably sorts weights?
The cortical dish failed this contrast while silent; Sinapayen's
spontaneously active model passes it. The STG dish is spontaneously
active *and* now known to sort. This is the scoping the peer asked for:
protocol, cost, and a registration proposal to be approved before
anything runs.

### The blocker, found before spending anything

Shahaf-Marom's criterion is "did the target respond within a window."
That works because their dish was near-silent between stimuli, so a
response carried information. Measured on PM-4's own spikes, this
preparation fires with duty cycle 0.32:

    chance a target 'responds' within  25 ms of an arbitrary probe: 0.26
    chance within  50 ms: 0.31      within 100 ms: 0.34
    within 250 ms: 0.43

**A third of probes succeed with no network doing anything.** Inheriting
the cortical criterion here would have produced a learning curve made of
burst phase. That is the same class of error as the inherited `I_M` and
`I_KCa` defaults, caught the same way: by measuring first.

The fix is to probe phase-locked to the inter-burst silence (median
935 ms, 10th percentile 780 ms), at a fixed delay `D` after burst
offset:

    D ms     W=25    W=50   W=100   W=200
      200   0.001   0.004   0.004   0.004
      300   0.002   0.003   0.003   0.004
      400   0.001   0.001   0.003   0.008
      800   0.028   0.029   0.034   0.857   <- next burst arrives

`D` in 200-400 ms gives a chance floor of 0.1-0.8 per cent, recovering
the quiet-dish condition by phase rather than by silence. `D = 300 ms`
has the widest margin on both sides. Cycles whose silence is shorter
than `D + W + 100 ms` must be skipped, not squeezed — the minimum
observed gap is negative (merged windows).

### Protocol

Stimulate two E cells with a brief current pulse at `D = 300 ms` after
each burst offset; the target is a *distant* E cell (no direct synapse
from either stimulated cell), fixed per wiring before the run. Three
arms, identical stimulus statistics:

- **contingent** — the pulse train stops when the target fires within
  the criterion window; the measure is stimuli-to-criterion over trials.
- **yoked** — stop times replayed from a *different seed's* contingent
  record, so the stimulation is statistically matched and causally
  unrelated. Different seed, per this line's standing rule.
- **frozen** — plasticity off, contingent schedule replayed. This sits
  beside any responsiveness claim: it separates learning from adaptation
  and from the dish's own drift, both of which can produce a curve.

`STGDish.run(duration, learn, extra)` already accepts an additive
per-cell current, so the stimulation surface exists. Contingent stopping
needs `extra` to accept a callable of `t` — a small change in
`stg_dish.py`, which this lane owns. At 4.8M calls per 240 s run the
Python overhead is ~2 s against ~15 min: irrelevant.

### Cost, at the measured 3.7 s wall per biological second

One burst cycle is 1.39 s, so probes arrive at ~43 per biological
minute. A 240 s run gives ~172 trials and costs **14.8 min** — one run
is exactly one legal unit under the 15-minute line, and the unit cannot
be made bigger, only more numerous.

    calibration              3 wirings x 60 s        11.1 min
    3 arms x 3 seeds         9 units x 14.8 min       2.2 h  (1.1 h on 2 workers)
    3 arms x 5 seeds        15 units x 14.8 min       3.7 h  (1.9 h on 2 workers)

**That is a real cost and I am not going to start it on my own
authority.** It is one to two orders more wall clock than anything this
lane has spent per question so far.

### Two gates before the expensive part

The calibration run (11.1 min) must clear both, or SM-2 does not proceed
in this form:

1. **Response latency is measured, not inherited.** The criterion window
   `W` is fixed from the observed distribution of target-firing latency
   after a test pulse — a number this preparation has never been asked
   for. `W` is *not* chosen from the cortical dish and not from the
   numbers above, which are chance floors, not latencies.
2. **Baseline responsiveness is intermediate.** If the pulse already
   evokes the target on ~every trial, or on ~none, there is nothing to
   learn in either direction and no contrast can show anything. This is
   a go/no-go, and it is the cheapest way SM-2 can die.

### Registration proposal, for the orchestrator's pass

Predictions, on `N` fresh wirings, one file, wiring persisted, per-seed
rows:

- **P1** stimuli-to-criterion falls over trials in the contingent arm,
  in >= N-1 of N wirings.
- **P2** the contingent fall exceeds the yoked fall in >= N-1 of N.
  *This is the load-bearing one* — P1 alone is satisfied by adaptation.
- **P3** the frozen arm shows no fall, or one significantly smaller than
  contingent, in every wiring.

Falsifiers, stated so they cannot be renegotiated later:

- yoked falls as fast as contingent -> **no contingency**; the effect is
  stimulation, not its timing.
- frozen falls as fast as contingent -> the curve is dynamics, not
  plasticity, and the weight-sorting result buys nothing here.
- baseline response probability outside 0.15-0.85 at calibration ->
  **SM-2 does not run** in this form.
- fewer than 120 usable trials per run after the short-gap skip ->
  underpowered, report and stop rather than extend the runs past the
  line.

Magnitude of any fall is descriptive. PM-4's lesson applies directly: if
these measures turn out to share their inputs, they are one result
reported as one, not three.

Open and not yet decided: `N`, and whether to buy the 3-seed or 5-seed
version. My recommendation is calibration first, then the 3-seed
contrast as a pilot whose job is to measure whether 172 trials is even
enough to move the curve — not to test the hypothesis. Extending to five
seeds is a separate decision made after seeing that, and the founder
should see the hour count before it is spent.

## Run: SM-2 calibration — both gates fail, SM-2 does not run in this form

Two passes, 22.7 min, plus a 40 s single-cell probe. **The pre-written
falsifier fires: baseline response probability is ~0.01, far outside
0.15-0.85.** The contrast is not worth an hour of anyone's machine, and
the founder's decision on it is moot.

Both passes printed `PASS`. Both were wrong, and the second was only
legible because it carried the control the first lacked.

### What the sham bought

Pass 1 measured response probability against a chance floor computed
from *other wirings'* unstimulated spikes, and reported 6/19 = 0.30 at
one amplitude. Pass 2 interleaved no-pulse trials at the same detected
phase in the same run: **sham responds 0/35 at every window**, so the
floor really is zero — and the same run still reported 0.27 at
0.1 uA/cm2. A clean floor was not enough to make the number true.

What made it legible was the other thing pass 2 recorded: whether the
pulse fired the stimulated cells at all.

    pulsed trials                       77
      pulse fired the stimulated cells  10   -> target responded  1/10
      pulse did NOT fire them           67   -> target responded  4/67

Four of the five "responses" occurred on trials where the pulse did
nothing whatever. They are burst-phase artifacts, at 264-382 ms. **One**
trial is a real evoked response: seed 25, 0.1 uA/cm2, stimulated cells
firing 21 spikes each, target at **7.1 ms** — the right modality, a
plausible synaptic latency. So the true rate is 1/77, and the printed
0.27 was four non-events and one signal.

### The mechanism, and why it is structural rather than a bad parameter

The pulse fires the stimulated cells more often at the *weakest*
amplitude — 0.27 at 0.1 uA/cm2, 0.00 at 0.8, and 0/82 across 1.0-16.0
in pass 1. One isolated AB/PD cell, probed the same way, explains it:

     amp uA/cm2   spikes evoked
        0.0020 .. 0.4559        0
        0.6548                  7
        0.9404, 1.3505          0
        1.9394 .. 4.0000        2

Threshold sits near 0.65 uA/cm2 and depolarization block begins
essentially on top of it. The usable band is narrow, non-monotone, and
sensitive to where in the slow wave the pulse lands.

**That is the finding, and it is not a parameter problem.** The phase
that makes the criterion clean is the phase that makes the cell
unresponsive. Probing 300 ms after burst offset drove the chance floor
from 31% to zero — and it is also the interburst hyperpolarized trough,
where these cells are least excitable and where the current that fires
them is already the current that blocks them. Shahaf-Marom's cortical
dish had silence *without* refractoriness. An intrinsically bursting
preparation cannot separate the two: its silence IS its refractory
phase. That is the difference between the two preparations, and no
choice of amplitude, pulse width or window dissolves it.

### What survives

Pass 1 measured something real that this protocol was not asking for:
a weak pulse **advances the next burst**, from 935 ms after offset to
~650 ms, on 32% of trials against a 2.0% spontaneous rate at that
window. Burst-level, not spike-level — which is why it broke the
criterion — but causal and measurable.

So the honest options, none of them a patch to this design:

1. Probe near burst onset, where excitability is rising, and pay for it
   with a higher chance floor **controlled by the sham** rather than
   avoided by phase. The sham machinery now exists and works.
2. Stimulate far more than 2 of 48 cells.
3. **Change the response variable to burst phase.** Ask whether
   contingency can shift the phase of the next burst, rather than
   whether it can evoke a spike. This is the one that fits what the
   preparation can actually do, and pass 1 already shows the effect is
   there to be moved.

Option 3 is a different experiment and needs its own registration with
its own falsifiers. It is not a rescue of this one, and it must not be
written as though SM-2 merely needed tuning.

### Cost and record

Pass 1 11.3 min, pass 2 11.4 min (11.1 estimated each), single-cell
probe 40 s. Usable trials after the short-gap skip: 112/207 (54%), which
extrapolates to ~135 per 240 s run — above the 120 floor, but that floor
never became load-bearing. Nothing beyond the authorized calibration was
spent, and the 2.2-3.7 h contrast was not started.

## Run: SM-3 (scoping) — the premise does not hold; the prerequisite is the real unit

SM-3 was framed on a substrate that "demonstrably lets a weak pulse
advance the next burst, 32% against 2%". **Neither number survives.**
Both come from the free re-analysis below, no new simulation.

### The 32% was the artifact-contaminated measure

That figure is pass 1's *target-cell firing* rate — the measure pass 2
showed was four non-events and one signal. Pass 2 persisted its spikes,
so the burst onsets can be read directly instead:

    amp    n   median next-burst onset   advance rate (<=700 ms)
    0.00  34            936 ms                 0.03    (sham)
    0.10  13            931 ms                 0.31
    0.25  11            930 ms                 0.00
    0.40   7            934 ms                 0.00
    0.65  13            934 ms                 0.00
    0.80  20            934 ms                 0.10

**The median is unmoved at every amplitude**, 930-936 ms against a sham
of 936. Pulsed trials pooled: median 934, advance 0.08, n=72, against
sham 936 and 0.03. The pulse does not shift the phase, because 67 of 72
pulses never fired a cell.

### The 2% was one wiring diluted by nine zeros

    seed  17 -> 0.227      all nine others -> 0.000
    floor across wirings: 0.023 +/- 0.072, max 0.227

The pooled 2.0% is not a floor any wiring has. Seed 17 is the off-regime
wiring — shorter period, median gap 754 ms — so a fixed 700 ms threshold
catches its ordinary cycles. **An absolute phase criterion is invalid
across wirings.** Any criterion has to be a percentile of that wiring's
own sham distribution, measured in the same run. This is the inherited-
number lesson again, one level down: not inherited from another
preparation, inherited from a pooled average across wirings.

### What is real, on n = 5

Split the pulsed trials by whether the stimulus actually landed:

    pulse fired the cells   n= 5   median 1326 ms   advance 0.40   min 307
    pulse fired nothing     n=67   median  933 ms   advance 0.06   min 564
    sham                    n=34   median  936 ms   advance 0.03   min 605

When the pulse lands, the phase moves **hard**: two trials at 307 and
379 ms against a spontaneous minimum of 605 ms over 34 sham trials.
It also moves the other way — the landed median is 1326 ms, so some
cycles are delayed. A landed pulse resets the rhythm; the direction is
not yet predictable. That is a genuine, large effect resting on five
events, and it is the honest version of what SM-3 wanted to build on.

### The prerequisite is the same wall SM-2 hit

SM-2 died because the pulse does not reliably fire cells at the trough.
SM-3 dies of the identical cause one level up. **Both reduce to one
unanswered question: can the stimulus be made to land?** Nothing about
contingency — spike-level or phase-level — is testable until it does.

So the next unit is not SM-3's contrast. It is a stimulation-efficacy
characterization, and it is cheap:

    single cell, phase x amplitude x pulse-width grid   ~2 min
    best setting verified in one dish wiring, 60 s      ~3.7 min

Go/no-go: a setting that lands on **>= 80%** of trials at a phase whose
sham floor is near zero. If no such setting exists, the stimulation line
closes honestly and neither SM-2 nor SM-3 is a question this preparation
can be asked — which is itself worth knowing, and costs six minutes to
establish rather than hours.

### Design held in reserve, for if the prerequisite clears

Recorded now so it is not reinvented, and explicitly not authorized.

The response variable is the next burst's onset, expressed as a
percentile of that wiring's own sham distribution. Two stimulus sites
**A** and **B**, disjoint cell pairs, alternating pseudo-randomly across
cycles; only A is contingent — an advance on an A cycle pauses
stimulation for k cycles, B never gates anything. The scored quantity is
`p(advance | A) - p(advance | B)` across blocks.

This within-run difference is why yoked becomes optional rather than
mandatory: A and B receive equal stimulation in the same dish, same
wiring, same drift, and any pause schedule affects both equally. Frozen
(plasticity off, schedule replayed) stays mandatory.

**Shared-input caveat, stated in advance because PM-4 taught it:** in the
contingent arm the pause schedule is computed from the very advances the
outcome counts, so the raw advance rate shares events with the
manipulation. The A-minus-B difference is what breaks that — B trials
are counted in the outcome and gate nothing. Only the difference gets
scored; the raw rate is descriptive.

Cost if it ever runs: 240 s per unit is 172 cycles, 86 A and 86 B, at
14.8 min — one legal unit, giving ~21 A trials per block over four
blocks. That is thin per wiring and the direction test across wirings is
what carries it; the power should be stated, not hidden. Two arms x
three wirings is 1.5 h, 45 min on two workers — a smaller ask than SM-2,
but still the founder's call and not mine.

---


<!-- private working record: commit b82a4ca, 2026-09-03 -->

## Run: DUR-1 (scoping) — the brief's lever is from the other preparation

Scoping only. Nothing run, nothing registered.

### 1. N8b's lever cannot be used here, and a null from it would mean nothing

The brief points at N8b's finding that `U` moved duration 6/6 (0.10
against 0.20 at matched drive). **That was measured on the cortical
dish**, where N8b established duration is *recruitment-limited*: the
burst ends when the recurrent loop stops recruiting, which is why a
weaker synapse ends it sooner.

**PM-4's anecdote is on the STG dish, and that preparation works the
other way round.** `## Run: STIM-1` established its rhythm is
intrinsically generated and that its synapses do not recruit at rest —
a stimulus that lands fires the cells it touches and *nobody else*, 2 in
2, 12 in 10.8. Burst duration there is set by the AB/PD cell's own
conductances, not by a recurrent loop.

So `U` — a synaptic parameter — is the wrong lever on this substrate.
It would very likely move nothing, and **a null produced by an inert
lever is indistinguishable from a null that refutes the hypothesis.**
That is the failure the whole `g_ie` episode was about: a parameter
producing a null that looks like a result.

### 2. A rate-matched lever does exist here, and it is already measured

From this line's own `g_ie` sweep, inside the measured bursting band:

    g_ie    burst rate    duration    guard
    0.03      0.73 Hz      437 ms     ok
    0.06      (default)      —        ok
    0.10      0.73 Hz      364 ms     ok
    0.25      0.13 Hz      238 ms     FAILS

**Duration moves 437 → 364 ms at an identical 0.73 Hz.** The mechanism is
coherent with STIM-1: inhibition truncates a burst whose length is set
intrinsically, so it shortens duration without touching the pacemaker's
rate.

That solves prospectively the covariate problem the brief anticipated for
the `U` lever — the two quantities that move together on the cortical
dish do **not** have to move together here.

### 3. But it is narrow, and it is confounded

- **Narrow.** 437 → 364 ms is a 17% range. PM-4's observed spread was
  388–451, about 16%. The lever's range is comparable to the variation
  whose meaning it is meant to test — enough for a direction, not enough
  for a dose-response curve worth the name. Beyond 0.10 the guard fails,
  and a run outside the guard is not evidence about bursts.
- **Confounded, and no unconfounded lever exists.** `g_ie` changes
  inhibition, which can reach STDP sorting by paths that have nothing to
  do with duration. Changing the cells' conductances instead would change
  excitability too. **There is no knob in this preparation that moves
  burst duration and nothing else**, and any registration claiming clean
  causation would be overclaiming.

### The design that is available

Two arms that fail differently, which is the strongest thing on offer:

- **Experimental, paired, confounded.** Each wiring run at `g_ie` ∈
  {0.03, 0.06, 0.10}. The comparison is *within* wiring, so the regime
  bifurcation PM-3 found cannot contaminate it and nothing is pooled.
- **Observational, unconfounded, correlational.** Across wirings at the
  fixed default `g_ie`, natural duration variation against sorting
  strength — PM-4's observation, now pre-registered as a prediction. It
  reuses the middle-level runs, so it costs nothing extra.

If a confounded intervention and an unconfounded observation agree on the
direction, that is real evidence. If they disagree, the confound is doing
the work and we learn that instead.

### The direction, stated because it is counterintuitive

PM-4 measured ρ(magnitude, duration) = **−0.47**: *shorter* bursts sorted
*harder*. So the prediction is that **more** inhibition — shorter bursts —
gives **stronger** sorting, which is backwards from "more time to
accumulate".

The mechanism that would explain it: sorting strength is carried by the
causal *fraction*, not the causal count (`## Run: PM-4` measured the
effect as a counting phenomenon with the STDP kernel's fine timing adding
nothing). A shorter, tighter burst orders its spikes more sharply, so a
larger share of pairs are causal even though there are fewer of them.
That is the claim the test would support, and naming it now stops it
being invented afterwards to fit whichever way the numbers fall.

### Price

At the measured 3.7 s wall per biological second, 60 s per run:

    3 levels x 5 wirings x 3.7 min = 55 min, in 4 batches under the line
    observational arm: reuses the g_ie = 0.06 runs, no extra cost

### Registration proposal, for the pass

- **P1** within each wiring, the causal-count excess is larger at
  `g_ie = 0.10` (364 ms) than at `0.03` (437 ms), in ≥ 4 of 5 wirings.
- **P2** the same holds for strong−weak `P/M`, in ≥ 4 of 5.
- **P3** across the five wirings at `g_ie = 0.06`, duration correlates
  negatively with sorting strength, ρ ≤ −0.4 — PM-4's anecdote, now
  pre-registered.

Covariate plan, named in advance:

- Burst rate is recorded per run. **If mean rate differs between levels
  by more than 10%, rate enters as a covariate and each level's runs are
  split at the median rate**, with the direction required to hold inside
  both strata. If it differs by more than 25%, the experimental arm is
  uninterpretable and only P3 is scored.
- Duration is recorded per run and reported per wiring. If a level does
  not actually move duration in a given wiring, that wiring contributes
  no evidence about duration and is reported as such rather than counted.

Falsifiers:

- P1 and P2 both fail → the anecdote dies as registered, and the
  duration story does not return without a new preparation.
- P1/P2 pass but P3 fails → the intervention moved sorting and natural
  duration does not track it, which points at the confound rather than at
  duration. Reported as evidence for `g_ie`, not for duration.
- The guard fails at any level in any wiring → that cell is void, not
  reinterpreted.

---


<!-- private working record: commit 59cf33b, 2026-09-03 -->

## Run: DUR-1 — the counterintuitive direction holds 5/5, and my scorer nearly buried it

### The scorer said the opposite, and it was wrong

The first thing this run produced was the line **"The anecdote dies as
registered"** on a 5/5 result.

The registration says a failed guard voids *that cell* — a (wiring,
level) pair. My scorer voided the whole **wiring**, requiring all three
levels to pass before a wiring could contribute. Two wirings lost their
`g_ie = 0.06` cell, so they were dropped from P1 and P2 — **comparisons
that contrast 0.03 against 0.10 and never touch the void level.** 5/5
was reported as 3/3, and 3 < 4 printed a refutation.

Nothing about the data was wrong. The instrument that reads the data was.
It is the same shape as the gate errors this line has been collecting all
week, and it is the most dangerous one yet, because **it failed in the
conservative direction** — a scorer that wrongly reports failure is one
nobody has any incentive to check.

Fixed to void per cell, as registered.

### The result

    seed   g_ie   dur ms  rate Hz  causal diff   P/M diff  guard
      30   0.03      451     0.73     +0.00051    +0.0100
      30   0.06        0     0.00          void      void   VOID
      30   0.10      290     0.60     +0.00957    +0.1215
      31   0.03      473     0.70     +0.00088    +0.0176
      31   0.06      449     0.73     +0.00160    +0.0318
      31   0.10      268     0.70     +0.00384    +0.0533
      32   0.03      465     0.73     +0.00192    +0.0443
      32   0.06      445     0.73     +0.00235    +0.0577
      32   0.10      253     0.73     +0.01322    +0.1909
      33   0.03      464     0.73     +0.00076    +0.0193
      33   0.06      444     0.73     +0.00147    +0.0376
      33   0.10      267     0.63     +0.01216    +0.1718
      34   0.03      450     0.73     +0.00068    +0.0147
      34   0.06        0     0.00          void      void   VOID
      34   0.10      280     0.67     +0.00911    +0.1228

    P1  causal excess larger at g_ie = 0.10   5/5   MET
    P2  P/M difference larger at 0.10         5/5   MET
    rate spread between levels 8.6% — under the registered 10%, no
    covariate needed

    duration        461 -> 272 ms   (-41%)
    causal excess   0.00095 -> 0.00958   10.1x
    P/M difference  0.02118 -> 0.13209    6.2x

**Shorter bursts sort harder, in every wiring, at matched rate.** The
direction was registered in advance and it is the backwards one: "more
time to accumulate" predicts the opposite and would already be refuted.
Three of the five wirings are monotone across all three levels.

### P3 is uninformative, not failed, and the difference matters

The observational arm needs natural duration variation at fixed `g_ie`.
There is none: the three surviving wirings at 0.06 span **444–449 ms,
1.0%**. There is nothing to correlate.

PM-4's ρ = −0.47 came from a 16% spread, and that spread existed only
because it included the two off-regime wirings (17 at 388 ms, 18 at 415).
Among on-regime wirings the natural duration variation is ~1%, which is
why the anecdote needed an intervention to test at all.

**So the registered falsifier does not fire.** It said *P1/P2 pass while
P3 fails → evidence for `g_ie`, not duration*. P3 did not fail; it did
not run. Applying a falsifier written for a different outcome would be
exactly the reinterpretation the falsifier existed to prevent.

### What this does and does not establish

**Does**: an intervention that shortens bursts by 41% strengthens sorting
by 6–10×, in 5 of 5 wirings, at matched burst rate, in the direction
predicted before the numbers existed.

**Does not**: separate duration from inhibition. `g_ie` was named a
confounded lever in scoping — *"there is no knob in this preparation that
moves burst duration and nothing else"* — and the arm meant to break the
confound produced no signal. **Duration versus inhibition is unresolved**,
and this run does not license the sentence "burst duration controls
sorting".

The next design has to move duration without moving inhibition, and
nothing in this preparation does that. That is the open problem, and it
is the same one scoping named.

### A property of the substrate, found on the way

At `g_ie = 0.06`, wirings 30 and 34 fire for ~2.8 s and then fall silent
**permanently** — 2,496 spikes in the first bin, nothing for the
remaining 57 s, with weights healthy (mean 0.0148, 2% at zero, no capped
cells). Both wirings live at 0.03 and at 0.10.

A network that dies at an intermediate inhibition and survives on either
side of it is not on a monotone path to silence: **the preparation has a
silent absorbing state that a transient can push it into.** It connects
to the monostability theme this line already has, from the other side.

---


<!-- private working record: commit 6a47184, 2026-09-03 -->

## Run: DUR-2 (scoping) — an intrinsic lever, and what triangulation can and cannot buy

Scoping only. Nothing registered, nothing scored.

### The lever, picked by measurement

Conductance order is `(Na, CaT, CaS, A, KCa, Kd, H, leak)`; the validated
AB/PD set is `[100, 2.5, 6.0, 50, 5.0, 100, 0.01, 0.0]`. Two candidates
terminate a plateau: raising `g_KCa` or lowering `g_CaS`. Measured on an
isolated cell rather than argued from the textbook:

    baseline (CaS 6.0)     495 ms   0.70 Hz
    g_KCa x2.0 (10.0)      156 ms   1.00 Hz     rate +43%
    g_CaS x0.6 (3.60)      209 ms   0.70 Hz     rate unchanged

**`g_CaS` is the lever.** `g_KCa` would drag rate up 43% and reintroduce
into the intrinsic arm exactly the covariate the triangulation exists to
avoid.

### Verified in the network before being proposed

Fourth preparation, same trap — the conductance is checked at the point
of use. Seed 31, `g_ie = 0.06`, 60 s:

    CaS x1.0   449 ms   0.73 Hz   1095 spikes/burst   guard ok
    CaS x0.8   327 ms   0.63 Hz    317               guard ok
    CaS x0.6   196 ms   0.57 Hz    191               guard ok

Guard, rate and regime all survive. **But in the network the lever does
move rate** — 0.73 → 0.57, a 25% spread, where on the isolated cell it
did not. The isolated-cell result would have been the wrong number to
plan on, which is why it was not.

### That rate behaviour is a feature, not a defect

DUR-1's `g_ie` lever held rate nearly fixed (8.6% spread) while cutting
duration 41%. This lever cuts duration 56% while dropping rate 25%.
**The two levers have different rate profiles**, so if sorting tracked
*rate* rather than duration, the two arms would disagree. Triangulating
with a lever that shares the confound would buy nothing; this one does
not share it.

### The limitation that must be registered up front

**Spikes per burst falls with duration under both levers**, and hard:
1095 → 191 here, and in DUR-1 wiring 30 the `g_ie` arm went 26,736 → 6,458
spikes while sorting rose 19×. The two quantities are collinear in *both*
arms.

So this design separates duration from **inhibition** (different lever)
and from **rate** (opposite rate behaviour). It does **not** separate
duration from **spikes per burst**, and no combination of these two
levers can. A pass would license "sorting tracks burst duration or the
spike count that falls with it" — and that sentence is the honest ceiling
of the experiment, stated before it runs rather than conceded after.

### Price: one arm, as asked

Run at `g_ie = 0.03`, where DUR-1 had all five wirings alive. **That makes
DUR-1's `g_ie = 0.03` runs the `CaS x1.0` baseline at zero marginal
cost.**

    new: CaS x0.8 and x0.6, 5 wirings, 60 s     = 37 min, batched
    baseline: already run                        = 0

### Registration proposal

- **P1** within each wiring, sorting is stronger at `CaS x0.6` than at
  `x1.0`, in ≥ 4 of 5 — the same direction DUR-1 found, through a lever
  that touches no synapse.
- **P2** pooling both arms against **realized** duration, sorting
  correlates ρ ≤ −0.5, **and the relationship holds within each arm
  separately** — lever identity as a registered covariate, so a pooled
  correlation created by two clusters at different offsets cannot pass.
- **P3** rate does not account for it: within the intrinsic arm, where
  rate falls *with* duration, sorting still rises as duration falls.

Falsifiers:

- **Intrinsic shortening does not strengthen sorting → "duration" dies**,
  and DUR-1's result belongs to inhibition. Stated in those words.
- P2 passes pooled but fails within an arm → the correlation is two
  clusters, not a relationship, and is reported as an artifact.
- Sorting tracks spikes per burst better than realized duration → nothing
  is separated, and the result is reported as unresolved rather than as
  support for either.
- Guard failure voids **that cell**, not that wiring — DUR-1's scorer bug,
  written into the registration so the next scorer cannot repeat it.

---


<!-- private working record: commit 5f22af2, 2026-09-03 -->

## Run: DUR-2 — the triangulation holds, and the registered ceiling holds with it

Intrinsic arm: `g_CaS` lowered inside the AB/PD cell, touching no
synapse, at `g_ie = 0.03` where DUR-1 had all five wirings alive. DUR-1's
runs are the `CaS x1.0` baseline, so this cost one arm.

    seed          arm   dur ms   rate  spk/brst     causal       P/M
      30     baseline      451   0.73      1215   +0.00051   +0.0100
      30     CaS x0.6      265   0.40       363   +0.00820   +0.1114
      31     baseline      473   0.70      1273   +0.00088   +0.0176
      31     CaS x0.8      428   0.73      1008   +0.00155   +0.0331
      31     CaS x0.6      279   0.47       385   +0.00825   +0.1730
      32     baseline      465   0.73      1259   +0.00192   +0.0443
      32     CaS x0.6      258   0.53       265   +0.00316   +0.0531
      33     baseline      464   0.73      1257   +0.00076   +0.0193
      33     both intrinsic levels VOID
      34     baseline      450   0.73      1215   +0.00068   +0.0147
      34     CaS x0.6      282   0.50       340   +0.00713   +0.0998

    P1  stronger at CaS x0.6 than x1.0        4/4 of the wirings with
                                              data, 4 of 5 overall   MET
    P2  rho(duration, sorting)   pooled -0.733
                                 intrinsic -0.648
                                 inhibitory -0.808                    MET
    P3  within the intrinsic arm, shortest beats longest              MET

**All three met.** Sorting tracks realized burst duration across two
levers with different confounds — one synaptic, one intracellular — and
the relationship holds *within* each arm as well as pooled, so it is not
two clusters at different offsets.

**The rate argument is now made rather than assumed.** In the inhibitory
arm rate held near-fixed (8.6% spread) while sorting rose; in the
intrinsic arm rate *falls* with duration (0.73 → 0.40–0.53) and sorting
rises anyway. Rate cannot drive both.

### The tie-check did not fire, and it should not be read as a separation

    rho(duration, sorting)        -0.733
    rho(spikes/burst, sorting)    -0.714

The registered clause fires when the spike correlation **matches or
beats** duration's. It does not, by **0.019** — on Spearman correlations
over twenty points, which is nothing.

**So the registered ceiling stands exactly as written before the run:**
the mediator is *"burst duration, or the spike count that falls with
it"*. Duration edging spike count by two hundredths of a rank
correlation is not evidence that it is the one, and reporting it as a
separation would be precisely the overclaim the ceiling was written to
forbid. The two quantities are collinear in both arms, as predicted, and
this design cannot separate them. Nothing here licenses dropping the
second half of that sentence.

### What discovery #2 gains

A mediator, at the strongest level this substrate permits: **shorter
bursts sort harder, through two independent interventions, monotonically
in realized duration, at 6–10× effect sizes, with rate excluded.** What
remains open is duration against spike count, and separating those needs
a lever that changes one without the other — which is a new preparation
problem, not a new run of this one.

### The void pattern is a finding, not attrition

    CaS x0.8   VOID in 4 of 5 wirings
    CaS x0.6   VOID in 1 of 5

**The intermediate level kills more than the extreme one**, and DUR-1
showed the same shape: `g_ie = 0.06` killed two wirings while 0.03 and
0.10 left all five alive. Two different parameters, both with a death
band in the middle rather than at the end.

That is not attrition to be tidied away — it says the silent absorbing
state of `## Run: DUR-1` is reachable from a *band* of parameter values,
not monotonically approached. `SILENT-1` now has a second instance and a
shape.

---


<!-- private working record: commit c2c9722, 2026-09-03 -->

## Run: BAND-1 (registration and price) — nothing run

The evidence, from DUR-1 and DUR-2, both levers, five wirings each:

    g_ie  0.030    0/5 dead
    g_ie  0.060    2/5 dead      <- interior
    g_ie  0.100    0/5 dead

    CaS   1.0      0/5 dead
    CaS   0.8      4/5 dead      <- interior
    CaS   0.6      1/5 dead

Two unrelated parameters — one synaptic, one intracellular — each kill
more in the middle of their range than at either end. Either the
absorbing state of `## Run: DUR-1` is reachable from a *band*, or 2/5 and
4/5 out of five wirings is small-n noise that happened twice.

### Registered prediction

**P1** On a finer grid with ten wirings per level, the death probability
is **non-monotone in the lever level, with an interior maximum**, for
**both** levers: some interior level kills strictly more than both the
lowest and the highest level tested.

### Falsifier, and it closes the row

**Death probability rises (or falls) monotonically across the finer
grid** → the band is small-n noise from five wirings, BAND-1 is wrong,
and the row closes. Reported in those words. A flat profile also closes
it: no interior peak, no band.

### The grid step, justified rather than assumed

The standing rule is that a grid coarser than the transition will step
over it — and **here the band *is* the transition**, so the step needs an
argument rather than a habit.

It has one: the grid is **anchored on a point already known to be inside
the band**. For `g_ie`, 0.060 kills and 0.030 and 0.100 do not, so the
band contains 0.060 and lies within (0.030, 0.100). Sampling
0.030 / 0.045 / 0.060 / 0.075 / 0.090 puts a known-interior point in the
grid by construction: the grid cannot miss the band entirely, whatever
the band's width. The same for `CaS` at 1.0 / 0.9 / 0.8 / 0.7 / 0.6, with
0.8 known interior.

What a step of 0.015 *cannot* do is resolve a band narrower than itself —
that would show as a single killing level with alive neighbours, which is
still non-monotone and still answers P1. **Resolving the band's width is
not the question; whether it is a band at all is.**

### The detector, and it gets validated before it is used

Death is fast: the one wiring examined in detail fired for 2,821 ms and
then never again in 57 s. So a **10-second run** should classify, at a
quarter of the cost of 60 s.

**That rests on a single observation, so it is validated first.** Stage 0
re-runs the five `g_ie = 0.060` wirings at 10 s and must reproduce the
60-second classification **5/5** (30 and 34 dead, 31/32/33 alive). If it
does not, the short detector is wrong and the map is redesigned rather
than run on it.

A late death misclassified as alive biases every level's rate *down*
uniformly, which cannot manufacture an interior peak — so the residual
risk after validation is masking a band, not inventing one.

    dead := the network spiked in the first half and produced no spike
            in the final 5 s. Distinguishes death from never having
            started, which is a different failure.

### Price, single-core, at the measured 3.7 s per simulated second

    stage 0  detector validation   5 runs x 10 s   =  3.1 min
    stage 1  g_ie map    5 levels x 10 wirings x 10 s  = 30.8 min
    stage 2  CaS map     5 levels x 10 wirings x 10 s  = 30.8 min

    total if both arms run                             = 64.7 min

**Stage 1 is a gate**: if `g_ie` comes back monotone, the falsifier fires
and stage 2 never runs — 34 minutes instead of 65.

**This exceeds the 15-minute line and is therefore not launched.** Stage 0
alone (3.1 min) is under it. Batched, stage 1 is three launches of about
10 minutes each. The decision on whether to spend it is the
orchestrator's, and this section is the price he asked for rather than
the run.

### BAND-1 amendment, recorded before any run

**Arm order swapped: `CaS` first.** It carries the stronger signal (4/5
interior against 1/5 extreme, versus `g_ie`'s 2/5 against 0/5), and the
original ordering was chronological habit, which is not a reason.

**The gate, rewritten for the new order and stated before the numbers:**

- **`CaS` monotone or flat on the finer grid → the falsifier fires and
  `g_ie` does not run.** The registered justification: if the *strong*
  evidence dissolves, the weak evidence cannot carry the row alone.
- **`CaS` confirms the band → `g_ie` runs anyway.** P1 asks for both
  levers, and a band on one lever is a different result from the one
  registered. A confirmation on `CaS` does not licence stopping early
  and claiming P1.

**Stage 0 is unchanged and its bar is hard.** The detector is
lever-agnostic, so validating it on `g_ie = 0.060` serves both arms.
**5/5 or the map is redesigned** — a 4/5 is a redesign, not something to
interpret.

### BAND-1 redesign, recorded before the map runs

The 10-second detector was wrong for this arm: `CaS` deaths are **late**.
Probe, seed 30 at CaS 0.8: fires normally for **22,135 ms** then stops,
against the `g_ie` deaths at 2,821 ms that stage 0 validated on.

**A claim of mine withdrawn.** I said ten wirings would be "one sample
printed ten times" because their per-second spike counts were
bit-identical. That over-read it. Identical counts are not identical
systems — the masks differ, the weights evolve differently, and the
aggregate rate merely coincides because intrinsic bursting dominates. The
probe proves it: seed 30 dies at 22 s and seed 31 runs to 60 s. Two
systems different all along, visibly different only when one fell over.
**The sampling was valid; only the window was wrong.**

Final design, estimated at 60.4 min single-core:

    levels   1.0 / 0.9 / 0.8 / 0.6   (known-interior 0.8, a shoulder each side)
    window   35 s, quiet tail 10 s   (sees deaths before 25 s)
    wirings  7 per level             (the falsifier names five as the problem)

**Sensitivity ceiling, in the required words:** deaths later than ~25 s
are missed; this biases every level's rate **uniformly downward**, which
**can mask a band but cannot invent one**. It is the safe direction and it
is a real limit on what a null would mean.

**Re-runs of externally killed batches do not count as new spend against
the cap** — a kill is not new science — **but every kill is reported.**

---


<!-- private working record: commit 72b1d0b, 2026-09-03 -->

## Run: BAND-1 CaS arm — the band is real, and the gate says `g_ie` still runs

Seven wirings, four levels, 35 s window with a 10 s quiet tail.

    seed     1.0     0.9     0.8     0.6
      30   alive   alive    DEAD   alive
      31   alive   alive   alive   alive
      32   alive   alive    DEAD   alive
      33   alive   alive    DEAD    DEAD
      34   alive   alive    DEAD   alive
      35   alive   alive    DEAD   alive
      36   alive   alive    DEAD    DEAD

    level   dead    rate
      1.0   0/7     0.00
      0.9   0/7     0.00
      0.8   6/7     0.86     <- interior maximum
      0.6   2/7     0.29

**P1 MET for this arm.** Death is non-monotone with an interior maximum:
0.86 at `CaS 0.8` against 0.00 at the top of the range and 0.29 at the
bottom. **The falsifier does not fire** — the profile is neither
monotone nor flat.

The upper edge is sharp: 0/7 at 0.9 and 6/7 at 0.8, with nothing in
between. The lower shoulder is soft: 2/7 still die at 0.6.

### Why the fresh wirings matter more than the others

Seeds 30–34 are the wirings the band was *noticed* in, during DUR-1 and
DUR-2. Finding it there again is consistency, not independent
confirmation. **Seeds 35 and 36 had never been run**, and both reproduce
the shape — 35 with the canonical pattern (dead only at 0.8), 36 dying at
0.8 and 0.6.

### The gate, applied as written rather than as convenient

The registration says: *`CaS` confirms the band → `g_ie` runs anyway,
because P1 asks for both levers and a band on one lever is a different
result from the registered one.*

So this is **not** BAND-1 answered. It is one arm of two, and the
registered claim is not available until the second arm runs. The
temptation at this point is to call an 0.86-against-0.00 result decisive
and stop; the gate was written before the numbers precisely to remove
that option.

### Cost, and the decision it forces

    stage 0 detector validation                3.0 min
    probe (death time + divergence)            7.2 min
    CaS map, 7 wirings                        59.0 min

    g_ie arm at the same window and n         60.4 min   (not yet run)

### Kills

**Six externally-killed runs this session**, all SIGTERM. No duration
threshold explains them: 8.6 min survived twice and died once, 12.3 died,
17.3 survived, 25.9 died two minutes in. The machine was under memory pressure at the time.

**The mitigation that worked was small batches**: one wiring is 8.6 min,
so a kill now costs at most that, and every completed wiring is on disk
before the next starts.

### BAND-1 `g_ie` arm — second exception, recorded before the numbers

**Second exception granted by the orchestrator**, 60.4 min, taking
BAND-1's total past 120 against an original 65-minute cap. Reason on the
record: a registered prediction answered half way is the worst kind of
waste — 59 minutes already spent license no claim at all until the second
lever speaks.

**Design symmetry is the point of the parameters**: same 35 s window,
same 10 s tail, same seven wirings, four levels arranged like the `CaS`
arm — a known-interior point, a shoulder either side, and one far
extreme.

    CaS   1.0 extreme / 0.9 shoulder / 0.8 interior / 0.6 far
    g_ie  0.030 extreme / 0.045 shoulder / 0.060 interior / 0.090 far

The two arms are then comparable without clauses.

**Seeds 35 and 36 are the independent replication for this arm too.** The
`g_ie` band was noticed on wirings 30–34 in DUR-1; those five cannot
confirm it, only agree with themselves. 35 and 36 have never been run at
any `g_ie` level.

---


<!-- private working record: commit 7c33c56, 2026-09-03 -->

## Run: BAND-1 — both levers, both non-monotone, both extremes clean

Seven wirings per arm, four levels each, 35 s window with a 10 s quiet
tail, detector validated 5/5 before use.

    CaS      1.0     0.9     0.8     0.6        g_ie    0.030   0.045   0.060   0.090
      30   alive   alive    DEAD   alive          30   alive    DEAD    DEAD   alive
      31   alive   alive   alive   alive          31   alive   alive   alive   alive
      32   alive   alive    DEAD   alive          32   alive   alive   alive   alive
      33   alive   alive    DEAD    DEAD          33   alive   alive   alive   alive
      34   alive   alive    DEAD   alive          34   alive    DEAD    DEAD   alive
      35   alive   alive    DEAD   alive          35   alive    DEAD    DEAD   alive
      36   alive   alive    DEAD    DEAD          36   alive    DEAD   alive   alive
    rate    0.00    0.00    0.86    0.29        rate    0.00    0.57    0.43    0.00

**P1 MET on both levers. The falsifier does not fire on either.**

The `g_ie` arm is the cleaner of the two: **zero deaths at both extremes**
and a majority in the middle. `CaS` has a soft lower shoulder — 2/7 still
die at 0.6 — while its upper edge is sharp, 0/7 at 0.9 against 6/7 at
0.8.

**The independent replication holds.** Seeds 35 and 36 had never been run
before this experiment; the band was noticed on 30–34, so only 35 and 36
can confirm rather than merely agree. Both show it, in **both** arms —
35 dies at both interior levels of each lever, 36 dies interior in each.

So the absorbing state of `## Run: DUR-1` is reachable from a **band of
parameter values in the interior of a range**, not approached
monotonically, and this is true of two unrelated parameters — one
synaptic, one intracellular.

### A wiring that never falls in

Seed 31 survives every level of both levers, and seeds 32 and 33 survive
every `g_ie` level. That is not part of P1 and is not claimed as a
result, but it is visible and worth recording: **the band is a property
of the parameter range, and susceptibility to it is a property of the
wiring.** Some networks simply do not enter it.

### Cost, and an overrun I should have flagged sooner

    stage 0 detector validation      3.0 min
    probe                            7.2 min
    CaS arm                         59.0 min   (estimated 65)
    g_ie arm                        85.5 min   (estimated 60.4)

The `g_ie` arm overran by 25 minutes, **entirely because one run took
34.9 minutes where every comparable run took 8.4**. The machine degraded
during the night: free memory fell from 87 MB to 65, nine runs were
killed externally, and the measured 3.7 s/s cost constant — reliable to
within 3% across dozens of runs all day — stopped holding.

**The process miss is mine**: I recorded the slowdown in the moment but
did not stop to tell the orchestrator that the authorisation was about to
be breached. A budget granted in minutes is a budget in *machine*
minutes, and when the machine changes underneath, the right move is to
report before spending the rest, not after.

**A measured cost constant is valid only while the hardware under it is
stable.** A cost measured before a change in the machine underneath is a
plan against a machine that no longer exists.

### BAND-1: the 35-second detector, validated after the fact

The 10-second detector was validated 5/5 and then failed on the other
lever, which is why the window moved to 35 s. **The 35-second window was
never validated in turn** — it was set from a single death-time
observation (22,135 ms) and used.

It can be validated, from data already in hand. DUR-1 and DUR-2 ran
fifteen of the same (wiring, CaS) cells at 60 s under a different
criterion — the burst guard rather than a quiet tail:

    35 s detector vs 60 s ground truth, every cell both experiments ran
    15 / 15 agree

Two different windows, two different criteria, the same verdict on every
cell. That is stronger than the original stage 0, which was five cells at
one level of one lever.

**Stated here because the section did not claim it and should have.** The
evidence existed before the map ran and nobody looked; a referee would
have been right to attack the unvalidated window, and the answer would
have been sitting in a file the whole time.

### BAND-1: the uniform-bias assumption, tested rather than asserted

The registration conceded a sensitivity ceiling — deaths later than ~25 s
are missed — and argued the bias is **uniform across levels**, so it can
mask a band but not manufacture one. **That uniformity was asserted.**

It cannot be tested the way it first appeared: **the saved traces do not
contain it.** BAND-1 stored only spike *counts* per window, and DUR-1/2
stored spike times only from t = 30,000 ms, so a death at 22 s leaves an
empty trace and no recoverable death time. The "zero compute" route was
not available.

A different free test was. Fifteen cells were run at **both** 35 s and
60 s, and a late death — one occurring between 35 s and 60 s — would
appear as *alive at 35, dead at 60*:

    CaS 1.0  (extreme)   5 cells   0 changed verdict
    CaS 0.8  (interior)  5 cells   0 changed verdict
    CaS 0.6  (extreme)   5 cells   0 changed verdict

**No cell changed verdict, at either extreme or in the interior.** So in
the 35–60 s window there are no level-dependent late deaths to bias
anything, and the concern is answered for that interval with n = 5 per
level.

**What remains open, precisely**: nothing has ever been run past 60 s in
this preparation. A death at 90 s would be invisible to every measurement
on this board. The claim is now *"no late deaths between 35 and 60
seconds"*, which is a measurement, rather than *"the bias is uniform"*,
which was a hope.

---


<!-- private working record: commit 6d37942, 2026-09-03 -->

## Run: MEA-1 (registration, phase 1) — no download, no analysis

Testing §5's observational half on archived public recordings:
Wagenaar, Pine & Potter (2006), dissociated cortical cultures on MEAs,
spontaneous activity across development. **Nothing has been downloaded
and nothing has been looked at.**

### The circularity that would have invalidated everything

§5 predicts that *stronger* connections carry an excess of causally
ordered pairs. On an MEA there are no synaptic weights: connection
strength must be *estimated from spikes*. **If strength is estimated from
causally ordered pair counts — which is what a cross-correlogram peak
essentially is — then "strong connections have more causal pairs" is
true by construction and the analysis measures nothing.**

The model had no such problem: weights were known independently of the
spikes used to score them. The MEA version does, and it is fatal if
unnoticed.

**Registered solution: split-half by burst.** Burst windows in each
recording are split into two interleaved halves (odd- and even-numbered
bursts). Connection strength is estimated on half A; the causal-pair
excess is measured on half B; then swapped, and the two directions
averaged. **The spikes that define a connection's class are never the
spikes that score it.** Any analysis that scores strength and excess on
the same events is void, not adjusted.

### Instruments, with the inherited number refused

- **Strength proxy** (on the estimation half): the peak of the smoothed
  cross-correlogram at |lag| ≤ 25 ms, normalized by the product of the
  two units' firing rates, computed over burst windows only. Ordered
  pairs are ranked and split at the **median**, exactly as the model
  splits synapses. Fixed before any data is seen.
- **Causal-pair statistic** (on the scoring half): the fraction of
  pre-before-post pairs among all pairs inside burst windows, a pure
  count with no kernel weighting — the quantity `## Run: PM-4`
  established as the carrier.
- **The lag window is NOT inherited as 250 ms.** That number was chosen
  for ~450 ms model bursts; cortical culture bursts differ, and this line
  has been burned four times by carrying a number across preparations.
  The primary window is **half the median burst duration of that
  recording**, computed per recording. The fixed 250 ms is reported
  alongside for comparability with the model, descriptively, never
  scored.
- **Burst detection** is the published method of the dataset's own
  literature rather than this board's 10 ms/20% detector, for the same
  reason. Which method is fixed at phase 2 from the dataset
  documentation, before data is opened.

### Registered predictions — directional only

  **P1** Within a recording, the strong half of connections shows a
  higher causal-pair fraction than the weak half, in **≥ N−1 of N**
  recordings analysed.

  **P2** Across recordings, sorting strength (strong-minus-weak causal
  excess) correlates **negatively** with median burst duration —
  **within an age stratum**, not pooled (see confounders).

  **P3** First-spike-latency order carries no signal: Spearman ρ between
  a unit's median within-burst first-spike latency and its estimated
  incoming strength is not consistently negative — the refuted temporal
  mechanism should be absent in tissue too.

No magnitude is predicted. §5 already states the model-derived
tenths-of-a-percentage-point scale does not export, and this registration
does not smuggle it back in.

### Falsifiers, in closed words

- **P1 fails** → the counting mechanism does not appear in cortical
  tissue under an independent strength proxy, and §5's first prediction
  is refuted on observational data. Reported as a refutation of our own
  headline, not as a proxy problem.
- **P2 fails within strata but holds pooled** → Simpson artifact, and the
  duration relationship is a developmental confound. Reported as such.
- **P3 shows a consistent negative correlation** → the temporal mechanism
  we refuted in simulation is present in tissue, and the model is wrong
  about the mechanism even if right about the outcome.
- **Any recording whose strength estimate and causal score share spikes**
  → void, not adjusted.

### What this design cannot distinguish, registered now

Observational data on a developing culture confounds almost everything
with **age**. Burst duration, connection density, firing rate and
synaptic maturation all change together across DIV.

- **Age (DIV) is a registered stratifying variable**, not a covariate to
  regress away: P2 is scored **within** age strata and the pooled value
  is reported descriptively only. The Simpson rule applies prospectively.
- **Active-unit count and mean firing rate** are recorded per recording
  and reported beside every result. If sorting strength tracks unit count
  at least as well as burst duration, **nothing is separated** and the
  result is reported as unresolved — the same tie-clause that governed
  `## Run: DUR-2`.
- **No intervention is possible here.** §5's predictions 2 and 3 are
  interventional and **cannot be tested on this dataset at all**. MEA-1
  addresses the observational half only, and a pass does not license the
  interventional claims.

### Download perimeter, fixed before the download

    cultures                 6, chosen by recording quality criteria
                             stated in the dataset documentation, not by
                             looking at their activity
    ages                     3 DIV points per culture, spanning the
                             dataset's developmental range
    recordings               18 total
    size                     UNKNOWN — to be measured on ONE file first
                             and reported before the rest is fetched
    analysis cost            UNKNOWN — timing probe on one recording
                             before the batch, priced to the orchestrator
                             if the total exceeds the 15-minute line

Both unknowns are registered as unknowns rather than guessed. This board
has twice been burned by a cost constant that stopped holding, and once
by a "zero-compute" re-analysis whose data did not exist.

### MEA-1 amendment: the per-pair null, and why it is not jitter

The requirement is right and the analysis needs it: class assignment
correlates with firing rate, and a pair's baseline causal fraction need
not be 0.5. Scoring against a nominal 0.5 would let rate back in through
the window after it was shown out of the door.

**But jitter cannot serve as that null here, and this can be shown
before any data is opened.** The statistic is the causal fraction among
pairs with |lag| ≤ W, where W is half the recording's median burst
duration — a *coarse* quantity, at the scale of the burst itself.

- **Jitter wide enough to matter (uniform within the burst window)**
  makes both spike trains uniform over that window. Two uniform trains
  give an expected causal fraction of exactly **0.5, for every pair,
  regardless of rate** — the null collapses into the nominal 0.5 it was
  introduced to replace, and buys nothing.
- **Jitter small compared to W** leaves the lag distribution essentially
  unchanged at scale W, so the null lands on top of the data and every
  pair's excess goes to ~0 — the opposite failure.

Jitter fails at both ends because it is a *fine-timing* null applied to a
*coarse-ordering* statistic. It is the right tool for correlogram peaks
at millisecond lags; it is the wrong one here.

**Registered instead: the burst-shuffle null.** For each pair, unit A's
spikes from burst *i* are paired against unit B's spikes from burst *j*,
i ≠ j, over N = 100 random derangements of the burst index. This

- **preserves** each unit's spike count, firing rate, and its typical
  position and shape within a burst — every property that could create a
  baseline away from 0.5;
- **destroys** only same-burst coordination, which is exactly the
  quantity a synapse is supposed to produce.

The reported excess is `observed − mean(null)` per pair, and each pair's
null spread is retained so that a pair whose excess is inside its own
null spread contributes as such rather than as a small positive number.

**What the shuffle deliberately keeps** is a static lead–lag between two
units — if A tends to fire before B in every burst, the shuffle keeps
that. This is intentional: `## Run: stg-P/M` refuted first-spike order as
the carrier, so a static lead is a confound here, not the signal, and a
null that removed it would credit us for the very thing we showed does
not carry the effect.

If the orchestrator prefers, a small-δ jitter can be reported alongside
as a descriptive check, but it is not the scored null and, for the reason
above, cannot be.

**Amendment addendum (orchestrator, minuted).** Keeping the static
lead–lag inside the null is *conservative* with respect to our own claim:
if a real synapse produced a static phase shift of the postsynaptic unit
on every burst, that contribution is removed along with the confound. In
the registered wording: **the null subtracts any synapse-induced static
phase shift; the measured excess is therefore a lower bound on same-burst
coordination.** The null can miss true signal but cannot manufacture it,
which is the correct direction of error for a test we intend to be able
to lose. Small-δ jitter may be reported as a descriptive check and is
never scored.

---


<!-- private working record: commit ccc5c1a, 2026-09-03 -->

## Run: MEA-1 phase 2 — one file measured, five corrections to the design

Downloaded exactly one recording (`1-1-25`, dense, spontaneous), measured
its size and the analysis cost, and read the dataset's own paper for the
burst detector the registration deferred to it. **The scientific outcome
of this recording was written to JSON and not displayed** — pricing a run
must not let a first glimpse of the result bend a frozen design.

### The numbers asked for

    one recording (1-1-25)   1.14 MB bz2 -> 4.43 MB text, 330,675 spikes,
                             56 active electrodes, 2712 s
    download                 1.8 s
    analysis, single core    30.0 s  (load 0.0, detect 0.1, slice 0.0,
                             pair scoring 29.8)
    full perimeter           18 files, 20.5 MB, all present
    projected analysis       ~10-15 min single core, one file (6-1-24,
                             6.8 MB) carrying ~3 min of it

Inside the 15-minute line without a second worker.

### The burst detector, taken from the dataset's literature as registered

SIMMUX, quoted from the paper's Methods: a *burstlet* is ≥4 spikes on one
electrode with every ISI below a threshold of ¼ of that electrode's
inverse mean spike rate, or 100 ms if that rate is below 10 Hz; any group
of burstlets overlapping in time across electrodes is a burst. Implemented
verbatim; our board's 10 ms/20% detector was not used.

### Correction 1 — duration is measured on the ASDR, not on the burstlet extent

Taking the burstlet extent as the burst duration gave a **1713 ms** median
at 25 div against the paper's published **<200 ms after 20 div** — a 9×
disagreement. The paper measures duration on the ASDR profile smoothed
with a 10 ms Gaussian, between the 20% points either side of the peak.
Implemented that way the median is **186 ms**, against the paper's figure.
**This is the first external validation on this board**: an implementation
agreeing with a published measurement made by other people on the same
data. Every previous check has been internal.

### Correction 2 — the single-linkage merge was gluing successive bursts

"Any group of burstlets that overlapped in time" is single-linkage, and a
long-tailed burstlet bridges two bursts: **24 of 84** detected bursts had
more than one ASDR peak. Burst windows are now the 20%-to-20% span around
*each* peak, which is the paper's own machinery applied consistently.

### Correction 3 — the lag window is 92 ms, not 250 ms

W = half the median burst duration = **92 ms**. The registration's refusal
to inherit the model's 250 ms was worth its ink: the correct value is 2.7×
smaller, and inheriting would have counted pairs across a window nearly
three times the structure it was meant to sit inside.

### Correction 4 — the first five minutes of every recording are discarded

The paper reports that moving a dish into the rig makes cultures fire
volleys of bursts, and that "mechanical perturbation increased the
synchronicity between neurons almost without increasing total firing
rates" — a transient that raises synchrony specifically, which is the very
quantity we measure, with the effect "mostly limited to the first 5
minutes". The first 300 s of every recording is dropped before anything
else. Registered now, before any pair was scored.

### Correction 5 — the registered DIV points would have given N=1 dressed as N=6

The perimeter registered 6 cultures at 3 DIV points. Only **five** cultures
have recordings at {12, 20, 28}, and **all five are from plating batch 2**.
The paper's own conclusion forbids exactly this: cultures from one batch
"developed along strikingly parallel lines", and it "emphasizes the
importance of using multiple preparations — not just multiple cultures
from one preparation". Six sister cultures would have been pseudo-
replication, and P1's "≥ N−1 of N" would have counted one preparation six
times.

**Amended perimeter, chosen structurally and before any activity was
looked at:** DIV **{10, 17, 24}**, the triple with ≥5 div spacing covering
the most plating batches (6 batches, 21 eligible cultures), and one
culture per batch — the lowest-numbered in each: **1-1, 2-1, 3-1, 4-1,
6-1, 8-1**. Six cultures, six preparations, 18 recordings. File sizes were
measured for cost only and played no part in this selection; the cost
probe `1-1-25` falls outside the scored perimeter.

### Two limits of this dataset, registered before results

- **The data is multi-unit and cannot be sorted.** The paper is explicit:
  "sorting was not attempted, and all results in this paper are based on
  multiunit data", because during bursts overlapping waveforms make
  sorting problematic. A "unit" here is an **electrode** carrying many
  cells, 200 μm from its neighbours. A strong cross-correlogram peak can
  therefore reflect a shared population rather than a synapse. The
  split-half design still blocks the circularity it was built for, but
  **P1 tests ordering between electrode populations, not between neurons**,
  and will be reported in those words.
- **Burst duration is largely determined by DIV** — the paper's Figure 6B
  has it falling from ~1 s at onset to <200 ms after 20 div. P2 is scored
  within age strata, so it lives on whatever duration spread survives
  inside a stratum. **Registered in advance:** if the within-stratum
  interquartile spread of median burst duration is under 25% of the
  across-stratum spread, P2 is reported as **underpowered on this
  dataset**, not as a refutation, and the pooled value stays descriptive.

---


<!-- private working record: commit d07e529, 2026-09-03 -->

## Run: MEA-1 — the prediction meets real tissue, and does not survive

18 recordings, 6 preparations (1-1, 2-1, 3-1, 4-1, 6-1, 8-1), DIV
{10, 17, 24}, Wagenaar–Pine–Potter dissociated cortical cultures.
Split-half by burst, burst-shuffle null, W computed per recording, first
300 s discarded. 2 recordings void on the registered guard (fewer than 8
bursts): `2-1-10`, `2-1-17`. **N = 16.**

### Verdicts, as registered

    P1  strong > weak in  6/16     (pass required >= 15)   FAIL
    P2  rho(sorting, duration) within strata  +0.08        FAIL
    P3  latency rho +0.35, negative in 1/14                PASS

**P1 is refuted.** Per the registered falsifier, in its own words: the
counting mechanism does not appear in cortical tissue under an independent
strength proxy, and §5's first prediction is refuted on observational
data. This is a refutation of our own headline, and it is not a proxy
problem.

**P2 is refuted and cannot hide behind either escape hatch.** The
registered underpower clause did not fire — within-stratum duration spread
was 2.32× the across-stratum spread, nine times the 0.25 floor — and the
DUR-2 tie-clause did not fire either, since sorting tracked unit count
(ρ = +0.05) *less* well than duration (ρ = +0.08). Both were available and
neither applied. Sorting does not track burst duration negatively in
tissue; the per-stratum values are +0.10, +0.40, −0.26.

**P3 passes.** First-spike latency order carries no consistently negative
signal, as the model predicted after `## Run: stg-P/M` refuted it there.
The relation is in fact consistently *positive* (ρ = +0.35, 13 of 14),
which this design cannot attribute — P3 controlled nothing for firing
rate — and which is recorded as descriptive only.

### An instrument fault found after the first result, and that result withdrawn

The first pass scored P1 at 7/17 and would have been reported as the same
refutation. It was wrong. `durations()` was validated against the paper's
published <200 ms at 25 div; the run then called `windows()`, a **second,
unvalidated function**, which fragmented one SIMMUX burst into ten and
collapsed W from ~90 ms to ~13 ms. The two agreed on the probe recording
— 186 vs 185 ms — which is exactly why it went unseen. **The external
validation announced in phase 2 covered code that was not the code that
ran.** Corrected to one window per SIMMUX region spanning the 20% points
around its dominant peak, reproducing the validated function on every
control file, and the whole perimeter re-run. Both verdicts are from the
corrected run; the first is void, not reported.

It survived only because a broken instrument printed 2575 bursts of 26 ms
where the paper documents about one second.

### What the data actually constrains — and the falsifier we mis-wrote

The sorting statistic across 16 recordings: mean **−0.105 pp**, sd 0.80 pp,
95% CI **[−0.53, +0.32] pp**, t = −0.52. So the data excludes any tissue
sorting effect above **0.32 percentage points**.

**The model's carrier was ~0.2 pp, which lies inside that interval.** And
P1's own pass rule — 15 of 16 recordings positive — required a
per-recording effect of **1.0 pp against a between-recording sd of 0.80 pp,
five times the magnitude the model itself produced.**

So P1 was **unpassable at the effect size we had reason to expect, and we
registered it without checking that.** This does not soften the verdict:
P1 is refuted as registered, and no goalpost moves. But it does bound what
the refutation means. The honest statement is a **bounded null**: sorting
in cortical tissue is below 0.32 pp, and an effect at the model's own scale
cannot be resolved by this dataset. "The mechanism is absent" and "the
mechanism is present at 0.2 pp" are not separated here, and the design
could never have separated them.

P2 carries no such caveat. It had power by its own registered test and
failed on the merits.

### The interventional half remains untested

§5's predictions 2 and 3 are interventional and were declared untestable on
this dataset before the download. They stay untested. A P1 failure does not
touch them, and neither would a pass have licensed them.

---


<!-- private working record: commit 6a9c893, 2026-09-03 -->

### DUR-2's published tie-check numbers are computed on a set with duplicates

The referee asked for a bootstrap CI on the 0.019 gap between ρ(duration)
and ρ(spikes/burst). Re-running the **registered scorer** on the persisted
JSONs reproduces the published numbers exactly — ρ(duration) = −0.733,
ρ(spikes/burst) = −0.714, gap 0.019 — so this is not a reproducibility
failure. It is a defect inside the number itself.

`dur2.score()` builds `pool = intr + inh`, where `intr` includes the
CaS×1.0 baseline rows and `inh` is every guard-passing DUR-1 row. **The
baseline rows ARE DUR-1's g_ie = 0.03 rows** — the shared baseline that
made the second arm cost one arm instead of two. They therefore enter the
pooled correlation **twice**: 5 of the 23 pooled points are exact
duplicates of 5 others.

    scorer pool, as run   n = 23   rho_dur −0.7332  rho_spk −0.7144  gap 0.0188
    de-duplicated         n = 18   rho_dur −0.7895  rho_spk −0.7564  gap 0.0330

**The verdict does not change.** The tie-clause fires when
|ρ_spikes| ≥ |ρ_duration|, and it stays silent in both — duration wins
narrowly either way, and P2 is MET either way. What changes are three
published figures and the n behind them.

The inconsistency is already visible inside our own outputs: `fig4()`
draws the shared baseline **once**, in its own marker, with a comment
saying it does so "so a reader who counts points counts correctly". The
figure de-duplicates and the scorer does not, and the paper reports both.

Not corrected in the preprint from here — the numbers are the fork's and
the founder's. Recorded with the de-duplicated values so whoever fixes it
has them.

**Ground rule this earns:** *a row that legitimately belongs to two arms
still belongs to a pooled set once.* The shared baseline was a good design
— it halved the cost of DUR-2 — and the double-count is the price of not
asking what "pooled" means when two arms overlap by construction.

---


<!-- private working record: commit e088ac5, 2026-09-03 -->

### DUR-2's tie-check, bootstrapped: indistinguishable, and the gate was a coin-flip away

On the de-duplicated pool (n = 18), 10,000 cell resamples:

    rho(duration, sorting)              −0.7895
    rho(spikes/burst, sorting)          −0.7564
    rho(duration, spikes/burst)         +0.9319
    |rho_spk| − |rho_dur|               −0.0330   95% CI [−0.1858, +0.0888]

**Indistinguishable — the interval spans zero.** The two predictors correlate
at +0.93 with each other; there was never enough independent variation in 18
cells to separate them, and no amount of care in the scoring could have
manufactured it.

The number that matters more than the interval: **in 27.4% of resamples the
registered tie-clause would have fired**, declaring the result UNRESOLVED
instead of crediting duration. The gate that let DUR-2 through was better
than a coin flip but not by much, and it passed on a margin of 0.033 with a
half-width of 0.14.

This does not retract DUR-2. The registered clause was applied honestly to
the data as it came, and it did not fire. But "duration wins over spike
count" was never established by this experiment — it was *not refuted* by
it, which is a weaker claim, and the honest wording is that the two
candidates cannot be told apart here.

**Ground rule:** *a gate that passes narrowly has told you the experiment
lacked the power to run it.* Report the margin next to the verdict, or the
verdict reads as a finding when it is a near-miss.

---


<!-- private working record: commit 44eaa45, 2026-09-03 -->

### The 0.8299 provenance sentence in Methods is false

Found while checking the bibliography, which is not where it was expected.

§2 states: *"The depression/potentiation area ratio A₋/A₊ was set to the
measured neutral point of the network's own lag distribution (0.8299) — the
value at which the rule neither systematically grows nor shrinks weights on
this preparation's activity — determined in a separate registered
experiment."*

Three things are wrong with that sentence.

1. **0.8299 was never measured.** It follows in closed form from the fixed
   integral: `A₋ = (A₊·τ₊ − integral)/τ₋` with integral = −0.134 gives
   A₋/A₊ = 0.829871 exactly. No data enters `rule_with_integral`. This is
   the same fact the (e) check established, and the correct sentence
   already sits **two lines below it in the same paragraph** — the
   paragraph now gives two incompatible provenances for one number.
2. **The separate registered experiment measured something else.** N7
   (`neutral_point.py`, `neutral_point.json`) measured the neutral ratio at
   **0.605 ± 0.046** (predicted 0.598, r = +0.738 across seeds) — and on the
   **cortical** dish, not this preparation.
3. **What N7 actually established was a principle, not this value**: that
   the dividing line is the rule's own A₋/A₊ rather than zero area. Every
   other place in this repository calls 0.8299 "the rule's own ratio" or
   "N7's correction", and both are right. Only the preprint calls it a
   measurement.

The coincidence that made this hard to see is recorded in our own text: STG
P/M values lie between 0.776 and 0.888 with 0.8299 in the middle of that
range. 0.8299 *does* sit near this preparation's operating point — which is
an observation made afterwards, not the way the value was chosen.

**Not corrected here.** Recorded for the fork with the true provenance:
the integral was fixed at −0.134 and the ratio follows from it; N7 justifies
using the rule's own ratio as the dividing line rather than zero.

**Ground rule:** *a number's provenance is a claim, and decays like any
other.* "Measured in a separate experiment" reads as a strength and is
rarely re-checked, which is exactly why it needs the same second route as
the number itself.

---


<!-- private working record: commit ecfedcf, 2026-08-29 -->

### Ground rules learned the hard way

These cost real time. Do not rediscover them.

- **Measure on the right preparation.** Measuring an EPSP on an intact cell
  contaminates it with the cell's own spikes. Block sodium (`block_sodium`) for
  subthreshold quantities; use `voltage_clamp` for conductances. Six separate
  measurements went wrong this way before the pattern was noticed.
- **Never change two things and attribute the difference to one.** A 2.9×
  speedup was credited to adaptive stepping; isolating the two changes showed
  adaptive stepping contributed 1.00× and lazy synapses contributed all of it.
- **Choose tolerances from a convergence curve, never by intuition**, and always
  check the firing rate before comparing two conditions. A benchmark run in the
  wrong regime produces a confident, wrong number.
- **Verify surrogate controls actually control.** A uniform spike-shuffle passed
  completely unconnected neurons at z = 5.3. Use `surrogate='shift'`.
- **A detached parameter does not error, it just stops learning.** Check that
  every trainable quantity is used as a `Var` inside the computation.
- **State objectives as relative properties, not absolute numbers.** Absolute
  targets the system cannot reach make the optimizer find degenerate solutions.
  Broken again in D2's first version, and caught by its own registered control.
- **Probe a trained network on a fresh copy, never on the one that just
  trained.** The residual activity carries the result, and what it produces
  looks exactly like success.
- **A scan that stops before the effect does invents a finding.** Four active
  dendrite configurations were recorded as having "no passive equivalent"; the
  passive scan had simply ended one step short of them. The lag cutoff in N7
  nearly did the same: 60 ms read 0.66, converged value 0.586.
- **Let the cell settle before measuring it.** Extra channels move the resting
  potential.
- **Build the experiment so that a failure produces an impossible number, not a
  plausible one.** A broken run that reads 12.7 gets written up; one that reads
  exactly 1.0 gets caught.
- **A measure taken over a whole run is an average in disguise.** Anchor the
  measurement window to the stimulus, not to the simulation.
- **When a system's own time constant changes between conditions, a fixed
  window compares an equilibrium with a transient.** Measure at the half-way
  point and at the end.
- **Any number that reaches an abstract, a title, or a figure gets two
  independent routes to it before the PDF.** Not one route checked twice --
  two: re-run the instrument, AND recompute the quantity from the persisted
  data. Their agreement is worth little; their DISAGREEMENT is the whole
  point. DUR-2's published gap of 0.019 reproduced perfectly when the
  registered scorer was re-run, and only the second route showed that its
  pooled set counted five of twenty-three rows twice. Nobody was looking for
  that error, and no amount of re-running the same code would ever have
  found it. Promoted to lane practice by the orchestrator on 2026-09-03.
- **A shared row belongs to a pooled statistic once**, however many arms it
  legitimately belongs to. DUR-2's baseline was a good design that halved
  the experiment's cost; the double-count was the price of not asking what
  "pooled" means when two arms overlap by construction.
- **A gate that passes narrowly has told you the experiment lacked the power
  to run it.** DUR-2's tie-clause did not fire, but it would have fired in
  27.4% of bootstrap resamples. Report the margin beside the verdict, or the
  verdict reads as a finding when it is a near-miss.
- **Validate the code that runs, not a sibling of it.** MEA-1's duration
  measurement was checked against a figure other people published on the
  same data and agreed to 1 ms — then the run called a second, unvalidated
  function for the same job, which fragmented one burst into ten. The two
  agreed on the one recording used to validate and disagreed by 10×
  everywhere else. A validation is only worth the code path it touches.
- **Register a falsifier you could pass at the effect size you believe in.**
  MEA-1's P1 required 15 of 16 recordings positive, which needed a 1.0 pp
  effect against 0.80 pp of between-recording spread — five times the
  ~0.2 pp the model itself produced. It was unpassable before the first
  file was downloaded. The verdict still stands as registered, but a
  falsifier that cannot fire distinguishes nothing: work out its minimum
  detectable effect at registration time, and say so there.
- **Replication counts preparations, not siblings.** MEA-1's registered
  perimeter asked for six cultures at three ages; exactly five cultures had
  all three, and all five came from one plating batch. Six sister cultures
  would have entered the arithmetic as N=6 and been N=1, and the pass rule
  "≥ N−1 of N" would have counted a single preparation six times. The
  dataset's own authors had already written the warning — cultures from one
  batch "developed along strikingly parallel lines" — and it was found by
  checking which files exist, not by reading results. Before any N is used
  as an N, ask what the independent unit actually is.
- **Register the covariate as a split before looking, or Simpson decides
  for you** (pure example: the mean over seeds in N8's calibration). And
  its stronger companion, from phase-h's correction: **several registered
  checks failing together is evidence about the instrument, not about the
  world.** A sign, a pooling structure and a confound's direction all broke
  at once under the undiscounted predictor and all healed under the
  discounted one — a real effect rarely breaks three registered checks
  simultaneously.
- **Registered falsifiers catch wrong theories; invariants catch wrong
  arithmetic.** A minimum cannot exceed its population, a probability
  cannot exceed one, a subset cannot outnumber its set — T1d's first
  number (a shared step of 0.2525 against a per-neuron 0.1333 it
  minimises over) announced its own wrongness at no cost (step-mean vs
  time-mean; the harmonic time-weighted mean is the step-count quantity).
  In this document wrong arithmetic has been the more common failure than
  wrong theory. Companion warning: a measurement can answer the wrong
  question convincingly — T1d's number is clean, monotone, reproducible,
  and comes from a network whose synchronous firing excludes the effect
  asked about; the tell was in the spike count, not the number.
- **Verify the manipulation landed before interpreting the response.** A
  sham-controlled calibration with a genuinely zero floor still reported
  0.27: four of five "responses" sat on trials where the pulse never fired
  the stimulated cells. The sham checks the world without the manipulation;
  this checks the manipulation happened — they are different controls and
  this line now carries both.
- **A variable you discover is not a variable you can control.** PM-3's
  regime looked like a blocking factor because it was bimodal and
  separable; it was downstream of the process under study, and one wiring
  converted regimes mid-run. The cheap test that distinguishes a factor
  from an outcome: can you assign it before the outcome exists?
- **A criterion with two clauses will eventually split — report the split
  as a split.** PM-2's "positive in ≥ 8/10 AND exceeding the spread" was
  written precisely so a consistent-but-variable effect could not be called
  a win; it then produced exactly that, and the verdict led with the
  failure. Its companion: the double-barrelled registration earned its keep
  in the designed direction — when the mechanism barrel failed while the
  outcome held, a single barrel would have concluded "no order, so no
  sorting" and been wrong.
- **A gate that fails conservatively is the one nobody checks.** DUR-1's
  scorer voided wirings instead of cells and printed "the anecdote dies"
  on a 5/5 result — caught only because the conclusion contradicted a
  table read row by row. A wrong PASS invites scrutiny; a wrong FAIL
  invites a shrug. Validate the scorer against a hand-read table before
  trusting its verdict in either direction.
- **Verify the artifact's behaviour, not the inputs you handed it.** A
  verification that reads the values fed to a routine, rather than what
  the routine emits, will print "fixed" with the defect still in place —
  correct inputs read by wrong code. The acceptance is running the thing
  and checking a known value. (Encountered on a product-side defect in an
  export path, reported separately and outside this release.)
- **A shared baseline enters a pooled statistic once** (DUR-2 correction,
  2026-09-03). Sharing DUR-1's g_ie = 0.03 runs as DUR-2's CaS ×1.0
  baseline halved the cost of the second arm — and the scorer, built as
  `pool = intr + inh`, counted those five cells twice (n = 23 instead of
  18). The figure de-duplicated; the scorer did not; the paper reported
  both for half a day. Verdict unchanged, three published numbers not.
  Whenever a design reuses runs across arms, the pooled estimator must
  de-duplicate by run identity, and the n must be printed next to every
  pooled number so a reader can check it against the design.
- **A validated instrument is validated for the regime it was validated
  in** (BAND-1, 2026-09-03). A 10 s death detector, validated 5/5 on the
  `g_ie` arm where deaths are early (~2.8 s), was applied to the CaS arm
  where deaths are late, plasticity-driven — and read a cell as alive
  that DUR-2 had recorded as VOID at 60 s. The orchestrator wrote "the
  detector is lever-agnostic", the executor accepted it into the
  registration, and neither caught that the detector is agnostic about
  the knob but not about death *timing* — the only thing a short window
  is sensitive to. Cost: one invalid arm caught at 8/20, zero data lost.
  The transfer question to ask is not "is the instrument validated?" but
  "was it validated where I am about to point it?"
- **A number transferred across mechanisms is the inherited-number trap
  across methods.** A quantity measured under one process and applied
  under another can be catastrophically wrong while every summary column
  still looks healthy. Check not just where a number came from, but what
  process produced it. (First met outside this line of work.)
- **A prediction guaranteed by arithmetic is not a prediction.** N2-main's
  "degeneracy ≥ 5" was implied by 31 parameters minus 15 features before
  any measurement; the agent annotated that in the JSON before the SVD, so
  the held verdict could not be sold as support. Check what a registered
  number is forced to be by dimension counting alone, and register past it.
- **Registering a criterion protects against forgetting it; only a test
  protects against not applying it.** dish-slow2a's registration stated the
  duty criterion in prose and the code enforced only rate and duration — a
  cell firing 77% of the time passed as a target hit, and a judgement was
  formed on it before the spikes-per-event check. A criterion that exists
  only in prose is a criterion that will be waived exactly once, silently.
- **A measurement taken for one purpose is worth re-reading for the
  question you are about to ask.** dish-slow's single-cell scan contained
  the whole failure (I_M cuts firing fivefold at low drive, a fifth at
  high) fifteen seconds in, read down the columns for a usable range and
  not along the rows for the mechanism.
- **Profile before optimising — a whole-system vector computed inside a
  per-element loop is correct, readable and quadratic.** Nothing about
  `syn_current(t, V)[i]` looks expensive; the symptom is cost growing with
  system size where the physics says it should not, which is what an
  n-sweep exists to expose. A single-size benchmark reads 1.00× and
  concludes the method does not help — which is exactly what happened the
  first time.
- **A normalised system's state is an exponentially weighted sum of past
  learning, not a running total.** Comparing accumulated updates against a
  normalised network's state requires discounting by the normaliser's time
  constant — worst exactly when the input is structured on that timescale,
  which is when anyone is looking. Same τ/dt constant phase-link made a
  tuning law from.
- **When two analyses of one dataset disagree, suspect the definitions
  before the data.** phase-link's registered reproduction check failed on
  first pass because one analysis differenced against constructed weights
  and the other against the stored post-warmup vector — two intervals, one
  name. Against the same definition they agree to every printed digit.
- **Score the trajectory, never the endpoint.** D1's segregation was reported
  absent; it was present at the peak and gone by the end. Written down two
  sections earlier and broken anyway.
- **A control must differ from the experiment in exactly one thing, and the
  yoked control must replay a *different* network.** Yoking a run to its own
  counts reproduces it stimulus for stimulus and controls nothing.
- **A wrong number that finds its explanation in the literature is the hardest
  kind to drop.** Distrust any biological explanation discovered after the
  number.
- **"It vanishes at zero window integral" is not an ablation** (N7). Use
  `rule_with_expected_change(0, lags)`.
- **A quantity that can move without plasticity proves nothing about
  plasticity.** One wiring in ten showed a sevenfold responsiveness rise
  with the weights immobile to 1 part in 10¹¹. Print the frozen arm beside
  any claim on such a quantity — and check that the frozen arm *can* move
  before calling it a control: `WeightNormalization` cannot move weights at
  all with the rule off, which made the frozen ladder a stronger control
  than designed, by luck.
- **A batch's endpoint must be snapshotted inside each program, and an
  equality test certifies only the quantities it compares.** The 16-check
  exactness test verified trajectories; the final weights were read after
  the longest program finished, so early-finishing dishes kept integrating
  and their endpoints drifted. "Verified exact" never extends past the
  columns the test actually reads.
- **A scan whose grid is coarser than the transition invents a step
  function** — the mirror image of the scan that stops before the effect.
  And knowing the rule does not tell you when the grid is fine enough:
  N8's executor named this rule and violated it in the same message,
  arguing from "the character of the first non-silent point" on a grid
  that had not resolved the transition. **The only test is refining until
  the answer stops moving**; theirs had not stopped. Corollary from the
  same incident: **never average across seeds before storing** — the mean
  0.196 at the sliver could be one bursting wiring averaged with one
  silent one; it is the float-keyed-json mistake in another dress.
  *Confirmed by the per-seed test, and promoted: **a mean over seeds can
  invent a regime no seed occupies.** At the sliver the per-seed values
  were 0.000, 0.392, 0.000, 0.481 — every wiring steps; the "ramp" was an
  average across dishes on opposite sides of their own thresholds. Three
  passes to get one number right (2 → 3 → 2), settled by neither argument
  but by the measurement both parties had pre-registered as decisive.*
  N8's drive curve went silent→saturated inside one 0.05 step and flat for
  the next 0.5; only a 10× finer bracket can say whether an intermediate
  regime exists. And its companion: **a single criterion saturates
  silently** — burst rate read 0.10 Hz (inside the registered band) at
  every saturated drive because one unbroken 10 s burst counts as one
  event; the pre-registered second criterion (active fraction < 0.2) is
  what exposed it.
- **Register the premise as a falsifiable branch, never as setup.** N8b's
  retune rested on "lower U drains slower, so bursts last longer" — the
  textbook reading — registered as testable rather than assumed. It was
  false (6/6 inverted): in a recurrent network the parameter that slows the
  drain also weakens the loop. Registering it is what turned a failed
  tuning run into a stated result.
- **Clear the question that decides interpretability before the question
  that decides the claim.** N8's step 2 asked "do these bursts have
  internal structure?" before "do they sort the weights?" — so when the
  sorting came back null, the null meant something. Asked in the other
  order, the same 9 minutes would have produced an uninterpretable result
  and a 35-minute protocol run on top of it.
- **A fix is not a fix until something reads it back.** A patch that adds
  storage can fail its anchor match silently and fix only the printing —
  which happened, on the per-seed rows, while the fix was being reported
  as done. Apply the equality-test standard to your own edits: after
  writing, assert the new field is present in the file and round-trip it
  through a real run. Companion error from the same hour, the
  a "re-analysis of stored spikes, zero compute" was planned without
  checking that spike times had ever been persisted — they had not; every
  calibration json held summaries only. **Verify what is on disk before
  pricing work against it.**
- **A registered contrast may not name a figure taken from a different
  script.** "Chains ≥ 8/10 where assemblies are 3/10" imported the 3/10
  from another reconstruction scored by a different verdict on a different
  quantity; measured in its own script the baseline was 7/10 and the
  contrast did not exist. Import the procedure or measure the baseline —
  never the number.
- **Batch by schedule, not just by rule.** A batch costs its longest
  program; mixing a 25 ms train with a 250 ms one pays 250 for all. "Same
  arm" must mean same rule, same budget, same plastic flag, *and* same
  stimulus schedule.
- **Every catch comes from a measurement, none from reasoning.** Six results
  in one day were clean, reportable and wrong — caught by a convergence
  curve, a trajectory, a target neither arm could reach, the original paper,
  a visited-g list, a spread guard. In each case the fix was to measure something the
  author believed they already knew. When a result is about to be written
  up, measure one such thing.
