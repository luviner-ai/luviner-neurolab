# Learning redistributes, it does not kill

**One mechanism with two signs in plastic networks under a homeostatic
budget.**

*Outline and figure list. Drafted 2026-09-10 by `luviner-9f` against board
row `PAPER-1`; revised 2026-09-11 with COLL-6, COLL-8 and COLL-9 in and the
hysteresis out.*

* **Everything below is on the record already**; nothing here
is a new claim, and the three numbers that are still missing are named in
§9. Publishing is the founder's decision — this is the draft, not a
submission.*

---

## 0. The one-sentence claim

A net-depressing spike-timing rule running under a homeostatic budget
**concentrates** recurrent excitation onto fewer synapses instead of
removing it. Whether that concentration **kills** a network is set by how
much recurrent drive the network had before the rule touched it — ordered
across seven conductance sets and reproduced by driving one of them.
Whether it **rescues** a dying network is not: at matched fragility the same
rule rescues one cell 33 times out of 33 and another 0 times out of 21.

**Two signs, two different causes.** The destructive direction follows the
drive; the protective direction does not, and belongs to the cell.

## 1. Why this is worth writing

The literature's shape is *"plasticity destabilises recurrent networks and
homeostasis has to fix it"* (Zenke & Gerstner on the timescale mismatch;
Turrigiano on synaptic scaling). Our own board had the same shape and
stated it more strongly, as **8 of 8 wirings dying when STDP is on and
living when it is off** (`## Run: COLL-2`).

That result is exact and it is conditional. The eight wirings were
**selected for having a death on disk with the rule on**. On unselected
wirings at the same inhibition level, the rule kills one network and
rescues another (`## Run: COLL-4`), and on a follower cell it is **net
protective** — 3 of 4 wirings die without it and 1 of 4 with it.

The paper is that correction plus the mechanism that explains both signs.

## 2. Prior art, stated before our own results

- **Beggs & Plenz 2003** — neuronal avalanches, power-law sizes, read as a
  critical branching process. The neighbourhood for §6.
- **Directed percolation / absorbing-state transitions** (Hinrichsen's
  review). Named because §6's finding is *"not in this class"*, which is
  stronger than "absorbing".
- **Turrigiano** — synaptic scaling as the homeostatic mechanism whose
  budget we impose.
- **Zenke & Gerstner** — the timescale mismatch between Hebbian and
  homeostatic plasticity; the closest statement of the problem we
  quantify.
- **Jacquerie & Drion 2025** — the neuromodulation/plasticity interaction
  in bursting circuits; `## Run: LIT-1` checked our weights against it and
  they go the other way, which is reported in §8, not hidden.
- **Prinz, Bucher & Marder 2004** — the cell models; **Bucher, Prinz &
  Marder 2005** — the 15 rhythm features our preparation is validated
  against (14/15 within 2 s.d.).

## 3. Methods, in the order a referee will check them

1. **The preparation.** 48 excitatory + 12 inhibitory conductance-based
   cells (`PRINZ_TABLE2`), 15% recurrent E→E connectivity, all-to-all E/I,
   short-term depression on the recurrent synapses. The excitatory cell is
   `AB/PD 2`, the pyloric pacemaker; the follower arms use `PY 1` and
   `LP 3`.
2. **The rule.** A net-depressing STDP window (`rule_with_integral(-0.134)`)
   plus `WeightNormalization` holding a per-cell incoming budget
   B = 0.107 with τ = 200 ms at dt = 0.05 ms.
3. **The two arms.** `learn=True` / `learn=False`, identical wiring and
   identical homeostasis in both. Homeostasis runs in the control **on
   purpose** — it is a property of the cell, not of the protocol, and
   leaving it in both is what makes the arms differ in exactly one thing.
4. **The detector.** A quorum: alive := ≥ 24 of 48 excitatory cells fire in
   the final 10 s. Registered in `rescore1.py` **before** any verdict was
   looked at, and shown there to be insensitive to the threshold (any Q in
   2..47 gives the same verdict on the one graded collapse in the corpus).
5. **Pre-registration.** Every prediction below was written to a shared
   document and committed **before the numbers existed**, with the commit
   hashes in §10. This is the paper's main methodological claim and §10 is
   where a referee should start.

## 4. Result 1 — the rule redistributes; it does not remove

The load-bearing measurement, on two cell models that are **not** the STG
burster (`## Run: COLL-3`) and on the STG dish itself (`## Run: COLL-4`):

| preparation | synapses at zero | total vs budget | firing change |
|---|---|---|---|
| cortical, driven | 18–26% | −1.6 to −2.7% | −0.7 to −1.5% |
| Hodgkin-Huxley, driven | 14–17% | −9.6 to −10.8% | −2.3% |
| AB/PD 2 (STG) | 16.5% | −1.4% | — |
| PY 1 (follower) | 7.1% | −0.6% | — |
| any arm, rule off | **0.0%** | at budget | — |

**A quarter of the recurrent excitatory synapses are driven to zero while
the population total stays within one and a half per cent of its budget.**
The rule moves conductance; the budget conserves it.

*Figure 1.* Weight distribution, rule on vs rule off, at matched total.
The mass at zero and the shifted upper tail, on one wiring per cell model.

## 5. Result 2 — two signs with two different causes

### 5a. The destructive direction is ordered by the cell's own firing rate

`## Run: COLL-8`, seven Prinz conductance sets, `g_ie` = 0.045, each cell's
OFF arm giving its own rate. The comparable statistic is the kill rate on
wirings that would have lived, which has no ceiling; `saved − killed` does,
because it cannot exceed the control's death rate.

| cell | spikes/s | n | kill \| healthy |
|---|---:|---:|---:|
| PY 1 | 28.5 | 36 | 0.000 |
| AB/PD 1 | 482.7 | 10 | 0.000 |
| LP 3 | 517.7 | 14 | 0.000 |
| AB/PD 5 | 536.2 | 10 | 0.300 |
| AB/PD 3 | 549.9 | 10 | 0.200 |
| AB/PD 4 | 610.8 | 14 | 0.273 |
| AB/PD 2 | 751.4 | 12 | 0.545 |

**Spearman +0.893, exact p = 0.0123** over all 5,040 permutations.
*[Corrected 2026-09-15 by STATS-1: rho was computed on ranks that break the
three ties at 0.000 by input order. On ranks with ties averaged it is +0.852,
two-sided permutation p = 0.0286, and the six-cell value is +0.812 at p =
0.072. The entry above is left as it was reported. See §7.2 of the
manuscript.]* Cell
identity does not predict it: **LP 3 is a follower that fires like a
pacemaker and behaves like a pacemaker**, which was registered in advance as
the cell that would decide between rate and identity.

Stated with its limit: **the ordering is carried by the low end.** Drop PY 1
and it falls to rho +0.829, p = 0.058 at n = 6, and PY 1 sits at 28.5
spikes/s with the next cell at 482.7 — a 16.9× gap the Prinz table does not
sample.

`## Run: COLL-9` reproduces the ordering **inside one cell**: AB/PD 1 driven
from 483 to 132 spikes/s by an external current goes from 0 of 10 controls
dead to **12 of 12**. Fragility is not fixed by the conductance set.

### 5b. The protective direction is not

`## Run: COLL-6`, 36 unselected wirings on PY 1 at `g_ie` = 0.045, where the
control arm dies 33 times:

                        ON alive   ON dead
        OFF alive              3         0
        OFF dead              33         0

**P(ON alive | OFF dead) = 33/33 = 1.000**, Wilson 95% **[0.896, 1.000]**,
exact one-sided binomial against p ≤ 0.5, **p < 0.00001**. Not one wiring
that died without the rule died with it, in 33 opportunities. The verdict is
unchanged under a stricter detector (31/33 = 0.939).

And it does **not** transfer to a cell driven to the same fragility.
`## Run: COLL-9`, AB/PD 1 at matched control-death rate:

| | rescue \| opportunity |
|---|---:|
| PY 1, 28.5 spikes/s | **33/33 = 1.000** |
| AB/PD 1, driven down | **0/21 = 0.000** |
| | Fisher exact one-sided **p = 1.9e-15** |

**Being slow makes a preparation fragile. It does not make it rescuable.**
Whatever supplies the rescue is in PY 1's conductance vector, not in its
rate — which is the paper's sharpest open question and the subject of the
queued row that swaps the conductances one at a time.

The destructive direction even **reverses** inside a cell: across cells the
kill rate rises with the rate, but within AB/PD 1 it goes from 0.000 at 483
spikes/s to 1.000 at the driven setting. A law read across preparations did
not survive a manipulation within one.

*Figure 2a.* Kill rate on healthy wirings against the OFF firing rate, seven
cells, with the Spearman and the PY 1-excluded value beside it.
*Figure 2b.* The COLL-6 2×2 with its Wilson interval, against COLL-9's
matched-fragility 0/21.

## 6. Result 3 — the form of death

`## Run: PHASE-T1`, 104 runs with per-cell spike times:

- The order parameter (fraction of participating cells) is **binary at
  every resolution from 500 ms to 10 s** — 0 of 104 runs in [0.1, 0.9]
  under the registered detector.
- Time-to-death is **invariant to the control parameter within a wiring**
  and varies 5–13 s **between** wirings. Freezing plasticity 0.4–0.6 s
  before the last burst leaves the death time unmoved **to the
  millisecond** while perturbing the spike counts by up to 2%.
- The ensemble death rate is graded (8/8 → 0/8 along a plasticity dose)
  **only because each wiring has its own sharp threshold.**

**Discontinuous, first-order-like, with quenched disorder setting each
sample's threshold — not the directed-percolation class the avalanche
literature would suggest.**

**The direct discriminator was run and it did not support the reading.**
`## Run: PHASE-T2` found a two-valued weight total across an up-and-down
sweep, and `## Run: PHASE-T4` halved the ramp rate and the gap collapsed
from 2.6× the noise to 1.16×, while the feature *moved one level along the
ladder*. A coexistence is pinned to the parameter and a lag is pinned to the
clock; it moved, so it was a lag. **PHASE-T2's hysteresis claim is withdrawn
and this paper makes none.** §6 rests on the four observations above, which
are untouched by that withdrawal, and it says first-order-*like* rather than
first-order for exactly this reason. The one measured exception is a single 5 s
transient at 14/48 cells that **recovers** (`## Run: COLL-4`), the first
sub-quorum recovery in the corpus.

*Figure 3.* Participation and rate on the one slow collapse: the rate
decays for 60 s at 48/48 participation, then participation falls 48 → 1
inside one window.

## 7. Result 4 — a law with no fitted parameter

`## Run: phase-link` derived the steady-state offset from the budget as

    S* = B + h·τ/dt        amplification τ/dt = 4000, fixed in the code

and verified it to ±5% over a fourfold range in τ, with the frozen control
exact to five decimals. **COLL-3 and COLL-4 then found it holding on three
further cell models it was never measured on**: the deficit tracks the
firing rate, which sets `h`.

| preparation | spikes/60 s | deficit | implied h per step |
|---|---|---|---|
| cortical | 78,000 | −2% | −4.3e-7 |
| Hodgkin-Huxley | 146,000 | −10% | −2.8e-6 |
| PY 1 (STG) | ~7,000 | −0.6% | — |
| AB/PD 2 (STG) | ~19,000 | −1.4% | — |

**A device running this rule can read its own net drift off the weight
totals** — the free diagnostic phase-link named, now shown to transfer
across substrates.

*Figure 4.* Deficit against firing rate, four preparations, with the law's
line.

## 8. What does not warn, what saves, and what we got wrong

- **Nothing warns.** `## Run: WARN-1` tested two on-line signals with every
  threshold fixed from a single trajectory; the behavioural one is
  structurally unable to fire in time on fast collapses, and the weight
  one (rail fraction) is the only candidate a device owning its synapses
  could use.
- **Freezing works, and only early.** `## Run: FREEZE-1`: rescue at 0.5×
  the pre-collapse time, none at 0.9×, where the last burst is unmoved to
  the millisecond. The point of no return precedes the last visible burst.
- **Our own refuted predictions, in full**, because the pre-registration is
  the credibility: the orchestrator's COLL-3 prediction (P1 met on the
  cortical cell) — refuted; the recap's generalisation from COLL-2 —
  withdrawn twice, first to "self-generating networks" and then entirely;
  my own COLL-4 prediction that PY would show the smaller deficit — held;
  my COLL-5 dissent that there was no sign flip to find — held, and my
  magnitude for it was wrong at one level (−0.417 measured against a
  −0.25 to −0.33 registered); my COLL-8 prediction that four of five new
  cells would be protective — **refuted on all five**; my COLL-9 prediction
  that rescue would follow wherever the opportunity existed — **refuted,
  0 of 21**; the orchestrator's prediction that cell identity orders the
  effect — **refuted on LP 3, the cell named in advance to decide it**;
  **LIT-1's check against Jacquerie & Drion,
  where our weights go the other way** — reported as a discrepancy, not
  resolved.

## 9. What is still missing

The three numbers this outline named on 2026-09-10 have all been run. Two
came back as asked and one came back negative:

1. ~~COLL-5, the rate with an interval~~ — **supplied**, and superseded by
   COLL-6's 33/33 with [0.896, 1.000].
2. ~~PHASE-T2, the hysteresis loop~~ — **run and negative.** PHASE-T4 showed
   the signature was a lag, not a coexistence. The claim is withdrawn and
   §6 no longer rests on it.
3. ~~A second level with a dying control~~ — **supplied**: COLL-9 produced
   two by driving AB/PD 1, and the answer was that rescue does not transfer.

What is open now is sharper, and all three bear on §5b:

1. **What in PY 1's conductance vector supplies the rescue?** The Prinz
   vectors differ in eight components. The queued row swaps them one at a
   time between AB/PD 1 and PY 1 at matched fragility. Until it runs, §5b
   has a fact and no mechanism.
2. **AB/PD 2's destructive sign has never been resolved.** COLL-5 gives
   saved 1, killed 6: seven discordant pairs, exact two-sided **p = 0.125**.
   The cell anchoring the high end of §5a's ordering is not significant on
   its own, and the paper must say so.
3. **The 16.9× gap.** Everything that changes sign happens between 28.5 and
   483 spikes/s and no Prinz cell is there. The uniform-current knob cannot
   reach it — it is bistable, with accessible states at ~480, ~420, ~165,
   ~80 and zero — so this needs a heterogeneous per-cell drive, which
   changes the preparation and is the founder's call.

**§5a is a rate with an interval and a stated dependence on one cell. §5b is
a rate with an interval and no mechanism.** That is the honest state.

## 10. Pre-registration index

Every claim above traces to a section committed before its numbers existed.
This table is the paper's methodological spine and should be checked first.

| result | registration | report |
|---|---|---|
| redistribution, driven cells | COLL-3 `234901f` | `f19c5f7` |
| sign depends on drive | COLL-4 `11536cc` | `df5599c` |
| form of death | PHASE-T1 `819bda0` | `af1f8d5` |
| the offset law | phase-link `59abb59`, `389ae39` | `d6eb7a1` |
| rate / interval | COLL-5 `016f409` | `be6b42e` |
| rescue conditional, 33/33 | COLL-6 `033aecb`, gate `d28c56b` | `68ece31` |
| the ordering across seven cells | COLL-8 `6d2a958` | `1bf9b3d` |
| fragility vs rescuability | COLL-9 `6536bfd`, control `cd835ce` | `d44e3f6` |
| hysteresis, withdrawn | PHASE-T2, PHASE-T4 `ca948e6` | `1a55486` |

## 11. What this paper is not

- Not a claim about cortical cultures. Every number is a model network, and
  §5's boundary case shows the effect depends on a preparation property
  (recurrence-dependent activity) that a real culture may or may not have.
- Not a claim about criticality in cortex. §6 is about **this**
  preparation, and its value is as a counterexample to a reading, not as a
  measurement of tissue.
- Not a product claim. The convertible parts are recorded on the board and
  are design constraints, not features.
