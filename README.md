<p align="center"><img src="docs/assets/banner.svg" alt="PipeCompare" width="640"></p>

**Data-driven comparison of EEG preprocessing pipelines in EEGLAB.**

[![MATLAB tests](https://github.com/abwoo/PipeCompare/actions/workflows/matlab-tests.yml/badge.svg)](https://github.com/abwoo/PipeCompare/actions/workflows/matlab-tests.yml)
![Version](https://img.shields.io/badge/version-0.9.5-2f6fed)
![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-orange)
![EEGLAB](https://img.shields.io/badge/EEGLAB-2024.2%2B-blueviolet)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

PipeCompare is an [EEGLAB](https://github.com/sccn/eeglab) plugin that chooses how to preprocess
an EEG recording, based on the recording itself. It compares complete preprocessing pipelines
(filters, reference, bad channels, ICA with ICLabel, epoch rejection) and recommends the one that
measures your ERP or band power most precisely, without distorting it.

## Download

<p align="center">
  <a href="https://github.com/abwoo/PipeCompare/releases/latest"><b>Download the latest PipeCompare (zip)</b></a>
  &nbsp;·&nbsp;
  <a href="https://github.com/abwoo/PipeCompare/releases">All versions</a>
</p>

PipeCompare is used like ERPLAB Studio: installed once next to EEGLAB, then opened with one
command, without EEGLAB's own window.

1. Install [EEGLAB](https://sccn.ucsd.edu/eeglab/download.php) (2022 or later) and add its folder to
   the MATLAB path (**Home > Set Path > Add Folder…**, or `addpath('/path/to/eeglab')`).
2. Download `PipeCompare<version>.zip` from the
   [latest release](https://github.com/abwoo/PipeCompare/releases/latest) (under *Assets*).
3. Unzip it into the `plugins` folder of EEGLAB, so that you have
   `eeglab/plugins/PipeCompare<version>/`.
4. Save the path, once: in MATLAB's Command Window type `eeglab`, then, when EEGLAB has started,
   `savepath` (MATLAB then finds PipeCompare and EEGLAB every time it starts; close EEGLAB's window
   afterwards).
5. From then on, type `pcomp` in MATLAB's Command Window. The PipeCompare window opens; it
   starts EEGLAB's functions in the background, without EEGLAB's window.

If `pcomp` is not found, step 4 was skipped or the path could not be saved (on a shared
computer, MATLAB may not be allowed to write it): type `eeglab nogui` first, then `pcomp`.

PipeCompare is also listed in EEGLAB's extension manager (**File > Manage EEGLAB extensions**),
which downloads the same zip. The details and the other ways are under
[Installation](#installation). In every PipeCompare window, **Help > PipeCompare on GitHub** opens
this page, and **Help > Download the latest version** the latest release.

## Contents

- [Download](#download)
- [What PipeCompare does](#what-pipecompare-does)
- [What happens to your data](#what-happens-to-your-data)
- [Features](#features)
- [Where PipeCompare fits in an analysis](#where-pipecompare-fits-in-an-analysis)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Scripting interface](#scripting-interface)
- [How a pipeline is chosen](#how-a-pipeline-is-chosen)
- [Outputs](#outputs)
- [Limitations](#limitations)
- [Testing](#testing)
- [Documentation](#documentation)
- [Citation](#citation)
- [License](#license)

## What PipeCompare does

Before an ERP amplitude or a band power can be measured, the raw recording has to be cleaned: it
is filtered, re-referenced, bad channels are found and repaired, eye and muscle artifacts are
removed with ICA, the data are cut into epochs and noisy epochs are rejected. Every one of these
steps has settings, such as the high-pass cutoff, the ICLabel threshold for removing a component
or the amplitude limit for rejecting an epoch. One complete set of steps and settings is a
*pipeline*. Published studies use many different pipelines, and the choice changes how much noise
is left in the final measurement. In practice the settings are usually copied from an earlier
paper or kept out of habit, without checking how well they work on the data at hand.

PipeCompare makes this choice from the data, in these steps:

1. **It reads the dataset open in EEGLAB**, including what has already been done to it (from
   `EEG.history`), so earlier processing is neither repeated nor ignored.
2. **You say what will be measured**: an ERP component (mean amplitude, peak amplitude or peak
   latency in a time window, at the electrodes you choose, time-locked to the events you choose),
   or the power in a frequency band. Presets cover the ERP CORE components and the classic
   frequency bands.
3. **It builds every pipeline** from the steps and settings left open. The *Standard* recipe, for
   example, compares up to 12 combinations of high-pass and low-pass filters, and in each of them
   finds the best ICLabel threshold and epoch-rejection limit from the data; you can also decide
   the steps and values yourself.
4. **It runs each pipeline** with EEGLAB's own functions, on your data.
5. **It removes pipelines that harm the data.** A pipeline is excluded when it keeps too few
   trials, interpolates too many channels, or changes the brain signal. To check the last point,
   a noise-free copy carrying a known artificial signal replays the fitted pipeline decisions;
   if it comes out smaller, shifted or reshaped, the pipeline is excluded.
6. **It ranks the rest by precision.** The measure is the standardized measurement error (SME;
   Luck et al., 2021): how much your averaged value would vary if the experiment were repeated.
   It is corrected for how much each pipeline shrinks the signal, so a pipeline cannot look good
   just by making everything smaller. A bootstrap test then finds the pipelines that cannot be
   told apart from the best.
7. **It recommends one pipeline** and explains why. Each pipeline's error is first estimated
   together with the others that share its settings (how much a 0.5 Hz high-pass helps is seen in
   every pipeline with it), which picks the truly better pipeline more often than its own error
   alone. Among pipelines the data show to be equally good, the one that keeps the most trials and
   changes the signal least is used. You can adopt it as a new EEGLAB dataset,
   already epoched and cleaned and ready to average, or save it as a MATLAB script to process your
   other recordings the same way.

The comparison never uses an experimental effect (condition differences, p-values), so selecting a
pipeline does not bias the statistical test you run afterwards.

PipeCompare is useful when you start working with a new dataset, paradigm or recording setup and
want preprocessing settings with a documented, data-based reason, or want to check whether your
usual settings suit these data. It works on one recording at a time, and covers ERP and band-power
measures; time-frequency, connectivity and source analysis are not covered (see
[Limitations](#limitations)).

## What happens to your data

With the *Standard* recipe, every pipeline runs the steps below in this order. Steps marked
*compared* are tried with several settings; the others are done once, in the same way for every
pipeline, so that the comparison is about the settings that matter.

1. **Channels removed earlier are put back** (only with the average reference). If you deleted EEG
   channels before PipeCompare, they are interpolated back first, so that the average is taken
   over the whole montage.
2. **Bad channels are found and repaired.** A channel is marked bad when, more than 5 standard
   deviations away from the other channels, its signal is much spikier (kurtosis), its values are
   unusually distributed (joint probability), or it carries **noise of its own** that the other
   channels do not explain. The last test is PipeCompare's: brain activity and artifacts reach
   neighbouring electrodes too, so the other channels predict most of a good channel (a
   least-squares fit of each channel on all the others); a poorly connected or broken electrode
   adds noise they cannot predict. Kurtosis and joint probability look only at the shape of a
   channel's values, not their size, so a channel that is just two to five times noisier passes
   them (in simulations they found almost none of these; the new test finds most of them, with few
   few false alarms), and one such channel can cost half the epochs at a 100 µV rejection limit.
   Size differences the channels share (strong alpha, blinks, a nearby reference) are predicted,
   so they do not count, and the test needs no channel locations. Kurtosis and joint probability
   are measured without the stretches with blinks (where an eye channel, or the blink's own
   course read from all channels, is beyond 5 robust SDs, or where a hidden Markov model of the
   eye channels puts a blink, see [Blinks and bad stretches found by a hidden Markov
   model](#blinks-and-bad-stretches-found-by-a-hidden-markov-model); 200 ms on either side): blinks make the frontal channels' values heavy-tailed, and in
   simulations several good frontal channels per recording were otherwise marked bad, the
   channels on which ICA sees the blinks best. It is left out (and said in the
   Command Window) when the channels are linearly dependent, as after an average reference. A
   **flat channel** (no signal, as when
   an electrode is disconnected: its spread is zero or under 1 % of the median spread of the
   channels tested) is always bad; it is set aside before the measures, which cannot be
   computed on a constant and would otherwise be thrown off for all the other channels too. The
   test looks at a copy of the data
   high-passed at 1 Hz, so slow drifts are not mistaken for bad channels, and low-passed at 40 Hz
   (for band power, 5 Hz above the highest band), the band the analysis keeps, so that line noise
   alone, which the low-pass removes, does not make a channel bad. Bad channels are then
   interpolated from their neighbours (spherical interpolation), so the dataset keeps all its
   channels, **unless interpolating one would make the measure worse**: PipeCompare computes the
   total error (see [How a pipeline is chosen](#how-a-pipeline-is-chosen)) with each bad channel
   interpolated or left as it is, through the reference, epochs, baseline and epoch rejection that
   follow (ICA and the filters compared are left out of this check, the filters replaced by the
   1–40 Hz band the channels were tested on), and leaves a bad channel as it is only when
   interpolating it would raise the error measurably. This matters most for a measured electrode:
   interpolating Pz from noisy neighbours brings their noise in and loses the P3 at Pz. The
   Command Window and the result name the channels left as they are. In scripts this is the
   badchannels option `'repair', 'error'` (`'all'`, the default outside the simple mode and the
   panel, interpolates every bad channel). On continuous data a channel can also be **bad for
   only part of the recording** (an electrode that came loose after 20 minutes): such stretches,
   10 s or longer, are found by a hidden Markov model of each channel's noise the others do not
   explain and interpolated there only (see [Blinks and bad stretches found by a hidden Markov
   model](#blinks-and-bad-stretches-found-by-a-hidden-markov-model); badchannels option
   `'stretches', true`, on in the simple mode and the panel). EOG, ECG and EMG channels (by channel type, or by name such as VEOG or ECG1) are never
   tested. A pipeline that interpolates more than 20 % of the channels is excluded.
3. **Re-reference** (if you chose the average reference). It comes after the bad channels are
   repaired, so a bad channel cannot spread its noise into every other channel; non-EEG channels
   are left out of the average. Data that were already average-referenced are averaged again,
   which removes the repaired channels' share of the earlier average.
4. **ICA is computed once** (extended Infomax, `pop_runica`) on a copy high-passed at 1 Hz, which
   gives a cleaner decomposition. The copy also leaves out the stretches of recording that are far
   noisier than the rest (movement, electrode bursts; see
   [ICA on clean data](#ica-on-clean-data)), so that they do not claim components of their own.
   The resulting decomposition is applied to all of the data. All pipelines share this one
   decomposition, which saves the largest part of the computing time.
5. **High-pass filter**, *compared*: 0.1, 0.3, 0.5 and 1 Hz.
6. **Low-pass filter**, *compared*: 20, 30 and 40 Hz. After the comparison, the cutoffs of the
   recommended pipeline are also tried on a finer grid between these values (0.1, 0.15, … 1 Hz and
   20, 22, … 40 Hz), and when a pair is better, that pipeline is run as well and ranked with the
   others (see [Filter cutoffs fine-tuned from the data](#filter-cutoffs-fine-tuned-from-the-data)).
7. **Artifact components are removed with ICLabel**, *threshold chosen from the data*: a
   component is removed when ICLabel gives it a probability above a threshold of being eye,
   muscle, heart, line noise or channel noise, and the threshold is found in each pipeline from its
   own data together with the rejection limit of step 9: every threshold from 0.5 up is tried, and
   the highest one that measures as precisely as the best one is used (see
   [ICLabel threshold chosen from the data](#iclabel-threshold-chosen-from-the-data)). For peak
   measures, or when epochs are repaired, the thresholds 0.7, 0.8 and 0.9 are compared instead.
8. **Epochs are cut** around the events you chose, and the **baseline** is subtracted.
9. **Noisy epochs are rejected**, *limit chosen from the data*: an epoch is dropped when any EEG
   channel exceeds ±*x* µV (non-EEG channels are ignored), and *x* is found in each pipeline from its
   own data: every limit is tried, not a few, and the highest limit whose SME is shown to be within
   5 % of the best one's is used (see [Rejection limit chosen from the data](#rejection-limit-chosen-from-the-data)).
   For peak measures (peak amplitude or latency, advanced panel) the limits 75, 100 and 150 µV
   are compared instead.
10. **Epochs weighted by their noise**, *compared with equal weights*: in half of the pipelines
   the average counts every remaining epoch the same, in the other half a noisier epoch counts
   less, and no epoch is removed for it (see [Epochs weighted by their noise](#epochs-weighted-by-their-noise)).
   ERP only.

Two more options of the rejection can be ticked; they are not part of *Standard*.
**Peak-to-peak limits** (ERP CORE's method): instead of testing whether any value is beyond ±*x*
µV, the rejection takes the highest minus the lowest value within a moving 200 ms window; the
limit is again chosen from the data (see [Peak-to-peak rejection](#peak-to-peak-rejection)).
**Repairing epochs instead of rejecting them**: an epoch in which only 1 to 3 channels exceed the rejection limit is kept, and those
channels are interpolated from the other channels within that epoch only; an epoch with more
channels over the limit is still rejected. With repair, the limit is not chosen from the data
(that choice assumes the epochs over the limit are removed): the fixed limits 75, 100 and 150 µV
(100, 150 and 200 µV peak-to-peak) and the ICLabel thresholds 0.7, 0.8 and 0.9 are compared, which
makes nine times as many pipelines. See
[Repairing epochs with a few bad channels](#repairing-epochs-with-a-few-bad-channels).
An option of ICA can also be ticked, not part of *Standard*:
**choosing the components by their effect on the noise**. Every pipeline is then also run with
the components to remove chosen not by what ICLabel says they look like, but by how much noise
each one adds compared with how much of the signal it carries; ICLabel is only used to protect the
components it takes for brain activity, and components that follow the stimulus are kept unless
they follow an eye channel. Both ways are compared, which doubles the number of pipelines. See
[ICA components chosen by their effect on the noise](#ica-components-chosen-by-their-effect-on-the-noise).

That makes 4 × 3 × 2 = 24 pipelines for an ERP, each with its own best ICLabel threshold and rejection limit. The *Filters only* recipe compares only
steps 5 and 6. For band power, the data are cut into 2 s segments instead of epochs, and only the
filter edges nearest the band, outside it, are used (a filter outside a band does not change its
power), which leaves 1 pipeline (its ICLabel threshold and rejection limit still chosen from the
data). In the advanced panel or a script you can change every list
above, change the order, allow a step to be skipped, and add other steps: line-noise removal (the
50 or 60 Hz mains frequency is detected from the recording), ASR (artifact subspace
reconstruction), resampling, rejection by joint probability or kurtosis, a reference to chosen
channels, or any operation from EEGLAB's menus, plugins included.

Along the way PipeCompare also takes care of the following, so you do not have to:

- **It names flat channels** (no signal, e.g. a disconnected electrode) before the run, in the
  dialog's first line and the panel's data warnings, and the bad-channel step always marks them
  bad. Without that step, remove or interpolate them yourself before an average reference.
- **It checks that the channel labels fit the electrode positions** before the run, and names
  channels whose signal does not look like their neighbours' (see [Montage check](#montage-check)).
  This is a warning only.
- **It tells you when ICA did nothing useful**: when ICLabel recognised almost none of the
  components, or no pipeline removed any component, and when the recording was too short for ICA
  (see [When ICA does nothing](#when-ica-does-nothing)).
- **It tells you when an electrode you measure was rebuilt** from its neighbours by the
  recommended pipeline (as a bad channel). A rebuilt electrode has less noise than a real one, so
  that pipeline's noise looks smaller than it is (see [Reading the result](#reading-the-result)).
- **It tells you when one channel caused most rejected epochs**, which usually means it is a bad
  channel the detection missed (see [Reading the result](#reading-the-result)).

- **It respects what was done before.** Processing already applied (read from `EEG.history`) is
  listed and not repeated. A filter the data already have is kept as one of the choices
  (e.g. *low-pass as in the data (30 Hz)*) instead of forcing a stricter one, and you are warned
  when that earlier filter alone already distorts the signal you measure. When ICA was already run
  and components removed, ICA is not compared again.
- **It leaves out what the data cannot support, and says why.** Without channel locations,
  interpolation and ICLabel are left out; without ICLabel installed, ICA is left out; on data that
  are already epoched, filters are not compared (they must run before epoching) and the data keep
  their own epochs as long as they hold the measurement window. A condition with too few events
  (fewer than 10) is flagged before the run.
- **It checks units.** Data stored in volts are recognised from their amplitude and compared in µV.
- **It measures lateralized components correctly.** N2pc and LRP are scored as ERP CORE does,
  contralateral minus ipsilateral, from the event types you give for each side.
- **It handles long runs.** A progress window shows the step running and the time so far and left;
  *Stop* keeps the pipelines already finished. A search can be saved to a checkpoint folder and
  resumed, and can run in parallel with the Parallel Computing Toolbox.
- **It keeps everything reproducible.** The adopted dataset's `EEG.history` contains every EEGLAB
  command that produced it, and *Save script…* writes the same steps as a MATLAB function that makes
  the data-driven decisions (bad channels, components, rejected epochs) anew on each recording you
  run it on.

## Features

- **One window for everything.** `pcomp` (or *Tools > PipeCompare > PipeCompare
  window…*) opens the PipeCompare window: open the recording, add channel locations, choose what
  to compare (simple or advanced), follow the progress, read the result and save the cleaned data
  and the script, all in one window and without EEGLAB's main window. The EEGLAB menus stay.
- **Works on the dataset you have open.** Reads the current EEGLAB state (continuous or epoched
  data, channel locations, ICA decomposition, reference, sampling rate) and parses `EEG.history`,
  so earlier processing is taken into account.
- **ERP and band-power measures.** Mean amplitude, peak amplitude and peak latency for event-related
  data. Log band power for continuous recordings such as resting state, estimated with multitaper
  spectra (3 Slepian tapers per segment, fewer for a band narrower than 4 Hz in 1 s segments),
  which in simulation lowers the error of the band power by about a quarter compared with one
  Hann taper.
- **Presets for common analyses.** ERP CORE parameters (Kappenman et al., 2021) for N170, MMN,
  N2pc, N400, P3, LRP and ERN, and standard delta, theta, alpha and beta bands.
- **User-defined search space.** Choose the steps and their order, fix some values and let others
  vary, allow steps to be skipped. Any EEGLAB menu operation that takes and returns the dataset,
  plugins included, can be added as a step and configured in its own EEGLAB dialog.
- **Exhaustive search.** Every admissible pipeline is run, with shared leading steps computed only
  once. The search is never sampled or truncated.
- **Signal-preservation check.** A known signal is carried through every pipeline with the same
  data-driven decisions as the real data (bad channels, ICA, rejected epochs, removed components,
  ASR reconstructions); steps whose decisions cannot be replayed are re-run on it and flagged (see
  [Limitations](#limitations)). Pipelines that distort it are excluded.
- **Exact best rejection limit and ICLabel threshold.** Within each pipeline, the epoch-rejection
  limit and the ICLabel threshold are chosen from the data together, among all values that give
  different results, not from a short list.
- **Filter cutoffs fine-tuned from the data.** Between the compared cutoffs, a finer grid is
  scored exactly from the measured electrodes alone, with the recommended pipeline's decisions held
  fixed; a better pair runs as one more pipeline.
- **ICA components chosen by their effect on the noise.** Optionally, the components to remove are
  chosen by how much removing each lowers the SME (never one ICLabel takes for brain activity),
  compared with ICLabel's threshold.
- **Signal check with your recording's own waveform.** In the simple mode, the known signal takes
  the shape of your data's own average waveform at each measured site, where that average is
  reliable, instead of a standard bump.
- **Epochs weighted by their noise.** Instead of only keeping or dropping epochs, the average can
  count each epoch by how clean it is (inverse-variance weights), compared with equal weights.
- **How stable the recommendation is.** The result says how often the same pipeline would be
  recommended if the trials were drawn again.
- **One pipeline for a group, and each recording's own.** The *6. Group* page of the PipeCompare
  window (or `pipecompare.Group.run`) compares the pipelines on every recording and recommends the
  one that is best for the group as a whole, and also each recording's own pipeline, estimated
  together with the others (two-level empirical Bayes).
- **Hidden Markov models instead of fixed thresholds.** Blinks, and channels bad for only part of
  the recording, are found by two-state hidden Markov models fitted to each recording (EM and
  Viterbi), which find weak blinks a fixed threshold misses and interpolate a channel only where
  it is bad.
- **Eye channels regressed out when ICA cannot remove the blinks.** The recommended pipeline is
  also tried with the EOG channels regressed out of the EEG, fitted on the eye activity without the
  ERP, and kept only when its total error is lower.
- **Statistically controlled ranking.** A paired bootstrap with simultaneous intervals identifies
  the pipelines that cannot be distinguished from the best, controlling the error rate across all
  candidates.
- **Reproducible output.** Each candidate's full list of EEGLAB commands, an exportable script, and
  adoption as a new dataset whose `EEG.history` reproduces it.
- **Practical for long searches.** A progress window that can stop a search and keep the finished
  pipelines, checkpointing and resume, and optional parallel execution.

## Where PipeCompare fits in an analysis

PipeCompare replaces the preprocessing part of a single-recording EEGLAB analysis:

1. **Before PipeCompare** (done by you): import the recording, add channel locations
   (*Edit > Channel locations*, or both on the *Import* page of the PipeCompare window), and, if
   needed, remove channels that are not EEG or that you know are broken.
2. **PipeCompare**: filtering, reference, bad-channel detection and interpolation, ICA with
   ICLabel, epoching, baseline correction and epoch rejection. These are the steps it compares, so
   they should not be applied beforehand.
3. **After PipeCompare**: average the epochs into ERPs (or compute spectra for band power) and
   measure. The adopted dataset is already filtered, re-referenced and cleaned, so these steps are
   not repeated.

Processing that was already applied to the data is read from `EEG.history` and taken into account:
the dialog lists it, does not compare it again, and checks whether it already distorts the
measured signal. The cleanest comparison starts from the raw continuous recording, because then
every step is part of the comparison.

PipeCompare compares pipelines on one recording at a time. To process a whole study the same way,
choose the pipeline on a representative recording and run the saved script (*Save script…*) on
every recording; the script makes the data-driven decisions (bad channels, components, rejected
epochs) from each recording's own data.

## Requirements

| Component | Requirement |
|---|---|
| MATLAB | R2021b or later (see the tested versions below). GNU Octave is not supported: the interface uses `uifigure`. |
| [EEGLAB](https://github.com/sccn/eeglab) | 2022.0 or later. The PipeCompare window needs EEGLAB on the MATLAB path but not its main window (it starts EEGLAB without it); the PipeCompare menu and *Add EEGLAB menu step…* need EEGLAB's main window; scripts also run after `eeglab nogui` |
| EEGLAB plugins | firfilt (bundled with EEGLAB); ICLabel for IC removal; clean_rawdata for ASR |
| Optional | Parallel Computing Toolbox, for parallel execution; Signal Processing Toolbox, for ASR at sampling rates other than 100, 128, 200, 256, 300, 500 and 512 Hz and for faster low-frequency high-pass filtering |

MATLAB R2026a with EEGLAB 2026.0.0 is the current tested combination
([docs/COVERAGE.md](docs/COVERAGE.md)). From PR #50 onward, continuous integration,
installation checks and benchmarks use this pinned pair. Numerical, native EEGLAB,
GUI and end-to-end suites are retained. Earlier releases were tested by the historical
compatibility matrix through PR #49; they are no longer checked on each change.

Two known problems of older EEGLAB releases, which PipeCompare works around or reports:

- EEGLAB 2024.0 and older, on recent MATLAB releases, stop with "Colon operands must be real
  scalars" when cutting time ranges out of continuous data (EEGLAB's own *Edit > Select data*).
  PipeCompare makes these cuts itself when that happens, with the same result.
- EEGLAB 2024.0 comes with ICLabel 1.5, which classifies components wrongly (its scalp maps are
  rotated, so eye components are mostly missed). When IC removal runs with that ICLabel, the
  Command Window shows a warning; update ICLabel to 1.6 or later in *File > Manage EEGLAB
  extensions* before comparing pipelines with IC removal.

PipeCompare reads the current dataset from EEGLAB's base-workspace variables `EEG`, `ALLEEG`
and `CURRENTSET`, and stores adopted pipelines there.

## Installation

PipeCompare is listed in EEGLAB's extension manager
([list of EEGLAB extensions](https://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_list_all.php)), so the simplest way to install it is from
inside EEGLAB. The other two ways are for computers without internet access in MATLAB and for
working with the source code.

Keep only one copy of PipeCompare in `eeglab/plugins/`. EEGLAB loads every plugin folder it finds,
so two copies (for example one from the extension manager and one unzipped by hand) both add
their menus and their functions shadow each other.

**From EEGLAB's extension manager (recommended)**

1. Start EEGLAB and choose **File > Manage EEGLAB extensions**.
2. Find **PipeCompare** in the list (typing its name in the search field narrows the list), tick
   it and press **Install/Update**.
3. The menu **Tools > PipeCompare** appears; if it does not, restart EEGLAB.

EEGLAB downloads the release zip and places it in `eeglab/plugins/PipeCompare<version>/`. When a
newer version is listed, the same window offers it, and installing it replaces the old folder.

**From a GitHub release**

1. Download `PipeCompare<version>.zip` from the
   [Releases](https://github.com/abwoo/PipeCompare/releases) page.
2. Unzip it into `eeglab/plugins/`. The folder must be named `PipeCompare` or
   `PipeCompare<version>`, because EEGLAB takes the plugin's name and version from the folder name.
3. Start or restart EEGLAB. The menu **Tools > PipeCompare** appears.

**From source**

```bash
git clone https://github.com/abwoo/PipeCompare.git
```

Either copy the cloned folder into `eeglab/plugins/` as above, or keep it elsewhere and register
it with the running EEGLAB:

```matlab
eeglab                                    % start EEGLAB first
pipecompare_setup                         % run from the PipeCompare folder
```

`pipecompare_setup` adds the menu to the running EEGLAB session only. Run it again after each
`eeglab` call, because EEGLAB rebuilds its menus.

## Quick start

There are three ways to start: the PipeCompare window (everything in one window), EEGLAB's menu
(the same dialogs on the dataset open in EEGLAB), or the command line. They use the same
comparison and give the same results.

### In the PipeCompare window

Type `pcomp` (or `pipecompare.start`, the same) in MATLAB's Command Window, or choose
**Tools > PipeCompare > PipeCompare window…** in EEGLAB. PipeCompare and EEGLAB must be on the
MATLAB path: once the path has been saved after installing (see [Download](#download): start
EEGLAB, then `savepath`; this works the same when PipeCompare was installed from EEGLAB's extension
manager), `pcomp` works in every MATLAB session. EEGLAB's main window does not have to be
open: when EEGLAB has not been started yet, the window starts its functions without its main
window. Without a saved path, type `eeglab nogui; pcomp`. Calling it again brings the open
window to the front.

The steps are listed on the left, in the order you do them; clicking one shows its page. Under
them, **Datasets** lists the datasets open (the recording, each pipeline used, in EEGLAB too);
clicking one makes it the current dataset, as EEGLAB's *Datasets* menu does. The window never
needs EEGLAB's own dialogs; only the EEGLAB buttons of the Advanced tab open them, for those who
want them.

1. **Import.** **Open data file…** reads a recording with EEGLAB's own import for its format:
   `.set` directly, EDF, BDF and GDF through EEGLAB's BIOSIG import, other formats (BrainVision,
   Neuroscan, EGI and more) through EEGLAB's FILE-IO import. If EEGLAB lacks the import for a
   format, it offers to install it. The recording becomes a new EEGLAB dataset, the current one,
   so EEGLAB's menus and the Command Window see it too; a dataset that is already current in
   EEGLAB is used as it is. The page describes the dataset (channels, sampling rate, length or
   epochs, event types) and says which channels have positions. **Channel locations** adds them:
   **Look up standard positions** finds each channel's position in the standard 10-05 system by
   its label (EEGLAB's channel lookup, with the standard file of its DIPFIT plugin); **From a
   file…** loads the locations file of your cap. Channels whose labels are not standard (for
   example EOG channels) stay without a position and are named. The table lists every channel
   with its label, type and whether it has a position: correct a label there (for example a
   misspelled one, then look up the positions again) or a type (EOG, ECG, EMG or REF for channels
   that are not scalp EEG, which the comparison then leaves out of bad-channel detection, the
   average reference and epoch rejection), and press **Apply label and type changes**. Each change
   replaces the current dataset and is recorded in its `EEG.history`. Interpolating bad channels
   and ICLabel need the positions. On the right, a plot of the data: continuous data as one trace
   per channel, 10 s at a time, with a slider to move through the recording (a flat or very noisy
   channel stands out); epoched data (for example a pipeline already used) as the average of the
   epochs of each event type, at the channel chosen.
2. **Settings**, in two tabs. **Simple** holds the simple mode's choices, the same as the dialog
   described [below](#from-the-eeglab-menu) (what to measure, the event types, the steps), with
   **Run** and **Advanced…** but no Cancel; **Electrodes…** opens the window's own list of the
   dataset's channels. **Advanced** holds the full advanced panel (every
   control is described in [docs/PANEL.md](docs/PANEL.md)); it is opened the first time the tab
   is shown. Its EEGLAB buttons (*Configure in EEGLAB…*, *EEGLAB pop_epoch…*, *View ERP
   (EEGLAB)…* and so on) open EEGLAB's dialogs and are optional: every setting can also be typed or
   set with *Edit values…*. *Add EEGLAB menu step…* lists EEGLAB's menus, so it is greyed out while
   EEGLAB's main window is closed (type `eeglab` to open it).
   **Advanced…** in the Simple tab opens a new panel in the Advanced tab with the choices made
   so far, and the Simple tab keeps them. When the current dataset changes (a file opened,
   locations added, a pipeline used, or a change made in EEGLAB), the Simple tab is made again
   for the new data the next time a page is shown.
3. **Run.** **Run** (Simple tab) or **Run search** (Advanced tab) shows this page: how many
   pipelines are done, the step running now, the time so far and the time left, and **Stop**,
   which ends the search and keeps the pipelines already finished, as the progress window of the
   menu does. While a comparison runs, the other pages wait. If it stops with an error, the
   message is shown and the window goes back to Settings.
4. **Results.** When the comparison is done, this page shows the result in plain words: which
   pipeline to use and why, the table of pipelines, the steps of the selected pipeline, **Show all
   pipelines**, **Use this pipeline** and **Save script…**, as the result window of the menu (see
   [Reading the result](#reading-the-result)). It also shows the result of a search run in the
   Advanced tab, whose own table keeps the details. **Use this pipeline** makes the pipeline the
   current dataset and goes to Export. On the right, **Plot selected pipeline** runs the selected
   pipeline (or the recommended one) again from the start of the comparison, reusing its ICA, and
   plots the average of each condition at the measured electrodes with the measurement window
   marked (band power: the spectrum at the electrodes, with the band marked). The lines of each
   pipeline plotted stay, so pipelines can be compared by eye; **Clear plot** starts again. For
   N2pc and LRP the plot averages the electrodes, while the comparison scores contralateral minus
   ipsilateral.
5. **Export.** **Save dataset (.set)…** saves the current dataset, for example the pipeline just
   used, as an EEGLAB `.set` file; until then it is in memory only. After **Use this pipeline**,
   the page says what to do next: the data are cut into epochs (segments for band power) and
   cleaned, ready to average into ERPs or to compute band power; do not filter, re-reference or
   reject epochs again. **Save script…** writes the recommended pipeline (or the one selected on
   Results) as a MATLAB function that runs the same steps on your other recordings.
6. **Group.** Compares the pipelines on several recordings at once (see
   [One pipeline for a group of recordings](#one-pipeline-for-a-group-of-recordings)). The list
   holds the datasets open in EEGLAB; **Add files…** adds `.set` files (they are read one at a
   time, not kept open). Select the recordings (all are selected at first) and press **Run for
   these recordings**: each is compared with the Simple tab's settings. The page then says which
   pipeline to use for the whole group and, in the table, each recording's own pipeline and
   whether it is the group's. **Save group script…** writes the group's pipeline as a MATLAB
   function, **Save script of selected recording…** the pipeline of the recording selected in the
   table.

The window and EEGLAB share the data: the window works on the dataset that is current in EEGLAB
and stores there what it opens or makes, so you can switch between the window, EEGLAB's menus and
the Command Window at any time. Closing the window stops a comparison that is running (the
pipelines already finished are kept in `pipecompare_result`).

The **Help** menu of the window opens this GitHub page (**PipeCompare on GitHub**), this guide
(**User guide (README)**), the latest release to download (**Download the latest version**) and
the page to report a problem (**Report a problem**), in your web browser; **About PipeCompare**
shows the version. The simple dialog, the results window and the advanced panel opened from
EEGLAB's menu have the same Help menu, and EEGLAB's menu has **Tools > PipeCompare > PipeCompare
on GitHub…**.

### From the EEGLAB menu

1. Load a dataset in EEGLAB.
2. Open **Tools > PipeCompare > Compare pipelines…**
3. Choose what to measure: an ERP component (with the event types it is time-locked to) or a
   frequency band, or your own time window or band with the electrodes you pick. The list only
   offers what the data support. The electrodes and time window of each ERP component are not
   taken from your data but from the published ERP CORE conventions (Kappenman et al., 2021),
   for example P3 at Pz, 300–600 ms. You can change a component's electrodes with
   **Electrodes…** next to its time window (for example Pz, CPz and POz for the P3); the window
   stays the convention's. N2pc and LRP keep their electrode pair, since they are scored as the
   difference between the two sides. When several electrodes are chosen, PipeCompare averages
   them first and scores that average waveform; for band power, the power of each chosen
   electrode is computed and then averaged. The simple mode scores one measure per run; to compare pipelines on several
   components at once, each with its own electrodes, use the advanced panel (*Add component…*).

   Then tick the steps you want. They are listed in the order they run, and that order is fixed,
   so you only decide which steps are done, never in which order:

   | Order | Step | What is compared |
   |---|---|---|
   | 1 | Bad channels: detect and interpolate | nothing (one fixed setting, see below) |
   | 2 | Reference: as recorded or average | nothing (your choice, the same in every pipeline) |
   | 3 | ICA: remove artifact components | the ICLabel probability above which a component is removed, chosen from the data in each pipeline (every threshold from 0.5 is tried; 0.7, 0.8, 0.9 compared with 7b) |
   | 4 | High-pass filter | the cutoff (0.1, 0.3, 0.5, 1 Hz) |
   | 5 | Low-pass filter | the cutoff (20, 30, 40 Hz) |
   | 6 | Epochs and baseline (segments for band power) | always done |
   | 7 | Reject epochs over an amplitude limit | the limit, chosen from the data in each pipeline (every limit is tried) |
   | 7a | Measure the limit peak-to-peak in 200 ms windows (optional) | the same, on peak-to-peak values |
   | 7b | Instead, repair epochs with up to 3 channels over the limit (optional) | the limit: 75, 100, 150 µV (100, 150, 200 µV peak-to-peak) |

   Two buttons tick a usual set at once: **Standard**, preselected, ticks steps 1 to 7 (every
   step except the options 7a and 7b); **Filters only** ticks the two filters. Lines 7a and 7b
   are options of the rejection and can only be ticked together with step 7 (see
   [Peak-to-peak rejection](#peak-to-peak-rejection) and
   [Repairing epochs with a few bad channels](#repairing-epochs-with-a-few-bad-channels)).
   Ticking 7a does not change the number of pipelines; ticking 7b multiplies it by nine (three
   fixed limits and, with ICA, three ICLabel thresholds instead of both chosen from the data). An unticked step is not done at all: for example,
   without the high-pass the data keep whatever high-pass they already had. Steps that these data
   cannot take are greyed out, with the reason next to them: on epoched data the filters (they
   must run before epoching), without channel locations bad-channel interpolation and ICA (ICLabel
   needs the locations), without the ICLabel plugin ICA, and ICA when the history shows that ICA
   components were already removed. When the average reference is chosen (or the data are
   already average-referenced), bad-channel detection is always ticked, because a bad channel in
   the average would spread into every channel. For band power, when ICA or epoch rejection is
   compared, each filter uses one cutoff, the one nearest the band outside it (outside the band a
   filter does not change its power), so the standard set gives 1 pipeline instead of 12; tick
   only the filters to compare their cutoffs.

   Bad channels are detected once and ICA is fitted once, before the filters, so every filter
   setting shares one decomposition. A channel is bad when its kurtosis (spiky), its joint
   probability (unusual values) or its noise the other channels do not explain (e.g. poor
   contact) is more than 5 SD from the other channels', or when it is flat (no signal); a channel
   bad for only part of the recording is interpolated in that part only. Both steps
   look at a 1 Hz high-passed copy (bad channels: also low-passed at 40 Hz), so slow drifts and
   line noise do not mislead them; the data themselves are
   filtered only by the settings being compared. EOG, ECG
   and EMG channels (by type or by name, such as VEOG) are not tested as bad channels and do not
   count in epoch rejection. The dialog shows how many pipelines will run before you
   start. Steps the data cannot support are left out, with the reason; for example, interpolation
   and ICLabel require channel locations, and for band power no filter edge inside the band is
   compared. Epoched data keep their own epochs when these hold the measurement window. N2pc and
   LRP are scored as ERP CORE measures them, contralateral minus ipsilateral: choose the event
   types of each side (target on the left and on the right; left-hand and right-hand responses).
   To compare ASR (artifact subspace reconstruction) or other steps, use the advanced panel or a
   script.

   Start from the raw continuous data. If your analysis uses the average reference, choose it in the
   **Reference** line of the step list rather than re-referencing beforehand: it is then applied in every pipeline after
   the bad channels are interpolated and before ICA. (Data already average-referenced are averaged
   again after the interpolation.) Steps already applied to the data are read from `EEG.history`
   and listed in the dialog; they are not compared again. A filter the data already have is kept
   as one of the choices (e.g. *low-pass as in the data (30 Hz)*), and the result warns when such
   a filter alone already distorts the measured signal. With the average reference, EEG channels
   removed before PipeCompare are interpolated back first, so the average covers the whole
   montage.
4. Press **Run**. A progress window shows how many pipelines are done, the step running now (for
   example *fitting ICA*), the time so far and the time left. Before the first pipeline is done
   (it runs the steps all pipelines share, such as ICA, which can take several minutes) the bar
   moves back and forth instead of filling. **Stop** ends the search and keeps the pipelines
   already finished. The Command Window prints each step and each finished pipeline as the
   search runs, then a summary; the same log is kept in `pipecompare_last_run.log` in MATLAB's
   `tempdir`.
5. The result window says which pipeline to use and why, naming each pipeline by its settings
   (for example *high-pass 0.5 Hz, low-pass 30 Hz*). **Use this pipeline** stores it as a new
   EEGLAB dataset; **Save script…** writes it as a MATLAB function that you can run on your other
   recordings; **Show all pipelines** lists every pipeline with the reason it was excluded.
   **Use this pipeline** runs the steps again, except ICA, whose decomposition comes from the
   comparison; save the new dataset afterwards with File > Save current dataset as. The new
   dataset is cut into epochs (segments for band power) and cleaned, ready to average into ERPs
   or to compute band power; do not filter, re-reference or reject epochs again.

### Reading the result

The recommended pipeline is marked with `*`, followed by the best of the others. The columns are:

| Column | Meaning |
|---|---|
| checks | *passed*: the pipeline met every constraint; *excluded*: it broke one (the reason is under *Show all pipelines*); *error*: a step failed; *not run*: the search was stopped first |
| error (noise + bias) | how far your averaged measure is from the truth on average: the noise (the standardized measurement error, SME, divided by the share of a known signal the pipeline keeps) and the bias (what blinks add to the mean, and what the filters would change for another signal shape) together, √(noise² + bias²). Lower is better |
| trials kept | the smallest share of trials kept in any condition |
| signal change | how much the pipeline changed the amplitude of the known signal |
| settings | the settings that differ between the pipelines |

The line under the table names the steps that every pipeline shares. When several pipelines cannot
be told apart from the best, PipeCompare recommends the one with the lowest error when each
setting's effect is estimated from all the pipelines that share it (see *Recommendation* under
[How a pipeline is chosen](#how-a-pipeline-is-chosen)); when some are shown to be equally good (within
5 %), the one of those that keeps the most trials, so that data are not cleaned more than they
need to be. The result then says whether the recommended
pipeline's noise is within 5 % of the best one's (a negligible difference) or whether the data are
too few to tell, with how much noisier it could be. It also reports a conditional split check:
some automatic limits, thresholds and component choices are reselected on one half and scored
on the other. Preprocessing, candidate eligibility and shrinkage were fitted using all trials.
The check therefore diagnoses part of the selection advantage; it does not establish that the
complete pipeline generalizes to independent data (see [validation limits](docs/VALIDATION.md)). Last, it
says how often the same pipeline would be recommended with the trials drawn again at random, for
example *this pipeline is recommended 85% of the time: a stable choice* (70 % or more), or *45% of
the time, No. 7 30% of the time: the data do not single out one pipeline clearly*; then the
pipelines it names are about equally good on these data and either can be used. When the
pipelines that decide the recommendation have ICA, they are run again with ICA started in 2 other
ways, and the result says *ICA's randomness does not change the choice*, or that another pipeline
would be recommended and *this recommendation depends on ICA's randomness* (see step 8 of
[How a pipeline is chosen](#how-a-pipeline-is-chosen)). When the settings compared make no difference on
these data, the result says so. When no pipeline passes because most lost too many epochs to
rejection, it names the channels most often over the rejection limit (likely bad channels to
remove or interpolate before running again) and, for data that keep their recorded reference,
suggests the average reference.

When the recommended pipeline interpolated an electrode you measure (it was detected as a bad
channel, or deleted before PipeCompare and put back), the result says so, for example *Pipeline 3
rebuilt Pz, which you measure, from the electrodes around it*. A rebuilt electrode is a weighted
average of its neighbours and has less noise than a real electrode, so that pipeline's noise
(SME) looks smaller than it is and it is favoured over pipelines that keep the electrode. The
real electrode's noise cannot be known (it was bad), so this is not corrected; if the electrode is
really bad, measure at other electrodes (**Electrodes…**) and run again.

When pipelines do pass, PipeCompare still looks at why the recommended pipeline rejected its
epochs. If one channel (or two or three) was over the limit in at least half of the rejected
epochs, and in at least 3 of them, the result names it, for example *Most rejected epochs are due
to one channel: T7 was over the limit in 40 of the 45 epochs rejected in pipeline 2*. Such a
channel is probably bad for long stretches without being bad enough over the whole recording for
the bad-channel test. Look at it in **Plot > Channel data (scroll)**; if it is bad, mark it as bad
or interpolate it (**Tools > Interpolate electrodes**) and run again, which keeps more trials. The
note appears after the headline, in the Command Window and in the advanced panel's notes.

The result also says how large a difference this recording can show, for example *With this many
trials and this noise, two conditions must differ in P3 by about 6 uV to be told apart; smaller
differences need more trials.* This comes from the recommended pipeline's own measurement error
(SME, the error of your measure in its own units, not the gain-corrected value in the table). For
two conditions with errors *a* and *b*, the difference between them has the error
sqrt(*a*² + *b*²); a true difference of 2.8 times that (1.96 + 0.84) is found by a two-sided test
at p < .05 in 80% of recordings like this one. With more than two conditions, the pair with the
largest error is used; with one condition, the value is compared with 0. When the largest
difference actually seen between the conditions (or the value itself, with one condition) is
smaller than this, the result adds that more trials (more events, or several recordings) would
help more than other preprocessing. That comparison is made in the measure's own units, so it
works the same for amplitudes (uV), latencies (ms) and band power (log10 uV²). The same sentence
appears in the Command Window and in the advanced panel's result notes. It describes this one
recording, trial by trial; a group study's power depends on the number of participants instead.

Under the table, the window lists what the recommended pipeline did, step by step in the order the
steps ran; select another row to see that pipeline's steps instead. Each line gives the step's
settings and what it decided on these data, for example:

```
1. Bad channels (kurtosis or joint probability or noise the other channels do not explain over 5 SD, found on a 1-40 Hz band-passed copy): O1, O2 interpolated (not tested: VEOG, HEOG)
2. Average reference (left out: VEOG, HEOG)
3. ICA (extended runica), fitted on a 1 Hz high-passed copy and applied to the data (the fit leaves out stretches far noisier than the rest)
4. High-pass filter 0.5 Hz
5. Low-pass filter 40 Hz
6. ICLabel: 4 of 60 components removed (Muscle, Eye, Heart, Line Noise, Channel Noise with probability 0.8 or more); 21 look like brain activity, 9 were labelled Other
7. Epochs -200 to 800 ms around event type(s) ...
8. Baseline -200 to 0 ms removed
9. Epochs beyond +/-150 uV on any channel rejected: 12 of 120
Trials kept per condition: ...
```

When epochs are repaired (line 8 of the step list), the rejection line adds how many epochs were
kept that way, for example *Epochs beyond +/-150 uV on any channel rejected: 9 of 120; 7 epoch(s)
kept by interpolating up to 3 channel(s) within the epoch*, and with the average reference a last
line *Average reference again (after the epochs repaired by interpolation)* follows.

When ICA did nothing useful, the result says so after the headline (see
[When ICA does nothing](#when-ica-does-nothing)).

This is a readable summary of the pipeline's `EEG.history`; the dataset you get with *Use this
pipeline* keeps the full history (every EEGLAB command) unchanged. The same list is printed in the
Command Window after a run from `pop_pipecompare`, and the advanced panel shows it under *Steps*
when you select a result row.

For any EEGLAB step, order search, several components or different constraints, open
**Advanced…** in the dialog, **Tools > PipeCompare > Advanced panel…**, or the **Advanced** tab of
the PipeCompare window. The panel defines the
same measures as the dialog (ERP components, N2pc and LRP contralateral minus ipsilateral, band
power). See [docs/PANEL.md](docs/PANEL.md) for a description of every control.

### Ear and mastoid channels

Electrodes on the earlobes or mastoids are reference sites, not scalp electrodes. PipeCompare
recognises them automatically, in any dataset, by their standard names: **A1** and **A2** (10-20
earlobes) and **M1** and **M2** (mastoids), in any upper or lower case and also with the `POL `
prefix some EDF exports add (`POL A1`). A channel whose type is set to `REF` in the channel
locations counts too. These channels:

- stay in the data, unchanged;
- are not tested as bad channels, so they are never interpolated (the scalp cannot predict them);
- do not count in epoch rejection;
- are left out of an average reference, as EOG channels are;
- are not interpolated back when they were removed before PipeCompare.

This also covers data already referenced to linked ears or mastoids, where A1 and A2 are flat or
mirror each other and a bad-channel test would wrongly flag them. The dialog (step 1, *Data*),
the Command Window log and the advanced panel's dataset summary list them, e.g. *Ear/mastoid
channels left out: A1, A2*; in the panel they are filled into each step's *exclude* list, where
you can edit them.

Some channels with these names are not ears, and are left as scalp channels:

- caps numbered by letter, such as BioSemi's A1-A32: when the data also have A3 (or M3), A1 and A2
  are scalp electrodes there;
- TP9 and TP10, which sit near the mastoids but are scalp electrodes in the 10-10 system;
- a channel you set to type `EEG` in **Edit > Channel locations** while other channels have other
  types. Use this to keep A1 and A2 as scalp channels. (When every channel has type `EEG`, that is
  the importer's default, not a choice, and the names decide.)

### Montage check

Before the run, PipeCompare checks whether each channel's signal looks like that of the channels
next to it. On the scalp, neighbouring electrodes record similar signals, so a channel that does
not resemble its neighbours is either a bad channel or not where its label says it is (for
example, two labels swapped when the recording was set up or exported). A wrong label matters
more than it seems: bad-channel interpolation and ICLabel both use the positions that come from
the labels.

How it works:

- It uses the EEG channels that have a location, leaving out EOG, ECG and similar channels and the
  ear and mastoid channels. It needs at least 8 of them; with fewer, it is skipped.
- It takes the first 10 minutes of a continuous recording (for epoched data, as many epochs as fit
  in 10 minutes, each with its mean removed), keeps 1-30 Hz and applies an average reference, all
  on a copy. Your data are not changed.
- For each channel it computes the mean correlation with its 3 nearest channels (by electrode
  position).
- A channel is flagged when this correlation is far below the other channels': more than 3 robust
  standard deviations below their median (median and MAD; the spread counts as at least 0.05), or
  below 0.1 while the median of the montage is at least 0.3. On caps with few, widely spaced
  electrodes the median is lower, and then only the first rule applies. The worst channel is
  flagged first and the check runs again without it, so a swapped channel does not also drag its
  neighbours down. At most a quarter of the channels are flagged.
- For each flagged channel it names the two channels it resembles most. If these are far away on
  the head, the label is probably wrong; if it resembles no channel (correlation below 0.3), it is
  more likely a bad channel.

The result appears in the dialog's first line (*Data*), in the Command Window log, and in the
advanced panel's list of data warnings, for example: *Montage check: these channels do not
resemble their neighbours ... AF3 (r = -0.50 with its neighbours; most like TP7 r = 0.61, CP6
r = 0.55) ...*. It is a warning only: nothing is changed and the comparison runs as usual. If
channels are flagged, check the montage with whoever recorded the data, correct the labels or
locations in **Edit > Channel locations**, and run again.

### When ICA does nothing

After the run, PipeCompare looks at what ICA and ICLabel did in the recommended pipeline (or, if
it has no ICA, the first pipeline with ICA), and says so when:

- ICLabel recognised almost none of the components: fewer than 2 components have a Brain
  probability of 50 % or more, or the median probability of *Other* is above 0.8. The message gives
  the numbers, for example *ICLabel recognised almost none of the 61 ICA components: 0 look like
  brain activity (Brain 50% or more) and 56 were labelled Other*. ICLabel judges components partly
  by their scalp maps, so this usually means that the channel labels do not match the electrode
  positions (see [Montage check](#montage-check)); too little clean recording for ICA can do the
  same. In that case the ICA step changed nothing useful.
- No pipeline removed any component. If ICLabel did recognise brain components, the message says
  that no component passed the artifact threshold, which is expected on clean data. When the
  threshold was chosen from the data and ICLabel did take some components for artifacts
  (probability over 0.5), it says instead that removing them did not make the measure measurably
  more precise, so they were kept.
- ICA had too little data. A usual rule is that ICA needs at least 20 × (number of components)²
  data points for a reliable decomposition, more being better: with 60 components, 72 000 points,
  that is 4.8 minutes at 250 Hz or 72 seconds at 1000 Hz. With fewer, the components can mix brain
  activity and artifacts, and ICLabel cannot sort them. The message gives the numbers, for example
  *ICA had too little data: 30000 data points for 61 components, while a reliable decomposition
  needs about 20 x 61 x 61 = 74420 or more*. This is said whatever ICLabel recognised. The fix is
  a longer recording; fewer channels (fewer components) also lower the need.

The note appears in the result window after the headline, in the Command Window, and in the
advanced panel's notes below the results. Each pipeline's step list also gives, on its ICLabel
line, how many components look like brain activity and how many were labelled *Other*.

### Eye channels regressed out

When ICA cannot separate the blinks (see above), they stay in the data. So after the ranking the
recommended pipeline is run once more with one step added: the eye channels regressed out of the
EEG channels (EOG regression). It goes just before the pipeline's first step that removes ICA
components, makes epochs, or rejects or weights epochs; the steps before it make the same
decisions as in the search (the same ICA, the same interpolated channels), the steps after it
decide again on the new data, and the signal check replays the same regression on its copy. The
new pipeline joins the ranking, which is made again: it is recommended only when its total error
(noise and what blinks add, see [How a pipeline is chosen](#how-a-pipeline-is-chosen)) is lower.

How the regression works:

1. The eye channels are the recording's own EOG channels: type EOG (or HEOG, VEOG), or a channel
   that is not scalp EEG whose label contains EOG or EYE (or is LO1/LO2, IO1/IO2, SO1/SO2). Fp1 and Fp2 are not used: they are EEG channels and would
   remove brain activity with them. Without an eye channel of its own, the step is not tried, and
   the Command Window says so.
2. How much of each eye channel reaches each EEG channel is fitted by least squares, only on the
   samples where an eye channel is beyond 4 robust standard deviations of its median, or where
   the hidden Markov model of the eye channels finds a blink (and 100 ms on each side), where
   blinks and eye movements dominate. With less than 1 s of such samples, all
   samples are used.
3. Before fitting, the average of the epochs around each condition's events is taken out of every
   channel. The ERP reaches the eye channels too; without this the fit would remove part of it.
4. The fitted part is subtracted from every EEG channel; the eye channels are left as they were.

Its step line says *Eye channels EOG1, EOG2 regressed out of the EEG channels (fitted on 40 s of
eye activity)*, and its name ends with *eye channels regressed out*. In a saved script it is one
line, `EEG = pipecompare.run.Steps.eogRegress(EEG, {...}, {...}, [...]);`. It is on in the simple
mode; in the advanced panel tick *Options… > Try regressing the eyes out*; in scripts pass
`'eogCheck', true` to `pipecompare.PipeCompare.optimize`. Its filter cutoffs are not fine-tuned
afterwards. The step can also be added to a plan as `eogregress`.

### Blinks and bad stretches found by a hidden Markov model

Two things PipeCompare needs to know are a matter of *when*: when the eyes blink or move (to
estimate what blinks add to each pipeline's result, to fit the EOG regression, and to leave the
blinks out when looking for bad channels), and when a channel is bad (an electrode can come
loose halfway through). A fixed rule (beyond so many standard deviations, plus a margin) misses
what is weak in any one sample but clear over a stretch. So both are found by a **two-state
hidden Markov model**: the data are taken to switch between a quiet and an active state, how
large each state's values are and how long each state usually lasts are learned from the
recording itself (the EM, or Baum-Welch, algorithm), and the most probable sequence of states is
found exactly (the Viterbi algorithm, by dynamic programming).

**Blinks and eye movements.** Each eye channel (and the left minus right pair) minus its moving
median over 1 s, smoothed over 1/16 s and taken at about 25 Hz; the largest over the channels of
how far its size is from usual (robust z). A stretch in the active state counts as a blink when it
reaches 3.5 robust SDs and lasts 100 ms or more; without that, the model would split the noise of
a recording without blinks into two states. The samples found this way are added to those the
fixed rule finds (beyond 5 robust SDs), so nothing the rule found is lost. In simulations (eye
channel noise 10 µV, blinks of 15 to 200 µV), blinks of 15–30 µV were found 64 % of the time
instead of 11 % (with slower noise, 28 % instead of 3 %), those of 30–50 µV all instead of
83 %, at about 2 extra short stretches a minute in a recording without blinks.

**Channels bad for part of the recording.** For each channel, its noise the other channels do
not explain (as in the bad-channel test) is measured in 1 s windows, and what all channels share
in a window (a movement) is taken out. On that sequence the model finds the stretches in which
the channel is noisier than usual. A stretch counts when it lasts 10 s or more and its noise is
on average 3 robust SDs above the channel's usual level (measured on the channel's quieter part,
so a channel bad for most of the recording is still measured against its good part). The
channel is then interpolated (spherical, from the channels not bad there) in that stretch only;
the rest of it is left as it was. The step line says, for example, *bad for part of the
recording, interpolated there only: P8 (1210-1500 s)*, and a saved script repeats it with
`pipecompare.run.Steps.interpStretches`. In a simulation (32 channels, 20 minutes, Pz from some
time on 1.5 to 3 times noisier), the whole-recording test found such a channel only when it was
bad for half the recording at 3 times the noise; the model found every bad stretch, and the
noise at Pz came to 1.02–1.10 times that of a clean recording, instead of 1.06–1.64 with the whole-recording test alone (interpolation assumed 1.2 times as noisy as a clean channel). It
needs continuous data, 8 channels or more and 30 s or more, and channels that are not linearly
dependent (it runs before the average reference).

### Peak-to-peak rejection

The usual rejection drops an epoch when any EEG channel goes beyond ±*x* µV. That test looks at
the distance from zero, which after baseline correction is the distance from the baseline. Slow
drifts can carry a clean epoch over the limit, while a blink that starts from a low point can stay
under it. ERP CORE (Kappenman et al., 2021) therefore uses a moving-window peak-to-peak test, as
ERPLAB's *moving window peak-to-peak threshold*:

- a window of 200 ms moves over the epoch in steps of 50 ms;
- in each window, for each EEG channel, the highest value minus the lowest value is taken;
- the epoch is rejected when, in any window, this exceeds the limit on any EEG channel (EOG and
  other non-EEG channels, and ear or mastoid channels, are not tested, as in the usual test).

A blink or a fast jump is a large change within 200 ms and is caught; a slow drift changes little
within 200 ms and is not. Because the limit is a range rather than a distance from zero, its
values are larger: simple mode chooses the peak-to-peak limit from the data, as for the usual
test, and compares 100, 150 and 200 µV peak-to-peak when epochs are repaired (tick **measure the
limit peak-to-peak in 200 ms windows** under the rejection; `'peaktopeak'` in `pop_pipecompare`).
It can be combined with repairing epochs: the channels over the limit are the ones repaired.

In the advanced panel, the amplitude rejection step (`reject_threshold`) has the parameters
`method` (`absolute`, the default, or `peaktopeak`; both can be compared, as `absolute |
peaktopeak`), `window` (the window length in ms, default 200) and `uv` (any limits, or `auto`). Each
pipeline's `EEG.history` records the test as `pipecompare.run.Steps.markPeakToPeak(...)`, which
marks the epochs where EEGLAB's own threshold test puts its marks (`EEG.reject.rejthresh`).

### Rejection limit chosen from the data

Which rejection limit is best depends on the recording: a limit that is too strict throws away
good trials, one that is too loose keeps artifacts, and both make the measure less precise.
Instead of trying three limits, PipeCompare finds the best one in each pipeline, from that
pipeline's data at the rejection step:

1. For every epoch it takes the value the limit is compared with: the largest absolute value over
   the tested channels and the whole epoch (peak-to-peak: the largest peak-to-peak value in the
   moving windows). An epoch is rejected when this value is above the limit.
2. A limit between two neighbouring values of this list keeps exactly the same epochs, so there
   are only as many different results as epochs, and all of them are tried: the limit is the exact
   best, not the best of a few tries. For each one the SME the ranking uses (RMS over the
   conditions, and over the measures for several components) is computed at once from sums over
   the kept trials.
3. Only limits that keep at least 10 trials and 50 % of the trials of every condition (the
   ranking's own limits, `minTrials` and `minRetention`) count.
4. Limits whose SME is shown to be within 5 % of the best limit's (the *equivalent* test of the
   ranking: the same paired bootstrap and simultaneous intervals, over up to 60 limits spread over
   the range, the best included) are equally good, and the **highest** of them is used: an epoch
   is kept unless keeping it makes the measure measurably less precise. A limit the data only
   cannot tell apart from the best is not enough: with few trials that includes limits that
   reject almost nothing even when rejecting clearly helps (the SME then gets larger). The number used is the roundest
   one among the limits that keep the same epochs (for example 120 µV rather than 117.43 µV).
5. The limit is then applied with EEGLAB's own test (`pop_eegthresh`, or the peak-to-peak marks),
   so `EEG.history` records the number, and the signal check removes the same epochs from its
   copy.

The Command Window shows the limit used and the best one, and the steps of each pipeline in the
results say *limit chosen from these data* with the value. It needs a score per trial, so it works
for mean amplitude and band power; with peak measures, or when epochs are repaired (the limit
then also decides which channels are interpolated), fixed limits are compared. The scores are
those of the data at the rejection step; it is exact when nothing after the rejection changes them
(the baseline is already removed from them). In the advanced panel and in scripts, set the
`reject_threshold` parameter `uv` to `auto` (it can be compared with fixed limits, as
`auto | 100`). A saved script keeps `auto`, so on another recording the limit is found from that
recording's data.

**Rejected for its P3?** A rejection looks at the largest value anywhere in the epoch, the measured
electrodes and window included. A low limit can therefore reject a trial because its P3 is large,
and the average of the trials kept is then smaller than the real one. The signal check cannot see
this, because the known signal is the same in every epoch. So every rejection (with a fixed limit
too) compares, at the measured electrodes in the window, the average of the trials kept with that
of all trials, as a z value: their difference divided by the spread it would have if the rejected
trials were a random draw. When the trials kept have clearly less of the component (z below −3),
the result says so for the recommended pipeline, with both averages. If the rejected trials are
artifacts that happen to push the window the same way (blinks at the P3's electrodes, for example),
nothing needs to change; otherwise a higher limit, or leaving the measured electrodes out of the
test, avoids it.

### ICLabel threshold chosen from the data

Which ICLabel threshold is best also depends on the recording: a low threshold removes components
that hold brain signal, a high one leaves eye and muscle activity in, and the epochs it leaves in
change which rejection limit is best. Instead of trying three thresholds, PipeCompare finds the
best one in each pipeline, together with the rejection limit that follows:

1. For every component it takes the largest ICLabel probability among the artifact classes (eye,
   muscle, heart, line noise, channel noise). EEGLAB's `pop_icflag` removes a component when this
   value is above the threshold (and below 1), so a threshold between two neighbouring values
   removes exactly the same components. There are only as many different results as components
   with a value over 0.5 (more likely that artifact than anything else), plus removing nothing,
   and all of them are tried: the threshold is the exact best, not the best of a few tries.
2. Removing components is a linear operation (the data minus the removed components' part, as
   `pop_subcomp` computes it), and so are epoching and baseline removal, which come between it
   and the rejection. So each set of components is taken out of the epoched data directly,
   without running the steps again.
3. For each set, the rejection limit is then found as in
   [Rejection limit chosen from the data](#rejection-limit-chosen-from-the-data) (every limit;
   a fixed limit when the step has one), and the SME at the best limit is that set's score. The
   SME is divided by the signal gain, which removing components can change: the known signal of
   the signal check goes through the same removal, and a set that fails the signal check (for
   example because it removes a component that carries the measured signal) does not count, nor
   does one that keeps too few trials.
4. Sets whose SME is shown to be within 5 % of the best set's (the *equivalent* test of the
   ranking: the same paired bootstrap and simultaneous intervals) are equally good, and the one
   removing the **fewest** components is used: a component is kept unless keeping it makes the
   measure measurably less precise. The number used is the roundest threshold that removes exactly that set
   (for example 0.8 rather than 0.8137).
5. The threshold is then applied with EEGLAB's own functions (`pop_icflag`, `pop_subcomp`), so
   `EEG.history` records the number, and the rejection step that follows finds its limit on the
   data that are left, as it would after any threshold.

The Command Window shows the threshold used and the best one, and the steps of each pipeline in
the results say *threshold chosen from these data* with the value. It needs a score per trial, so
it works for mean amplitude and band power; with peak measures, or when epochs are repaired, fixed
thresholds are compared. It is exact when only epoching, baseline removal and an amplitude
rejection that removes epochs come after it, which the plan checks; the signal check of each set
uses all epochs (the rejected ones change it only where epochs overlap), and the full pipeline's
own signal check is the one the ranking uses. In the advanced panel and in scripts, set the
`icremove` parameter `threshold` to `auto`; every pipeline below it must then have the same steps
after it. When the search adopts a pipeline, the threshold it chose is reused; a saved script keeps
`auto`, so on another recording the threshold is found from that recording's data.

### ICA components chosen by their effect on the noise

ICLabel decides from what a component looks like (its scalp map, spectrum and time course). On
some recordings it recognises little (for example when most components are labelled *Other*),
and then no threshold removes the eye or muscle activity that is still in the data. With this
option ticked (*also choose the components by their effect on the noise*), every pipeline is
also run with the components chosen by what removing them does to the measure:

1. The components that may be removed are all those whose likeliest ICLabel class is not
   *Brain*: a component ICLabel takes for brain activity is never removed. A component that
   carries clear time-locked activity into the measured window is not removed either, even when
   ICLabel calls it *Other*: removing it alone would move a condition's average in the window by
   more than half its SME, with a t value above 3 across trials. The noise measure counts real
   brain activity that varies from trial to trial as noise, so without this a P3 component that
   ICLabel did not recognise could be removed for lowering it. Eye components are not kept this
   way, since blinks can be time-locked to the stimulus too: those ICLabel calls *Eye*, and those
   whose activity follows an eye channel far more closely than the other components' (ICLabel
   often calls a blink component *Other*). The components removed together may not move the
   average this way either, so a P3 spread over several weak components is not removed piece by
   piece. The steps of the pipeline say how many were kept for these reasons.
2. The components are put in order by what each one does to the measure on its own: how much it
   makes the score vary from trial to trial, divided by how much of the known signal it carries
   there. The sets are the first one, the first two, and so on. (Adding each time the component
   that lowers the SME most, as before, fits the chance of a few dozen trials: in simulations at
   35 trials per condition its SME looked 38 % better than it was.) Each set is scored
   exactly as in [ICLabel threshold chosen from the data](#iclabel-threshold-chosen-from-the-data):
   the components are taken out of the epoched data directly (removal is linear), the rejection
   limit that follows is found at its best, the SME is divided by the signal gain, and a set that
   fails the signal check (for example because it removes the measured signal) or keeps too few
   trials does not count. This gives a path of sets, each one component larger, ordered by what
   the components do to the measure.
3. The path stops after 30 components, or when the SME has risen more than 10 % above the best
   on the path.
4. Among the sets on the path, those whose SME is shown to be within 5 % of the best one's are
   equally good, and the one removing the **fewest** components is used, as for the threshold.
5. The components are flagged (`EEG.reject.gcompreject`, recorded in `EEG.history`) and removed
   with EEGLAB's `pop_subcomp`.

Both ways, the ICLabel threshold and this one, are compared in the ranking. The steps of each
pipeline in the results say *ICA components chosen by their effect on the noise* with how many
were removed. It needs a score per trial (mean amplitude or band power) and epochs that are
removed, not repaired; it takes longer than the threshold (one score per component and step).
It chooses what lowers the noise of this measure in this recording, so a component that only
matters for another measure may stay. In the advanced panel and in scripts, set the `icremove`
parameter `threshold` to `sme` (or compare `auto | sme`); as with `auto`, every pipeline below it
must have the same steps after it. When the search adopts a pipeline, the components it chose are
reused; a saved script keeps `sme`, so on another recording they are chosen from its own data.

### Filter cutoffs fine-tuned from the data

Comparing 4 high-pass and 3 low-pass cutoffs leaves the best value between two of them untried
(for example 0.2 Hz between 0.1 and 0.3 Hz). Filtering the whole recording again for many more
cutoffs would take long, so PipeCompare uses a shortcut that gives the same numbers:

1. It takes the recommended pipeline and keeps what it decided on these data fixed: the
   components it removed and the epochs it rejected.
2. With those fixed, everything from the filters to the scores is adding and multiplying the data
   by fixed numbers. Filtering acts on each channel alone, and removing components is the same
   mixture of channels at every moment, so the order of the two does not matter. The scores
   therefore need only the electrodes you measure: the components are removed, the other channels
   are left out, and only the measured electrodes are filtered, epoched and scored, which gives
   exactly the numbers the whole data would give there.
3. Every pair of cutoffs on a finer grid between the lowest and highest compared ones is scored
   this way (with the default steps 0.1, 0.15, … 1 Hz and 20, 22, … 40 Hz: 209 pairs). The signal
   gain comes from the pipeline's own signal check, corrected by what the filters alone do to the
   known signal, and cutoffs whose filters alone distort the signal are left out. At the
   pipeline's own cutoffs the result must equal the pipeline's SME, as a check; otherwise nothing
   is fine-tuned.
4. When another pair is better, the pipeline with those cutoffs is run like every other one (its
   components and epochs are decided again on its own data, with the full signal check) and added
   to the results; its filters say *fine-tuned from these data*. The ranking then decides whether
   it is better.

The Command Window shows the pair found and the SME before and after. It is done for mean
amplitudes of an ERP on continuous data when only filters, IC removal, epoching, baseline removal
and an amplitude rejection (without repairing epochs) follow the first filter; it is on in the
simple mode. In the advanced panel tick *Options… > Fine-tune filter cutoffs*; in scripts pass
`'refineCutoffs', true` to `pipecompare.PipeCompare.optimize`.

### Epochs weighted by their noise

Rejection treats epochs as all or nothing: an epoch just over the limit is dropped, one just under
it counts as much as the cleanest epoch. Weighting instead counts each epoch by how clean it is,
and drops none:

1. For each epoch, its **noise** is measured at the electrodes you measure: the epoch minus the
   average of the *other* epochs of its condition (so an epoch cannot make itself look clean),
   squared and averaged over those electrodes and all samples of the epoch, leaving out the
   part of it that goes with the measured score (otherwise an epoch would count less for having a
   score far from the others', and the weighted average would look more precise than it is).
2. Its **weight** is 1 / noise, scaled so that the weights of a condition average 1. If epochs
   differ only in how much noise they have, these weights give the average with the smallest
   possible error (inverse-variance weighting). Plain averaging and rejection are both special
   cases of weighting (weights 1, or 1 and 0), so ideal weights are never worse than either;
   weights estimated from the data can be, which is why the ranking compares them. So that no epoch can
   count too much, a noise below a tenth of the condition's median counts as that tenth: no epoch
   counts more than 10 times the typical one.
3. The **average** of a condition is then sum(weight × epoch) / sum(weight), and its error (SME)
   is computed for a weighted mean (the robust "sandwich" formula: sqrt(n / (n − 1) ×
   sum(weight² × (score − average)²)) / sum(weight), which gives SD / sqrt(n) with equal weights
   and stays correct even when the weights are not ideal). The bootstrap of the ranking and the
   selection check resample the epochs with their weights.

In *Standard* this comes after the rejection: every pipeline is run once with equal weights and
once weighted, and the ranking decides. The rejection limit is chosen for equal weights. The
known signal is the same in every epoch, so weighting does not change the signal check. Weights
are for ERPs only: for band power, a weighted mean would favour the epochs with less power.

The weights are stored in the dataset (`EEG.etc.pipecompare.trialWeights`, by the event each
epoch is locked to, so they stay with their epochs) and listed in the steps as *Epochs weighted by
their noise*, with how much more the cleanest epoch counts than the noisiest. EEGLAB's and
ERPLAB's own averaging do not use them, so averaging such a dataset there would not give the ERP
the pipeline was scored on. Therefore, when the recommended pipeline is weighted, the result says
so and names the same pipeline with equal weights (use that one to average in EEGLAB or ERPLAB as
usual). When you use a weighted pipeline (*Use this pipeline*, `adopt`) or run its saved script,
the weighted average of each condition is stored in the dataset:

```matlab
A = EEG.etc.pipecompare.weightedERP;   % A.data: channels x time points x conditions
A.conditions, A.times                  % condition names; times in ms
plot(A.times, A.data(strcmpi(A.labels, 'Pz'), :, 1))   % e.g. the first condition at Pz
```

`pipecompare.eval.Measure.weightedAverages(EEG, contract)` computes it again. To average with
the weights yourself:

```matlab
w = pipecompare.eval.Measure.epochWeights(EEG);              % one weight per epoch
erp = sum(EEG.data .* reshape(w, 1, 1, []), 3) / sum(w);     % channels x time points
```

(for one condition, take its epochs and their weights only). The PipeCompare window's plots of a
pipeline use them. In the advanced panel the step is *weight_trials*; in a plan,
`p.addChoice('weight', {'weight_trials'}, 'none')` compares it with equal weights.

### Signal check with your recording's own waveform

The signal check carries a known signal through every pipeline (see
[How a pipeline is chosen](#how-a-pipeline-is-chosen)). By default its shape is a smooth bump in
each measurement window. Your real P3 may be broader, narrower or asymmetric, and a filter treats
those shapes differently. In the simple mode (and with *Options… > Signal check with own waveform*
in the advanced panel, `'signalTemplate', 'data'` in scripts) the shape is taken from your own
recording, before any pipeline runs:

1. The epochs of each condition are averaged at the measured electrodes, the baseline removed, and
   the conditions averaged with equal weight. The average is smoothed lightly (a 5 ms Gaussian,
   which keeps everything below about 25 Hz).
2. Only the component's own stretch is kept: the measurement window, faded to zero over half the
   window's width (at least 25 ms) on each side. What lies outside belongs to other components.
3. It is used only when it is **reliable**: the averages of the odd and of the even epochs must
   look alike over that stretch (split-half reliability, stepped up to all epochs with the
   Spearman-Brown formula, at least 0.8), and it must be one deflection of one sign in the window.
   A noisy average would make the check test noise instead of your component.
   It must also not look like a **blink**: in many tasks people blink right after a target, at
   about the P3's latency, and such a blink is as reliable as a component. So the eye channels
   (EOG by type or name, for example EYEL/EYER; without them Fp1, Fp2 and FPz) are averaged the
   same way, and when one of them has the same shape (similarity above 0.8) and is more than twice
   as large, the bump is used and the result says it looks like eye activity.
4. It is scaled to a peak of 5 µV in the window and checked exactly as the bump would be.

When the average is not reliable enough (in a simulation with a 6 µV P3 in 8 µV noise: about half
the time with 20 epochs per condition, nearly never with 90), the bump is used for that component,
and the result says which shape was used and why. Contralateral-minus-ipsilateral components (N2pc,
LRP) and band power always use the standard signal. The shape comes from the data before any
pipeline, unfiltered apart from the smoothing, so a slow drift that the averaging did not cancel can
remain in it. A saved script uses the same shape for its own choices, estimated from the data it
is run on.

### ICA on clean data

ICA finds the independent sources in the data. Stretches of the recording that are far noisier
than the rest, such as movement or a loose electrode for a few seconds, take up components of their
own and leave fewer for eye, muscle and brain sources, so the decomposition is worse for the
whole recording. PipeCompare therefore fits ICA on a copy without those stretches, and applies the
result to all of the data (nothing is cut from your data):

- continuous data are divided into 1 s windows (epoched data: each epoch is one window);
- for each window, the spread (standard deviation) of every EEG channel is taken, and the mean of
  their logarithms gives the window's noise level (EOG and other non-EEG channels do not count,
  so blinks stay in and ICA can learn them);
- a window is left out when its noise level is more than 3 robust standard deviations above the
  median of all windows (median and MAD), at most 20 % of the windows, the noisiest first.

On clean data nothing is left out. The Command Window and each pipeline's `EEG.history` say how
much was left out. In the advanced panel, the ICA step's parameter `fitClean` (1 by default)
turns this off with 0, and can be compared (`0 | 1`).

### Repairing epochs with a few bad channels

Sometimes a single electrode loses contact for a moment: in a few epochs one channel goes far over
the rejection limit while all the others are fine. Rejecting those epochs loses trials because of
one channel. With **repair epochs with up to 3 channels over the limit** ticked (line 8 of the
step list; `'epochinterp'` in `pop_pipecompare`), the rejection step works as follows, in each
epoch:

1. It tests every EEG channel against the limit, as usual (non-EEG and ear/mastoid channels are not
   tested).
2. If 1, 2 or 3 channels are over the limit, and they all have locations, those channels are
   replaced in that epoch only by spherical-spline interpolation from the other channels of the
   same epoch (the same computation as EEGLAB's `eeg_interp`), and the epoch is kept. An electrode
   you measure (for example Pz for the P3) is never repaired: if it is one of the channels over the
   limit, the epoch is rejected. A rebuilt electrode is an average of its neighbours and has less
   noise than a real one, so repairing it would make the pipeline look more precise than it is.
3. If more than 3 channels are over the limit, the epoch is rejected as before.

This is the idea of the epoch-level channel interpolation in FASTER (Nolan, Whelan & Reilly,
2010). It runs after ICA, filtering, epoching and baseline correction, as part of the rejection
step, so it sees the same data the limit is applied to. Repaired epochs are not tested again.
With the average reference, the data are averaged again afterwards: the replaced values were part
of that epoch's average, which would otherwise keep a share of the bad signal in every channel.
The signal check applies exactly the same repairs (same channels, same epochs) to its copy, so any
change to the measured signal is counted. Repaired epochs do not count towards the 20 % limit on
interpolated channels, which is about whole channels. Each pipeline's `EEG.history` records the
repairs as one command listing the channels and epochs.

In the advanced panel, every rejection step (amplitude, joint probability, kurtosis) has the
parameter `interpolate`: the largest number of flagged channels an epoch may have and still be
repaired (0, the default, turns it off). Set any number there, or several (for example `0 | 3`)
to compare pipelines with and without repairs. For joint probability and kurtosis, the flagged
channels are those over the per-channel limit.

### From the command line

Everything in the menu can also be started from MATLAB's Command Window. PipeCompare works on the
dataset that is current in EEGLAB, so start EEGLAB and load the dataset first; this creates the
variables `EEG`, `ALLEEG` and `CURRENTSET` that PipeCompare reads and writes:

```matlab
eeglab                                                       % start EEGLAB (adds PipeCompare to the path)
EEG = pop_loadset('filename', 'mydata.set', 'filepath', '/path/to/folder');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);       % make it the current EEGLAB dataset
eeglab redraw
```

When the dataset is already open in EEGLAB's window, skip these lines. If `EEG` in the Command
Window is not the dataset selected in EEGLAB (for example after `EEG = ...` on another dataset),
`pop_pipecompare` stops with *make this dataset current first*; `EEG = ALLEEG(CURRENTSET);` makes
them the same again.

Open the windows of the menu:

```matlab
pcomp                                                   % the PipeCompare window (or pipecompare.start; Tools > PipeCompare > PipeCompare window...)
[EEG, LASTCOM] = pop_pipecompare(EEG); if ~isempty(LASTCOM), eegh(LASTCOM); end   % the simple dialog (Tools > PipeCompare > Compare pipelines...)
pipecompare.PipeCompare.app();                          % the advanced panel (Tools > PipeCompare > Advanced panel...)
```

`eegh(LASTCOM)` adds the command that repeats your choices to EEGLAB's history, as the menu does.

Or skip the dialog and give the choices as options. The menu dialog records an equivalent
command in EEGLAB's command history (`ALLCOM`):

```matlab
% ERP: compare pipelines for the P3, one condition per event type
EEG = pop_pipecompare(EEG, 'measure', 'P3', 'events', {'target', 'standard'}, 'recipe', 'standard');

% N2pc, contralateral minus ipsilateral: the event types of each target side
EEG = pop_pipecompare(EEG, 'measure', 'N2pc', 'left', {'111', '112'}, 'right', {'121', '122'});

% Continuous data: compare pipelines for alpha-band power in 2 s segments
EEG = pop_pipecompare(EEG, 'measure', 'alpha', 'recipe', 'filters');

% Only some steps (any of 'badchannels', 'ica', 'highpass', 'lowpass', 'reject', 'epochinterp'),
% always run in that order
EEG = pop_pipecompare(EEG, 'measure', 'P3', 'events', {'target'}, 'steps', {'highpass', 'lowpass', 'reject'});

% Standard plus repairing epochs that have up to 3 channels over the rejection limit
EEG = pop_pipecompare(EEG, 'measure', 'P3', 'events', {'target'}, ...
    'steps', {'badchannels', 'ica', 'highpass', 'lowpass', 'reject', 'epochinterp'});

% A component with your own electrodes instead of its ERP CORE site
EEG = pop_pipecompare(EEG, 'measure', 'P3', 'channels', {'Pz', 'CPz', 'POz'}, 'events', {'target'});

% Your own window (s) and electrodes; or 'measure', 'band', 'band', [8 12]
EEG = pop_pipecompare(EEG, 'measure', 'custom', 'window', [0.25 0.5], 'channels', {'Cz', 'CPz'}, ...
    'events', {'target'}, 'recipe', 'standard');
```

Replace `'target'` and `'standard'` with the event types in your dataset. With options, the
progress window (with **Stop**) and the result window open as from the menu; add `'show', 'off'`
to run without any window, for example in a script that runs unattended (the Command Window
still prints each step):

```matlab
EEG = pop_pipecompare(EEG, 'measure', 'P3', 'events', {'target'}, 'reference', 'average', 'show', 'off');
```

The dataset is returned unchanged. The result is stored in the variable `pipecompare_result`
(also the third output, `[EEG, com, result] = pop_pipecompare(...)`), and these commands do what
the result window's buttons do:

```matlab
pipecompare.PipeCompare.adopt(pipecompare_result);           % Use this pipeline: new EEGLAB dataset, now current
pipecompare.PipeCompare.adopt(pipecompare_result, 12);       % ... or pipeline 12 (one that passed the checks)
pipecompare.PipeCompare.writeScript(pipecompare_result, [], 'my_pipeline.m');   % Save script...
pipecompare.PipeCompare.script(pipecompare_result);          % print the recommended pipeline's EEGLAB commands
```

After `adopt`, `EEG` in the Command Window is the new dataset (cut into epochs and cleaned); save
it with `pop_saveset`. The saved script is a function: run it on any recording with the same
channels and event types, for example `EEG = my_pipeline(EEG);` (the file must be in the current
folder or on the MATLAB path). The file name is the function's name, so it may contain only
English letters, digits and underscores, must start with a letter, and must not be the name of an
existing MATLAB or EEGLAB function (such as `pop_epoch`); PipeCompare refuses other names. Type
`help pop_pipecompare` for all options.

### One pipeline for a group of recordings

For a study you usually want **one** pipeline for every participant, not a different one each.
The *6. Group* page of the PipeCompare window does this for the recordings you select (see
[In the PipeCompare window](#in-the-pipecompare-window)); from the command line,
`pipecompare.Group.run` takes the recordings (a list of `.set` files, or datasets in EEGLAB) and
the same options as `pop_pipecompare`:

```matlab
files = {'/path/to/sub01.set', '/path/to/sub02.set', '/path/to/sub03.set'};
G = pipecompare.Group.run(files, 'measure', 'P3', 'events', {'target', 'standard'}, 'recipe', 'standard');
G.table                                        % one row per pipeline; G.recommended is the row to use
pipecompare.Group.writeScript(G, 'my_group_pipeline.m');   % then: EEG = my_group_pipeline(EEG); per recording
G.own                                          % each recording's own pipeline
pipecompare.Group.writeScript(G, 'sub02_pipeline.m', 2);  % the own pipeline of the 2nd recording
```

(or `G = pipecompare.Group.run(ALLEEG, ...)` for the datasets open in EEGLAB). It works like this:

1. Each recording is compared on its own, exactly as `pop_pipecompare` does (without windows). The
   steps still decide from each recording's own data: its bad channels, its ICA components, its
   rejected epochs and its automatic limits.
2. A pipeline counts for the group only when it **passes on every recording**. A recording on
   which no pipeline passes is left out and named (`G.leftOut`); look at it on its own.
3. The pipelines are ranked by the **root mean square of individual total errors over the
   recordings**. This summarizes the observed group's measurements; it is not the uncertainty
   of a group mean. Systematic bias shared across participants need not average away.
4. Which pipelines cannot be told apart from the best, and which are within 5 % of it, come from
   each recording's bootstrap of its trials, combined the same way. Among those that cannot be told
   apart, the group rule recommends the one that keeps the most trials in its worst recording.
5. **Each recording's own pipeline** is also given (`G.own`: recording, pipeline, and whether it
   is the group's). Picking on one recording alone is unreliable with a few dozen trials, and one
   pipeline for all ignores real differences between recordings. So each pipeline's error on a
   recording is estimated together with the others: a part common to the whole group (how good
   the pipeline is on average, and how noisy the recording is), plus what is particular to this
   recording, estimated with the full paired bootstrap covariance (two-level empirical Bayes).
   Correlated estimation noise can cancel between pipelines; it is not treated as independent
   evidence against individual differences. Missing candidate patterns use the actual design
   rank. The Command Window shows own pipelines and the average posterior gain (`G.share`);
   `G.ownShrunk` contains the estimates and `G.twoLevelModel` describes the fitted design.
   This fits the observed group; benefit for new participants still needs independent validation.

The Command Window shows each recording's comparison, then the group's recommendation. `G` is also
stored in the variable `pipecompare_group`; `G.perRecording` has each pipeline's noise on each
recording. The current dataset in EEGLAB is left as it was.

## Scripting interface

For full control, describe what is measured (an analysis *contract*) and what may vary (a *plan*):

```matlab
pipecompare.PipeCompare.state();          % summary of the current dataset and its history

% What is measured: conditions, epoch, baseline, and the measure (window, channels, type)
c = pipecompare.eval.Contract( ...
    'conditions', {'target', {'target'}; 'standard', {'standard'}}, ...
    'epoch',      [-0.2 1.0], ...
    'baseline',   [-0.2 0], ...
    'components', {'P3', [0.30 0.60], {'Pz', 'CPz', 'POz'}, 'mean'});

% What may vary
p = pipecompare.plan.Plan();
p = p.add('highpass');                                   % cutoff searched over 0.1, 0.3, 0.5, 1 Hz
p = p.add('lowpass', 'cutoff', 30);                      % fixed value
p = p.addEeglab('EEG = pop_reref(EEG, []);', 'reref');   % one EEG = pop_*(EEG, ...) call as a step
p = p.add('ica', 'fitHighpass', 1);
p = p.add('icremove', 'threshold', {0.8, 0.9});          % searched over two values
p = p.add('epoch');                                      % windows come from the contract
p = p.add('baseline');
p = p.add('reject_threshold', 'uv', {100, 150}, 'interpolate', {0, 3});   % with and without repairing epochs

r = pipecompare.PipeCompare.optimize(p, c, struct('checkpoint', 'pc_run1'));

pipecompare.PipeCompare.adopt(r);                          % recommended pipeline -> new dataset
pipecompare.PipeCompare.writeScript(r, 3, 'pipeline3.m');  % candidate 3 as a function for any recording
r = pipecompare.PipeCompare.resume('pc_run1');             % continue an interrupted search
```

Band power of continuous data uses a band-power contract and the same plan without the `baseline`
step (segments have no baseline window, so a `baseline` step would make every pipeline fail):

```matlab
c = pipecompare.eval.Contract('analysis', 'bandpower', 'segment', 2, ...
    'bands', {'alpha', [8 12], {'O1', 'Oz', 'O2'}});
```

By default, steps run in the order they were added. If that order is invalid, PipeCompare lists
the conflicts instead of rearranging the steps. Set `p.OrderMode = 'search'` to compare every
valid order; `p = p.pin(id)` keeps a step in place and `p = p.before(a, b)` makes step `a` run
before step `b` (a plan is a value object, so assign the result).

The available steps are listed in [docs/COVERAGE.md](docs/COVERAGE.md). Search options (for
example `parallel`, `checkpoint` and `maxLeaves`) are documented in `help pipecompare.run.Executor`,
and constraints and ranking options in `help pipecompare.eval.Rank`.

## How a pipeline is chosen

1. **Feasibility.** Pipelines that fail or violate a constraint are excluded, each with a stated
   reason. Default constraints:
   - at least 10 trials and 50 % retention per condition;
   - at most 20 % of channels interpolated;
   - preservation of the known signal: amplitude error ≤ 10 %, peak shift ≤ 10 ms, artifactual
     deflection ≤ 5 %, waveform correlation ≥ 0.95, topography correlation ≥ 0.90 (for band power
     only the amplitude error and topography correlation apply).

   The filters are checked first. Filters act on a signal the same way whatever the data, so
   what a pipeline's filters alone (high-pass, low-pass, line-noise filter) do to the known
   signal is computed before anything runs, on a short stretch of the signal copy around the
   first events (with more than a filter length on each side). A pipeline whose filters alone
   break one of these limits (the topography aside, which a filter does not change) is excluded
   without being run, with a reason that starts with *filters alone*; for example a 1 Hz
   high-pass for a broad P3. This saves the filtering of the whole recording and everything after
   it for those pipelines. It is done for ERP measures on continuous data, when no resampling or
   epoching comes before the last filter; `'filterCheck', false` turns it off.
2. **Total error.** Feasible pipelines are scored by how far their measure is from the truth on
   average: the noise and the bias together, √(noise² + bias²).
   - The **noise** is the SME divided by the gain each pipeline applies to the known signal. This
     makes the comparison one of signal-to-noise ratio, so a pipeline cannot appear more precise
     just by attenuating everything (Zhang, Garrett & Luck, 2024).
   - The **bias** is what moves the mean the same way in every set of trials, which the SME cannot
     see. It has two parts. First, what blinks add: blinks time-locked to the stimulus (for
     example after targets) fall in a share of the trials and move the mean much more than they
     raise the SME. PipeCompare finds the trials with a blink in the eye channels of the data as
     they are (EOG channels, else Fp1/Fp2/FPz when they are not measured: an epoch whose eye
     activity in the window, relative to the baseline, is more than 3 robust SDs above the
     others'). The same is done for two other artifacts that can be locked to the stimulus: eye
     movements, in the difference of a left and a right eye channel (left and right EOG by
     position or by names ending in L and R, e.g. EYEL/EYER; else F7 and F8), where a blink
     cancels and a horizontal movement doubles; and muscle activity, the fast activity (the SD of
     the sample-to-sample change) in the window at the measured electrodes. A trial with any of
     the three counts. For each pipeline PipeCompare computes what those trials it keeps add to
     its mean (their share times the difference of their mean from the other trials'). The
     Command Window shows the share of trials with each kind. A pipeline that
     removed the eye activity, by ICA or by rejecting those trials, adds about nothing. The
     estimate's own noise is taken off its square, so a bias the data do not show counts about 0.
     Second, what the signal's shape changes: the gain is measured with one known shape, and a
     filter changes the gain of a broad, slow component more than that of a narrow one. So each
     pipeline's filters are also run on bumps 0.5, 1 and 2 times as wide as the window, and the
     largest change of the gain (the worst case) times the measured mean counts as bias.
   - The rejection limit, ICLabel threshold and components chosen from the data minimize the same
     total error (without the shape part, which they do not change). The results table shows the
     total error (*error (noise + bias)*); the result says how much blinks, eye movements and muscle activity add in the recommended
     pipeline, how much its filters could change the measure for another shape, or that there is
     no eye channel to check blinks with. In the table of the ranking: `objective` (total error),
     `noise`, `bias`, and per measure `eyeBias_…` (µV, the condition where it is largest) and
     `shapeSpread_…`.
3. **Uncertainty.** A paired bootstrap over trials, with simultaneous intervals across all pairs
   of candidates (White, 2000; Romano & Wolf, 2005), finds the set of pipelines that the data
   cannot distinguish from the best.
4. **Equivalence.** *Not distinguished* can mean that the pipelines are equally good or that there
   are too few trials to tell. A pipeline is called *equivalent* to the best only when the upper
   bound of its difference from the best is within 5 % of the best's objective
   (`equivalenceMargin`, default 0.05); the bound is simultaneous over the pipelines not shown
   worse. The results say which of the two applies to the recommended pipeline.
5. **Selection check.** The best of many pipelines looks better than it is, because part of its
   advantage is chance. PipeCompare therefore also scores the choice on trials that were not used
   to make it: the trials of each condition are split into two halves at random, the pipeline
   that is best on one half is scored on the other half (rescaled to all trials), and the other
   way round, over 20 random splits (`nSplits`; consecutive halves for segments of one recording;
   0 turns the check off; at least 4 trials per condition are needed). The result reports how
   much better the best looks than it is on those trials. The pipelines' own decisions (bad
   channels, ICA, rejected epochs) were made on all trials; only the choice between pipelines is
   held out.
6. **Recommendation.** The pipelines are combinations of settings (high-pass × low-pass ×
   weights ...). Each setting's effect is estimated from every pipeline that has it, and what
   a combination adds beyond the effects of its settings is kept only as far as the data show it
   is more than noise (the rest is shrunk away: an empirical-Bayes estimate, with the noise from
   the bootstrap of step 3). Within the set of step 3, PipeCompare recommends the pipeline with
   the lowest of these estimates. With 30–40 trials almost every pipeline is in that set, and the
   lowest error of one pipeline alone is mostly chance: in simulation the estimate through the
   settings picked a pipeline with up to 60 % less excess error. Pipelines shown *equivalent* to the
   best (within 5 %) are equally good, and of those the least aggressive is used when it keeps
   more trials: the highest retention, then the least signal distortion. A pipeline whose bias is
   more than half its noise is recommended only when no pipeline with less bias is in that set:
   keeping more trials is no reason to keep the blinks they carry.
7. **Stability.** In each bootstrap resample of step 3, steps 3 and 6 are done again on the
   resampled values (with the same critical value), and PipeCompare counts how often each
   pipeline would be recommended (the bootstrap selection frequency; Efron, 2014). The table has
   it as `chosenShare`; the result reports it for the recommended pipeline. A low share means the
   data do not single that pipeline out.
8. **ICA's random start.** ICA (runica) starts from a fixed random state, so the same data
   always give the same components, but another start would give somewhat different ones, and
   with them other components removed. Every pipeline that has ICA is run again twice, with
   ICA seeing the samples in a different random order (seeds 1 and 2); every step after it
   decides again on its own. All pipelines use the same new decomposition, so pipelines that
   share the steps before ICA fit only one ICA per seed (*Standard*: two more ICAs in all).
   With each seed the ranking is made again and the recommendation compared. The Command
   Window shows how much their total error changed. On in the simple mode; in the advanced panel
   tick *Options… > Check ICA's random start*; in scripts pass `'icaSeeds', 2` to
   `pipecompare.PipeCompare.optimize`. The result is in `pipecompare_result.icaSeeds` (`ids`,
   `objective` with the search's value first, `recommended`, `same`). It costs one more ICA per
   seed for each different set of steps before ICA, plus the steps after it.

The ICLabel threshold and the epoch-rejection limit of the *Standard* recipe (and, with that option,
the ICA components chosen by their effect on the noise) are chosen inside each pipeline, before
this ranking, with the same SME, constraints and bootstrap. The check on trials not used to choose
makes these choices again on the half used to choose, so neither the choice between pipelines
nor the choices inside them flatter the result (see
[ICLabel threshold chosen from the data](#iclabel-threshold-chosen-from-the-data),
[ICA components chosen by their effect on the noise](#ica-components-chosen-by-their-effect-on-the-noise) and
[Rejection limit chosen from the data](#rejection-limit-chosen-from-the-data)).

Pipelines that differ in something that changes the measured quantity, such as the reference, are
ranked separately and never compared with each other. With several such groups there is no overall
recommendation: each group has its own, and `adopt` and `writeScript` then need a candidate id.
A multiverse summary reports, for each measure and condition, how much its value varies over the
feasible pipelines and which searched choice accounts for most of that variation (η²); it is
descriptive and not used in the ranking. The full derivations are in
[docs/METHODS.md](docs/METHODS.md).

## Outputs

- A ranked table of all pipelines: status, objective (total error: gain-corrected SME and bias), simultaneous interval
  of its difference from the best, trial retention, interpolation, signal-check results, and the
  reason for each exclusion.
- The recommended pipeline and why it was chosen.
- The full EEGLAB command sequence for every candidate.
- A MATLAB function that runs a chosen candidate's steps on any recording, deciding bad channels,
  components and rejected epochs from that recording's data (`writeScript`).
- Adoption of a candidate as a new EEGLAB dataset whose `EEG.history` replays it (`adopt`). For
  a pipeline that weights its epochs, the weights are in `EEG.etc.pipecompare.trialWeights`
  (see [Epochs weighted by their noise](#epochs-weighted-by-their-noise)).

## Limitations

- **The analysis contract is fixed by you.** Event types, epoch and baseline windows, regions of
  interest and measurement windows define what is measured, and are never searched.
- **The signal check is necessary, not sufficient.** The known signal has an assumed shape: a
  smooth bump as wide at half maximum as your measurement window (at least 50 ms), or in the simple
  mode your recording's own average over the window when it is reliable (see
  [Signal check with your recording's own waveform](#signal-check-with-your-recordings-own-waveform)),
  with a topography centred on the region of interest. Real components may be affected differently.
- **Some steps are re-run on the signal copy.** EEGLAB commands added as steps that make their
  own data-driven decisions, other than mark-and-remove workflows, are re-run on the copy carrying
  the known signal rather than replayed; this includes ASR added as an EEGLAB command. The built-in
  `asr` step replays the reconstructions ASR chose and is re-run only when they cannot be
  reproduced exactly. The result flags every re-run step.
- **Reference.** The dialog offers the recorded reference or the average reference. For another
  reference (e.g. linked mastoids), add a re-reference step with those channels in the panel or a
  script.
- **Channel locations are needed** for spherical interpolation and ICLabel. Without them those
  steps are excluded before the search, with the reason.
- **Peak latency is a weak criterion.** Its precision is hard to estimate with few trials, so it
  rarely separates pipelines; mean-amplitude measures are more informative.
- **Computation cost.** Every pipeline also runs on the signal copy, roughly doubling the EEGLAB
  computation. Memory use grows with plan depth; the expected peak is reported when it exceeds
  2 GB. Resampling early in the plan reduces both.
- **Single-dataset scope.** STUDY-level processing, time-frequency measures, connectivity and
  source analysis are not covered.

[docs/COVERAGE.md](docs/COVERAGE.md) lists every supported operation and how it has been
validated.

## Testing

```matlab
eeglab nogui                              % EEGLAB and its plugins must be on the path
addpath(fullfile(pwd, 'tests'));          % run from the PipeCompare folder
results = run_all();
```

The suite runs on synthetic data with known ground truth. Some of it is simulation, to show the
methods do not mislead: a P3 that differs from trial to trial in size (SD 50 %) and latency (SD
50 ms), which the SME counts as noise, is still measured without bias by every pipeline that
passes, weighted epochs included; a rejection limit chosen by chance is caught by the selection
check; a time-locked component ICLabel does not recognise is not removed; a stimulus-locked blink
is not taken for the P3's shape. It covers history parsing, plan
enumeration (checked against brute force), the statistics (validated by simulation), signal
preservation for every built-in step except `resample`, end-to-end searches with resume and
parallel execution, and EEGLAB and panel integration. Continuous integration runs every suite,
including those that open windows (on a virtual display); tests that need a toolbox CI does not
install, such as parallel execution (Parallel Computing Toolbox), are skipped there. Steps that
require clicking in EEGLAB's own dialogs cannot be automated and are listed in
[tests/MANUAL_GUI_CHECK.md](tests/MANUAL_GUI_CHECK.md).

**Ground-truth benchmark.** `benchmarks/pipecompare_benchmark.m` asks whether the pipeline
PipeCompare recommends measures the P3 more accurately than common fixed pipelines. It simulates
recordings (`benchmarks/bench_synth.m`: 32 channels, cortical background activity, alpha, drifts,
line noise, random blinks and blinks time-locked to the targets inside the P3 window, temporal
muscle bursts, movement artifacts and a bad channel) with a known P3, then measures the P3 (Pz,
300-600 ms) of targets and standards with the pipeline `pop_pipecompare` recommends, a minimal
pipeline (0.1-30 Hz, average reference, rejection at 100 µV) and a common one (the same with
bad-channel interpolation and ICA with ICLabel). It reports, per pipeline, the bias and the root
mean square error of each condition and of their difference against the P3 actually in the data,
the trials kept, and where the recommendation ranks by actual error among all the pipelines
PipeCompare compared. It also checks that PipeCompare's own error estimate tells the truth: for
every pipeline compared, the gain-corrected estimate's error against the population P3 is compared
with the ranked total error, including artifact and signal-shape bias. `conditionalError` also
retains the error against the realized trial mean, a distinct denoising target. Ratios above 1
mean reported error is understated. Records have equal weight, and confidence intervals resample
whole records. The selected pipelines are summarized separately; with/without learned-step
groups are descriptive, not controlled ablations. Small samples are marked exploratory and
missing recommendations are counted (`calibration.csv`, `benchmark_summary.txt`). It takes hours (an ICA per pipeline and recording), so on GitHub it runs
by hand only (Actions > *benchmark* > *Run workflow*; the results, written after every recording,
are kept as an artifact), or locally:

```matlab
eeglab nogui; addpath(fullfile(pwd, 'benchmarks'));
out = pipecompare_benchmark(8);           % 8 recordings; results also in tempdir
out = pipecompare_benchmark(30, tempdir, 1001); % disjoint validation seeds 1001..1030
```

Two optional checks use your own data, which is never committed: set the environment variable
`PIPECOMPARE_REAL_SET` to a `.set` file, or run `tests/realdata_validation.m`.

## Documentation

| Document | Contents |
|---|---|
| [docs/PANEL.md](docs/PANEL.md) | The dialog and the advanced panel, control by control |
| [docs/METHODS.md](docs/METHODS.md) | Measures, statistics and signal check, with derivations and references |
| [docs/COVERAGE.md](docs/COVERAGE.md) | Supported steps and their validation status |
| [CHANGELOG.md](CHANGELOG.md) | Release history |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Bug reports, documentation fixes and tests |

## Citation

If you use PipeCompare in published work, please cite it:

```bibtex
@software{pipecompare2026,
  author  = {abwoo},
  title   = {PipeCompare: data-driven comparison of EEG preprocessing pipelines in EEGLAB},
  year    = {2026},
  version = {0.9.5},
  url     = {https://github.com/abwoo/PipeCompare}
}
```

Please also cite the methods it builds on:

- Luck, S. J., Stewart, A. X., Simmons, A. M., & Rhemtulla, M. (2021). Standardized measurement
  error: A universal metric of data quality for averaged event-related potentials.
  *Psychophysiology*, 58(6), e13793.
- Kappenman, E. S., Farrens, J. L., Zhang, W., Stewart, A. X., & Luck, S. J. (2021). ERP CORE: An
  open resource for human event-related potential research. *NeuroImage*, 225, 117465.
- Zhang, G., Garrett, D. R., & Luck, S. J. (2024). Optimal filters for ERP research I: A general
  approach for selecting filter settings. *Psychophysiology*, 61(6), e14531.
- White, H. (2000). A reality check for data snooping. *Econometrica*, 68(5), 1097–1126.
- Romano, J. P., & Wolf, M. (2005). Stepwise multiple testing as formalized data snooping.
  *Econometrica*, 73(4), 1237–1282.

The pipelines run on EEGLAB (Delorme, A., & Makeig, S. (2004). *Journal of Neuroscience Methods*,
134(1), 9–21); when the chosen pipeline removes components with ICLabel, cite it too
(Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019). *NeuroImage*, 198, 181–197). The
complete reference list is in [docs/METHODS.md](docs/METHODS.md#references).

## License

PipeCompare is released under the [MIT License](LICENSE).
