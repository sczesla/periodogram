import 'dart:math';

import 'package:periodogram/periodogram.dart';

/// The periods of the two sines in [simulatedData], in days.
const truePeriods = [1.8, 2.9];

/// Simulated measurements: 40 points in 8 nights within 10 days (a gap from
/// day 5 to 7), each night within 12 hours, with individual uncertainties
/// between 0.1 and 0.5. The signal consists of two sines with periods of
/// 1.8 and 2.9 days ([truePeriods]) and amplitudes of 1.3 and 0.5.
({List<double> t, List<double> y, List<double> err}) simulatedData(
    {int seed = 15}) {
  final rnd = Random(seed);
  const nights = [0, 1, 2, 3, 4, 7, 8, 9];
  final t = [
    for (int i = 0; i < 40; i++)
      nights[rnd.nextInt(nights.length)] + 0.5 * rnd.nextDouble()
  ]..sort();
  final err = [for (final _ in t) 0.1 + 0.4 * rnd.nextDouble()];
  // Normally distributed noise (Box-Muller)
  double gauss() =>
      sqrt(-2 * log(1 - rnd.nextDouble())) * cos(2 * pi * rnd.nextDouble());
  final y = [
    for (int i = 0; i < t.length; i++)
      3.5 +
          1.3 * sin(2 * pi * t[i] / truePeriods[0] + 0.75) +
          0.5 * sin(2 * pi * t[i] / truePeriods[1] + 3.1) +
          err[i] * gauss()
  ];
  return (t: t, y: y, err: err);
}

void main() {
  final (:t, :y, :err) = simulatedData();

  // Generalized Lomb-Scargle periodogram on the default frequency grid
  // (up to the average Nyquist frequency, here about 2.1 per day)
  final p = gls(t, y, err: err);
  final best = p.best;
  print("Best period: ${best.period.toStringAsFixed(3)} "
      "± ${best.periodError.toStringAsFixed(3)} days");
  print("Amplitude: ${best.amplitude.toStringAsFixed(2)} "
      "± ${best.amplitudeError.toStringAsFixed(2)}");
  print("False-alarm probability: ${p.bestFap.toStringAsExponential(0)}");
  print("Power level for FAP = 1 %: ${p.powerLevel(0.01).toStringAsFixed(2)}");

  // The classical Lomb-Scargle periodogram ignores the uncertainties and
  // finds the signal with less significance
  final ls = lombScargle(t, y);
  print("Lomb-Scargle: best period ${ls.best.period.toStringAsFixed(3)} "
      "days, false-alarm probability "
      "${ls.bestFap.toStringAsExponential(0)}");

  // The nightly sampling causes a peak of the spectral window at 1/day,
  // so every signal at f has aliases at |f ± 1|; here at
  // 1 - 1/1.8 = 0.444/day (2.25 days)
  final w = spectralWindow(t, [1.0]);
  print("Spectral window at 1/day: ${w.single.toStringAsFixed(2)}");

  // Subtract the best-fit sine and search again: the weaker signal
  final model = best.evaluate(t);
  final residuals = [for (int i = 0; i < t.length; i++) y[i] - model[i]];
  final p2 = gls(t, residuals, err: err);
  print("Second period: ${p2.best.period.toStringAsFixed(2)} days, "
      "amplitude ${p2.best.amplitude.toStringAsFixed(2)}, "
      "false-alarm probability ${p2.bestFap.toStringAsExponential(0)}");
}
