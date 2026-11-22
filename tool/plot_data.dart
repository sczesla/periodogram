// Computes the data for the figures in doc/images with the package itself
// and writes them to tool/plot_data.json; tool/plots.py draws them.
//
//     dart run tool/plot_data.dart && python3 tool/plots.py

import 'dart:convert';
import 'dart:io';

import 'package:periodogram/periodogram.dart';

import '../example/periodogram_example.dart' show simulatedData, truePeriods;

void main() {
  final (:t, :y, :err) = simulatedData();
  // A finer grid than the default for smooth curves in the figures
  const oversampling = 30.0;
  final p = gls(t, y, err: err, oversampling: oversampling);
  final best = p.best;

  final tFine = [for (int i = 0; i <= 2000; i++) -0.2 + i * 10 / 2000];
  final model = best.evaluate(t);
  final residuals = [for (int i = 0; i < t.length; i++) y[i] - model[i]];
  final p2 = gls(t, residuals, err: err, oversampling: oversampling);

  final out = {
    "t": t,
    "y": y,
    "err": err,
    "residuals": residuals,
    "tFine": tFine,
    "modelFine": best.evaluate(tFine),
    "bestPeriod": best.period,
    "bestT0": best.t0,
    "trueFrequencies": [for (final period in truePeriods) 1 / period],
    "frequencies": p.frequencies,
    "gls": p.power,
    "glsLevels": [p.powerLevel(0.01), p.powerLevel(0.001)],
    "window": spectralWindow(t, p.frequencies),
    "residualGls": p2.power,
    "residualLevels": [p2.powerLevel(0.01), p2.powerLevel(0.001)],
    "secondPeriod": p2.best.period,
    "secondT0": p2.best.t0,
    "secondModelFine": p2.best.evaluate(tFine),
  };
  File("tool/plot_data.json").writeAsStringSync(jsonEncode(out));
  print("Written tool/plot_data.json");
}
