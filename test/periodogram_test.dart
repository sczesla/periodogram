import 'dart:math';

import 'package:periodogram/periodogram.dart';
import 'package:test/test.dart';

/// Unevenly sampled sine with offset and noise.
({List<double> t, List<double> y, List<double> err}) _data(
    {int n = 60, double period = 7.3, int seed = 1}) {
  final rnd = Random(seed);
  final t = [for (int i = 0; i < n; i++) 1000 + 50 * rnd.nextDouble()]..sort();
  final err = [for (int i = 0; i < n; i++) 0.2 + 0.3 * rnd.nextDouble()];
  final y = [
    for (int i = 0; i < n; i++)
      5 +
          1.5 * sin(2 * pi * t[i] / period + 0.4) +
          err[i] * (rnd.nextDouble() + rnd.nextDouble() - 1) * 2
  ];
  return (t: t, y: y, err: err);
}

/// Solves the linear system m x = r by Gaussian elimination.
List<double> _linSolve(List<List<double>> m, List<double> r) {
  final k = r.length;
  final a = [
    for (int i = 0; i < k; i++) [...m[i], r[i]]
  ];
  for (int c = 0; c < k; c++) {
    int piv = c;
    for (int i = c + 1; i < k; i++) {
      if (a[i][c].abs() > a[piv][c].abs()) piv = i;
    }
    final tmp = a[c];
    a[c] = a[piv];
    a[piv] = tmp;
    for (int i = 0; i < k; i++) {
      if (i == c) continue;
      final f = a[i][c] / a[c][c];
      for (int j = c; j <= k; j++) {
        a[i][j] -= f * a[c][j];
      }
    }
  }
  return [for (int i = 0; i < k; i++) a[i][k] / a[i][i]];
}

/// Reference: chi² of the weighted least-squares fit of the basis
/// functions to y.
double _chi2(List<double> y, List<double> w, List<List<double>> basis) {
  final k = basis.length, n = y.length;
  final m = [
    for (int i = 0; i < k; i++)
      [
        for (int j = 0; j < k; j++)
          [for (int l = 0; l < n; l++) w[l] * basis[i][l] * basis[j][l]]
              .reduce((a, b) => a + b)
      ]
  ];
  final r = [
    for (int i = 0; i < k; i++)
      [for (int l = 0; l < n; l++) w[l] * basis[i][l] * y[l]]
          .reduce((a, b) => a + b)
  ];
  final x = _linSolve(m, r);
  double chi2 = 0;
  for (int l = 0; l < n; l++) {
    double model = 0;
    for (int i = 0; i < k; i++) {
      model += x[i] * basis[i][l];
    }
    chi2 += w[l] * (y[l] - model) * (y[l] - model);
  }
  return chi2;
}

/// Reference power 1 - chi²(f) / chi²_0 by direct least squares.
double _refPower(List<double> t, List<double> y, List<double>? err, double f,
    {required bool fitMean}) {
  final n = t.length;
  final w = [
    for (int i = 0; i < n; i++) err == null ? 1.0 : 1 / (err[i] * err[i])
  ];
  final wSum = w.reduce((a, b) => a + b);
  double mean = 0;
  for (int i = 0; i < n; i++) {
    mean += w[i] * y[i] / wSum;
  }
  final yc = [for (final v in y) v - mean];
  double chi0 = 0;
  for (int i = 0; i < n; i++) {
    chi0 += w[i] * yc[i] * yc[i];
  }
  final cosB = [for (final ti in t) cos(2 * pi * f * ti)];
  final sinB = [for (final ti in t) sin(2 * pi * f * ti)];
  final chi = fitMean
      ? _chi2(y, w, [cosB, sinB, List.filled(n, 1.0)])
      : _chi2(yc, w, [cosB, sinB]);
  return 1 - chi / chi0;
}

double _chi2Of(List<double> y, List<double> model, List<double>? err) {
  double s = 0;
  for (int i = 0; i < y.length; i++) {
    final e = err == null ? 1.0 : err[i];
    s += (y[i] - model[i]) * (y[i] - model[i]) / (e * e);
  }
  return s;
}

void main() {
  final d = _data();
  final freqs = [for (int k = 1; k <= 40; k++) k * 0.012];

  group('power against direct least squares', () {
    for (final fitMean in [true, false]) {
      for (final withErr in [true, false]) {
        test('fitMean: $fitMean, errors: $withErr', () {
          final err = withErr ? d.err : null;
          final p =
              gls(d.t, d.y, err: err, fitMean: fitMean, frequencies: freqs);
          for (int k = 0; k < freqs.length; k++) {
            expect(
                p.power[k],
                closeTo(_refPower(d.t, d.y, err, freqs[k], fitMean: fitMean),
                    1e-9));
          }
        });
      }
    }

    test('lombScargle equals gls without errors and floating mean', () {
      final ls = lombScargle(d.t, d.y, frequencies: freqs);
      final g = gls(d.t, d.y, fitMean: false, frequencies: freqs);
      expect(ls.power, g.power);
      expect(ls.fitMean, isFalse);
    });

    test('lombScargle matches the textbook formula with tau', () {
      final ls = lombScargle(d.t, d.y, frequencies: freqs);
      final n = d.t.length;
      final mean = d.y.reduce((a, b) => a + b) / n;
      final yc = [for (final v in d.y) v - mean];
      final variance = yc.map((v) => v * v).reduce((a, b) => a + b) / n;
      for (int k = 0; k < freqs.length; k++) {
        final w = 2 * pi * freqs[k];
        double s2 = 0, c2 = 0;
        for (final ti in d.t) {
          s2 += sin(2 * w * ti);
          c2 += cos(2 * w * ti);
        }
        final tau = atan2(s2, c2) / (2 * w);
        double yc1 = 0, ys1 = 0, cc1 = 0, ss1 = 0;
        for (int i = 0; i < n; i++) {
          final x = w * (d.t[i] - tau);
          yc1 += yc[i] * cos(x);
          ys1 += yc[i] * sin(x);
          cc1 += cos(x) * cos(x);
          ss1 += sin(x) * sin(x);
        }
        final pScargle = (yc1 * yc1 / cc1 + ys1 * ys1 / ss1) / (2 * variance);
        expect(ls.power[k] * n / 2, closeTo(pScargle, 1e-8));
      }
    });
  });

  group('best-fit sine', () {
    test('the sine reproduces the chi² of the fit, with floating mean', () {
      final p = gls(d.t, d.y, err: d.err, frequencies: freqs);
      for (final k in [3, p.bestIndex, 30]) {
        final fit = p.sineFitAt(k);
        final chi2 = _chi2Of(d.y, fit.evaluate(d.t), d.err);
        expect(chi2, closeTo(p.normalized(Normalization.chiSquared)[k], 1e-8));
      }
    });

    test('without floating mean, the offset is the weighted mean', () {
      final p = gls(d.t, d.y, err: d.err, fitMean: false, frequencies: freqs);
      for (final k in [3, p.bestIndex, 30]) {
        final fit = p.sineFitAt(k);
        expect(fit.offset, p.weightedMean);
        final chi2 = _chi2Of(d.y, fit.evaluate(d.t), d.err);
        expect(chi2, closeTo(p.normalized(Normalization.chiSquared)[k], 1e-8));
      }
    });

    test('noiseless sine: power 1 and exact parameters', () {
      final t = d.t;
      final y = [for (final ti in t) 2 + 0.7 * sin(2 * pi * (ti - 1003.1) / 4)];
      final p = gls(t, y, frequencies: [0.2, 0.25, 0.3]);
      expect(p.bestIndex, 1);
      final b = p.best;
      expect(b.power, closeTo(1, 1e-12));
      expect(b.amplitude, closeTo(0.7, 1e-10));
      expect(b.offset, closeTo(2, 1e-10));
      expect(b.period, 4);
      // t0 is a zero crossing of the sine: 1003.1 modulo the period
      final cycles = (b.t0 - 1003.1) / 4;
      expect(cycles - cycles.round(), closeTo(0, 1e-10));
      expect(b.phase, inInclusiveRange(-0.5, 0.5));
      expect(b.residualRms, closeTo(0, 1e-6));
    });

    test('period found on the default grid, with error', () {
      final p = gls(d.t, d.y, err: d.err);
      final b = p.best;
      expect(b.period, closeTo(7.3, 0.1));
      expect(b.amplitude, closeTo(1.5, 0.15));
      expect(b.offset, closeTo(5, 0.1));
      expect(b.frequencyError, greaterThan(0));
      expect(b.periodError,
          closeTo(b.frequencyError * b.period * b.period, 1e-12));
      expect(b.amplitudeError, greaterThan(0));
      expect(p.bestFap, lessThan(1e-10));
    });

    test('frequency error is NaN at the edge of the grid', () {
      final p = gls(d.t, d.y, frequencies: freqs);
      expect(p.sineFitAt(0).frequencyError, isNaN);
      expect(p.sineFitAt(freqs.length - 1).frequencyError, isNaN);
      expect(() => p.sineFitAt(freqs.length), throwsRangeError);
    });
  });

  group('normalizations', () {
    final p = gls(d.t, d.y, err: d.err, frequencies: freqs);
    final n = d.t.length;
    final z = p.power;
    double chi0 = 0;
    for (int i = 0; i < n; i++) {
      final r = (d.y[i] - p.weightedMean) / d.err[i];
      chi0 += r * r;
    }

    test('all conversions', () {
      final pMax = z.reduce(max);
      final sumLn =
          [for (final e in d.err) log(2 * pi * e * e)].reduce((a, b) => a + b);
      for (int k = 0; k < z.length; k++) {
        expect(p.normalized(Normalization.zk)[k], z[k]);
        expect(p.normalized(Normalization.horneBaliunas)[k],
            closeTo((n - 1) / 2 * z[k], 1e-12));
        expect(p.normalized(Normalization.cumming)[k],
            closeTo((n - 3) / 2 * z[k] / (1 - pMax), 1e-9));
        final chi2 = p.normalized(Normalization.chiSquared)[k];
        expect(chi2, closeTo(chi0 * (1 - z[k]), 1e-8));
        expect(p.normalized(Normalization.lnL)[k],
            closeTo(-0.5 * chi2 - 0.5 * sumLn, 1e-8));
        expect(p.normalized(Normalization.dlnL)[k],
            closeTo(0.5 * chi0 * z[k], 1e-8));
        expect(p.normalized(Normalization.wrms)[k],
            closeTo(p.weightedRms * sqrt(1 - z[k]), 1e-12));
      }
    });

    test('without errors, chiSquared is the sum of squared residuals', () {
      final q = gls(d.t, d.y, frequencies: freqs);
      final fit = q.best;
      expect(q.normalized(Normalization.chiSquared)[q.bestIndex],
          closeTo(_chi2Of(d.y, fit.evaluate(d.t), null), 1e-8));
    });
  });

  group('false-alarm probability', () {
    final p = gls(d.t, d.y, err: d.err);
    final n = d.t.length;

    test('probability and fap formulas', () {
      final m = p.independentFrequencies;
      for (final z in [0.0, 0.01, 0.1, 0.3, 0.5]) {
        final prob = pow(1 - z, (n - 3) / 2).toDouble();
        expect(p.probability(z), closeTo(prob, 1e-14));
        expect(p.fap(z), closeTo(1 - pow(1 - prob, m), 1e-12));
      }
      expect(p.fap(0), 1);
    });

    test('powerLevel is the inverse of fap', () {
      for (final f in [0.5, 0.1, 0.01, 1e-3, 1e-8]) {
        expect(p.fap(p.powerLevel(f)), closeTo(f, f * 1e-9));
      }
      expect(() => p.powerLevel(0), throwsArgumentError);
      expect(() => p.powerLevel(1), throwsArgumentError);
    });

    test('tiny probabilities stay accurate', () {
      // M * prob is far below the double resolution of 1 - prob
      final z = 0.9;
      final prob = pow(0.1, (n - 3) / 2).toDouble();
      expect(
          p.fap(z),
          closeTo(p.independentFrequencies * prob,
              1e-6 * prob * p.independentFrequencies));
    });

    test('number of independent frequencies', () {
      final tSpan = d.t.last - d.t.first;
      final q = gls(d.t, d.y, frequencies: freqs);
      expect(q.independentFrequencies,
          closeTo((freqs.last - freqs.first) * tSpan, 1e-9));
      final single = gls(d.t, d.y, frequencies: [0.1]);
      expect(single.independentFrequencies, 1);
    });
  });

  group('frequency grid', () {
    test('default: from one step up to the average Nyquist frequency', () {
      final p = gls(d.t, d.y);
      final tSpan = d.t.last - d.t.first;
      final step = 1 / (tSpan * 10);
      final nyquist = d.t.length / (2 * tSpan);
      expect(p.frequencies.first, closeTo(step, 1e-15));
      expect(p.frequencies[1] - p.frequencies[0], closeTo(step, 1e-12));
      expect(p.frequencies.last, lessThan(nyquist));
      expect(p.frequencies.last + step, greaterThanOrEqualTo(nyquist - 1e-12));
      expect(p.independentFrequencies, closeTo((nyquist - step) * tSpan, 1e-9));
    });

    test('fMin, fMax, oversampling, nyquistFactor', () {
      final tSpan = d.t.last - d.t.first;
      final p = gls(d.t, d.y, fMin: 0.1, fMax: 0.2, oversampling: 5);
      expect(p.frequencies.first, 0.1);
      expect(p.frequencies.last, lessThan(0.2));
      expect(
          p.frequencies[1] - p.frequencies[0], closeTo(1 / (5 * tSpan), 1e-12));
      final q = gls(d.t, d.y, nyquistFactor: 2);
      expect(q.frequencies.last, greaterThan(d.t.length / (2 * tSpan)));
    });

    test('invalid grids throw', () {
      expect(() => gls(d.t, d.y, frequencies: []), throwsArgumentError);
      expect(() => gls(d.t, d.y, frequencies: [0.1, 0]), throwsArgumentError);
      expect(() => gls(d.t, d.y, frequencies: [0.1], fMin: 0.05),
          throwsArgumentError);
      expect(() => gls(d.t, d.y, fMin: 0.2, fMax: 0.1), throwsArgumentError);
      expect(() => gls(d.t, d.y, fMin: -1), throwsArgumentError);
      expect(() => gls(d.t, d.y, oversampling: 0), throwsArgumentError);
      expect(
          () => gls(d.t, d.y, nyquistFactor: double.nan), throwsArgumentError);
    });
  });

  group('degenerate frequencies with even sampling', () {
    final t = [for (int i = 0; i < 20; i++) i.toDouble()];
    final y = [for (final ti in t) sin(2 * pi * ti / 5) + 0.3 * cos(pi * ti)];

    test('multiple of the sampling frequency: power 0, no NaN', () {
      final p = gls(t, y, frequencies: [1.0, 2.0]);
      expect(p.power, [0.0, 0.0]);
      expect(p.best.amplitude, 0);
    });

    test('Nyquist frequency: fit with the cosine alone', () {
      final p = gls(t, y, frequencies: [0.5]);
      final n = t.length;
      final ones = List.filled(n, 1.0);
      final cosB = [for (final ti in t) cos(pi * ti)];
      final w = List.filled(n, 1.0);
      final mean = y.reduce((a, b) => a + b) / n;
      final chi0 =
          y.map((v) => (v - mean) * (v - mean)).reduce((a, b) => a + b);
      final expected = 1 - _chi2(y, w, [cosB, ones]) / chi0;
      expect(p.power.single, closeTo(expected, 1e-12));
      expect(_chi2Of(y, p.best.evaluate(t), null),
          closeTo(chi0 * (1 - expected), 1e-10));
    });
  });

  group('input validation', () {
    test('invalid data throw', () {
      expect(() => gls([1, 2, 3], [1, 2, 3]), throwsArgumentError);
      expect(() => gls([1, 2, 3, 4], [1, 2, 3]), throwsArgumentError);
      expect(() => gls([1, 2, 3, 4], [1, 2, 3, 4], err: [1, 1, 1]),
          throwsArgumentError);
      expect(
          () => gls([1, 2, 3, 4], [1, 2, double.nan, 4]), throwsArgumentError);
      expect(() => gls([1, 2, 3, 4], [1, 2, 3, 4], err: [1, 0, 1, 1]),
          throwsArgumentError);
      expect(() => gls([1, 1, 1, 1], [1, 2, 3, 4]), throwsArgumentError);
      expect(() => gls([1, 2, 3, 4], [2, 2, 2, 2]), throwsArgumentError);
    });

    test('four points are enough', () {
      final p = gls([0, 1, 2.5, 4], [1, 3, 2, 0]);
      expect(p.n, 4);
      expect(p.power.every((v) => v >= 0 && v <= 1 + 1e-12), isTrue);
    });
  });

  group('spectralWindow', () {
    test('W(0) = 1 and aliases of even sampling', () {
      final t = [for (int i = 0; i < 30; i++) i.toDouble()];
      final w = spectralWindow(t, [0, 1, 2, 0.5, 0.25]);
      expect(w[0], closeTo(1, 1e-12));
      expect(w[1], closeTo(1, 1e-12));
      expect(w[2], closeTo(1, 1e-12));
      expect(w[3], closeTo(0, 1e-12));
      expect(w[4], lessThan(0.01));
    });

    test('weights and the direct formula', () {
      final f = [0.05, 0.13, 0.4];
      final w = spectralWindow(d.t, f, err: d.err);
      for (int k = 0; k < f.length; k++) {
        double re = 0, im = 0, ws = 0;
        for (int i = 0; i < d.t.length; i++) {
          final wi = 1 / (d.err[i] * d.err[i]);
          re += wi * cos(2 * pi * f[k] * d.t[i]);
          im += wi * sin(2 * pi * f[k] * d.t[i]);
          ws += wi;
        }
        expect(w[k], closeTo((re * re + im * im) / (ws * ws), 1e-12));
      }
    });

    test('invalid input throws', () {
      expect(() => spectralWindow([], [0.1]), throwsArgumentError);
      expect(
          () => spectralWindow([1, 2], [0.1], err: [1]), throwsArgumentError);
      expect(
          () => spectralWindow([1, 2], [double.infinity]), throwsArgumentError);
    });
  });
}
