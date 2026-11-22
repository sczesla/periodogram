import 'dart:math';

/// Spectral window of the sampling times [t] at the given [frequencies]:
///
///     W(f) = | sum_i w_i exp(2 pi i f t_i) |²
///
/// with weights `w_i` proportional to `1 / err_i²` (all equal without
/// [err]) and normalized to sum 1, so that W(0) = 1.
///
/// It is the periodogram of the sampling alone. Peaks in W at a frequency
/// `f_w` (e.g., 1/day for nightly or daily observations, 1/week for
/// weekly readings) cause aliases at `f ± f_w` of every true signal at
/// frequency f. Compare W with the periodogram, e.g. on
/// `periodogram.frequencies`, to recognize them.
///
/// Throws an [ArgumentError] if the lengths differ, [t] is empty, or
/// values are not finite (errors also not > 0).
List<double> spectralWindow(
  List<double> t,
  List<double> frequencies, {
  List<double>? err,
}) {
  if (t.isEmpty) throw ArgumentError("t must not be empty");
  if (err != null && err.length != t.length) {
    throw ArgumentError(
      "t and err must have the same length "
      "(got ${t.length}, ${err.length})",
    );
  }
  for (int i = 0; i < t.length; i++) {
    if (!t[i].isFinite) throw ArgumentError("t must be finite (index $i)");
    if (err != null && !(err[i] > 0 && err[i].isFinite)) {
      throw ArgumentError("err must be finite and > 0 (index $i)");
    }
  }
  for (final f in frequencies) {
    if (!f.isFinite) throw ArgumentError("frequencies must be finite");
  }
  final w = [
    for (int i = 0; i < t.length; i++)
      err == null ? 1.0 : 1 / (err[i] * err[i]),
  ];
  final wSum = w.fold(0.0, (s, v) => s + v);
  final t0 = t.reduce(min);
  return List.unmodifiable([
    for (final f in frequencies)
      () {
        double re = 0, im = 0;
        for (int i = 0; i < t.length; i++) {
          final x = 2 * pi * f * (t[i] - t0);
          re += w[i] * cos(x);
          im += w[i] * sin(x);
        }
        re /= wSum;
        im /= wSum;
        return re * re + im * im;
      }(),
  ]);
}
