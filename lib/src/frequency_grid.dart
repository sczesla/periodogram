import 'dart:math';

/// Frequencies for a periodogram and the estimated number of independent
/// frequencies among them (internal; see `gls` for the meaning of the
/// arguments).
({List<double> frequencies, double independentFrequencies}) frequencyGrid({
  required int n,
  required double timeSpan,
  List<double>? frequencies,
  double? fMin,
  double? fMax,
  required double oversampling,
  required double nyquistFactor,
}) {
  double fBeg, fEnd;
  List<double> freq;
  if (frequencies != null) {
    if (fMin != null || fMax != null) {
      throw ArgumentError("Give either frequencies or fMin/fMax, not both");
    }
    if (frequencies.isEmpty) {
      throw ArgumentError("frequencies must not be empty");
    }
    for (final f in frequencies) {
      if (!(f > 0 && f.isFinite)) {
        throw ArgumentError("frequencies must be finite and > 0 (got $f)");
      }
    }
    freq = List.of(frequencies);
    fBeg = freq.reduce(min);
    fEnd = freq.reduce(max);
  } else {
    if (!(oversampling > 0 && oversampling.isFinite)) {
      throw ArgumentError.value(
        oversampling,
        "oversampling",
        "must be finite and > 0",
      );
    }
    if (!(nyquistFactor > 0 && nyquistFactor.isFinite)) {
      throw ArgumentError.value(
        nyquistFactor,
        "nyquistFactor",
        "must be finite and > 0",
      );
    }
    final step = 1 / (timeSpan * oversampling);
    fBeg = fMin ?? step;
    fEnd = fMax ?? nyquistFactor * 0.5 * n / timeSpan;
    if (!(fBeg > 0 && fBeg.isFinite)) {
      throw ArgumentError.value(fBeg, "fMin", "must be finite and > 0");
    }
    if (!(fEnd > fBeg && fEnd.isFinite)) {
      throw ArgumentError(
        "The maximum frequency ($fEnd) must be finite and "
        "larger than the minimum frequency ($fBeg)",
      );
    }
    final count = ((fEnd - fBeg) / step).ceil();
    freq = [for (int k = 0; k < count; k++) fBeg + k * step];
    if (freq.isNotEmpty && freq.last >= fEnd) freq.removeLast();
    if (freq.isEmpty) freq = [fBeg];
  }
  return (
    frequencies: freq,
    independentFrequencies: max(1.0, (fEnd - fBeg) * timeSpan),
  );
}
