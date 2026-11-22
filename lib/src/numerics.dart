import 'dart:math';

// Accurate ln(1 + x) and exp(x) - 1 for small x, which dart:math lacks
// (Kahan's tricks; see Goldberg 1991, "What every computer scientist
// should know about floating-point arithmetic").

/// ln(1 + x), accurate also for |x| << 1.
double log1p(double x) {
  final u = 1 + x;
  if (u == 1) return x;
  return log(u) * x / (u - 1);
}

/// exp(x) - 1, accurate also for |x| << 1.
double expm1(double x) {
  final u = exp(x);
  if (u == 1) return x;
  final um1 = u - 1;
  if (um1 == -1) return -1;
  return um1 * x / log(u);
}
