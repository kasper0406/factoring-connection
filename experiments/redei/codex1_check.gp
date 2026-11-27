\\ Run with gp -fq check.gp.
smallnorm(n) = {
  my(bound = logint(n, 2) + 1);
  for(y = 0, bound,
    if(issquare(n - y^2) || issquare(n + 2*y^2), return(1))
  );
  return(0);
};

two_part(d) = {
  my(orders = quadclassunit(d)[2]);
  my(field = bnfinit(x^2 - d/4));
  if(!bnfcertify(field) || orders != field.cyc, error("unverified group for ", d));
  return(select(x -> x > 1, apply(x -> 2^valuation(x, 2), orders)));
};

for(i = 1, 2, {
  my(n = [367*751, 383*991][i]);
  print("Q3 ", n, " mod 256 = ", n%256, "   ", two_part(-4*n), "   ", two_part(-8*n));
});

{
  my(counts = matrix(2, 2), first_miss = 0);
  forprime(p = 101, 499, forprime(q = p+1, 499,
    my(n = p*q);
    if(n%8 != 1, next);
    my(h = qfbclassno(-4*n), actual = (h%8 == 0), guessed = smallnorm(n));
    if(actual != (p%8 != 3), error("4-rank criterion failed for ", n));
    my(reps = 0);
    forstep(y = 1, sqrtint(n\8), 2, reps += issquare(n - 8*y^2));
    if(reps%2 != 1-actual, error("Delta^9 coefficient mismatch for ", n));
    counts[1+actual, 1+guessed]++;
    if(actual && !guessed && !first_miss,
      first_miss = n;
      print("Q4 first missed input ", p, " * ", q, " = ", n, "   h = ", h)
    )
  ));
  print("Q4 rows = true bit 0,1, columns = output 0,1   ", counts);
}
