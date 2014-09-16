import sys

INF = 2^31

# returns differential ring and variables up to order n
def diff_ring(n=30, F=QQ):
  R = PolynomialRing(F, n, 'y') #, order='invlex')
  return R, R.gens()

# Derivative of order n of the differential polynomial f
def df(f, n = 1):
  if n == 0:
    return f
  R = f.parent()
  var = R.gens()
  result = R(0)
  for i in xrange(len(var)-1):
    result += diff(f, var[i]) * var[i+1]
  if n > 1:
    return df(result, n-1)
  else:
    return result

# order of the differential polynomial
def diff_ord(f):
  R = f.parent()
  var = R.gens()
  for i in reversed(xrange(len(var))):
    if diff(f, var[i]) != 0:
      return i
  return -1

# separant
def S(f):
  n = diff_ord(f)
  if n < 0:
    return 0
  R = f.parent()
  var = R.gens()
  return diff(f, var[n])

# returns ideals G_0, ..., G_p
def FindGk(f):
  R = f.parent()
  var = R.gens()
  l = diff_ord(f)
  F = [f]
  max_k = 2
  for i in xrange(max_k):
    F.append(df(F[-1]))
  G = []
  k = 0
  while True:
    Gp = []
    for j in xrange(max(0, l-k), l+1):
      Gp.append(df(diff(f, var[j]), k-l+j))
    G.append(Gp)
    F += Gp
    if 1 in F*R:
      break
    k += 1
  return max_k, G

def compute_gkn(f, l, k, n):
  res = 0
  for j in xrange(l+1):
    if j + k - l < 0:
      continue
    res += binomial(n, j + k - l) * df(diff(f,y[j]), j + k - l)
  return res

def contains(G, ideal):
  for g in G:
    if g not in ideal:
      return False
  return True

def compute_psi(gn, ideal):
  G = ideal.groebner_basis()
  P = PolynomialRing(QQ, gn.parent().variable_names_recursive())
  g = P(gn).reduce(G)
  coef = gcd(gn.parent()(g).coefficients())
  if coef == 0:
    return INF
  elif coef in QQ:
    return -INF
  else:
    roots = filter(lambda x: x in ZZ, map(lambda x: x[0], coef.roots()))
    if roots:
      return max(roots)
    else:
      return -INF

def compute_Ik(f, k):
  F = [f]
  for i in xrange(k):
    F.append(df(F[-1]))
  return F*f.parent()

def IsTrivial(f, H):
  R = f.parent()
  l = diff_ord(f)
  k, G = FindGk(f)
  p = len(G)
  for d in xrange(p):
    Ik = compute_Ik(f, k) + H*R
    Ik = Ik.radical()
    Ikd = Ik.quotient(G[d]*R)
    if 1 in Ikd:
      continue
    P.<n> = QQ[]
    T = PolynomialRing(P, R.gens())
    gkn = compute_gkn(f, l, d, T(n))
    Pkd = Ikd.primary_decomposition_complete()
    K = k
    for primary, prime in Pkd:
      a = compute_psi(gkn, prime)
      if a < k:
        return False
      K = max(K, a)
    k = K + 1
  return True

R, y = diff_ring(10, QQ)

f = y[0]*y[1]^2 + y[1]^2 + y[0]*y[1] - 1
h = S(f)
print IsTrivial(f, [h])

f = y[0]*y[1]^2 + y[1]^2 + y[0]*y[1] + 1
h = S(f)
print IsTrivial(f, [h])

f = (y[1]-1)^2 - 1^2*y[0]^4*y[1]^2
h = S(f)
print IsTrivial(f, [h])

