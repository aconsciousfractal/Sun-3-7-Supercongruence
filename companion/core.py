"""Exact arithmetic for the Sun 3.7 companion (standard library).
Polynomials are ascending coefficient lists. No conjectural value is used here.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import comb, gcd, isqrt
import json
from pathlib import Path


def need(test, message):
    if not test:
        raise ArithmeticError(message)


def natural(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name+' must be an integer >= '+str(minimum))
    return value


def save(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8', newline='\n')


def primes(bound):
    natural(bound, 'prime bound')
    return [n for n in range(2, bound+1) if all(n % k for k in range(2, isqrt(n)+1))]


def add(*arrays):
    out = [0] * max(map(len, arrays))
    for a in arrays:
        for i, v in enumerate(a): out[i] += v
    return out


def scale(a, c): return [c*v for v in a]


def mul(a, b, n):
    out = [0]*(n+1)
    for i, x in enumerate(a[:n+1]):
        if x:
            for j, y in enumerate(b[:n-i+1]):
                if y: out[i+j] += x*y
    return out


def mpow(a, k, n):
    natural(k, 'series exponent'); natural(n, 'series precision')
    out = [1]+[0]*n
    while k:
        if k & 1: out = mul(out, a, n)
        k >>= 1
        if k: a = mul(a, a, n)
    return out


def inv(a, n):
    need(a[0] in (1, -1), 'inverse requires constant +/-1')
    out = [a[0]]+[0]*n
    for k in range(1, n+1):
        out[k] = -a[0]*sum(a[j]*out[k-j] for j in range(1, min(k+1, len(a))))
    return out


def theta(a): return [i*v for i,v in enumerate(a)]


def compose(a, t, n):
    need(t[0] == 0, 'composition not at zero')
    out = [0]*(n+1)
    for c in reversed(a):
        out = mul(out, t, n); out[0] += c
    return out


def ev(a, x):
    out = 0
    for c in reversed(a): out = out*x+c
    return out


def modmul(a, b, p): return [x % p for x in mul(a,b,len(a)+len(b)-2)]


def modpow(a, k, p):
    natural(k, 'polynomial exponent')
    out = [1]
    while k:
        if k&1: out=modmul(out,a,p)
        k >>= 1
        if k:a=modmul(a,a,p)
    return out


def residue(q, modulus):
    if type(q) is not int and not isinstance(q, F):
        raise ValueError('residue requires an integer or Fraction')
    natural(modulus, 'modulus', 2)
    q = F(q)
    return q.numerator*pow(q.denominator,-1,modulus)%modulus


def sequence(a,b,c,top):
    natural(top, 'sequence top index')
    if top == 0:return [1]
    out=[1,b]
    for n in range(1,top):
        numerator=(a*n*(n+1)+b)*out[n]-c*n*n*out[n-1]
        need(numerator%((n+1)**2)==0, 'nonintegral recurrence')
        out.append(numerator//((n+1)**2))
    return out


def central(u):return [comb(2*n,n)*v for n,v in enumerate(u)]


def sum_mod(g,p,t,exponent=2):
    modulus=p**exponent; t=residue(t,modulus); out=0
    for c in reversed(g[:p]):out=(out*t+c)%modulus
    return out


def leg(a,p):
    b=pow(a%p,(p-1)//2,p)
    return -1 if b==p-1 else b


def eta(rs,n,engine='product',remove_shift=False):
    shift=F(sum(d*r for d,r in rs.items()),24)
    if remove_shift:shift=F(0)
    need(shift.denominator==1 and shift>=0, 'nonnegative integral q shift required')
    shift=int(shift); N=n-shift
    if N<0:return [0]*(n+1)
    if engine=='product':
        out=[1]+[0]*N
        for k in range(1,N+1):
            e=sum(r for d,r in rs.items() if k%d==0)
            if not e:continue
            coeff=[((-1)**j*comb(e,j) if j<=e else 0) if e>=0 else comb(-e+j-1,j) for j in range(N//k+1)]
            nxt=[0]*(N+1)
            for j,v in enumerate(coeff):
                if v:
                    s=k*j
                    for i in range(N-s+1):nxt[i+s]+=out[i]*v
            out=nxt
    elif engine=='log':
        L=log_eta(rs,N);L[0]=0
        out=[1]+[0]*N
        for k in range(1,N+1):
            v=sum(L[j]*out[k-j] for j in range(1,k+1))
            need(v%k==0,'log eta exact division');out[k]=v//k
    else:raise ValueError('unknown eta engine')
    return [0]*shift+out


def log_eta(rs,n):
    out=[F(sum(d*r for d,r in rs.items()),24)]+[0]*n
    for k in range(1,n+1):
        e=sum(r for d,r in rs.items() if k%d==0)
        for j in range(k,n+1,k):out[j]-=k*e
    return [int(x) if F(x).denominator==1 else x for x in out]


def eta_metadata(rs,N=6):
    weight=F(sum(rs.values()),2)
    x=F(1)
    for d,r in rs.items():x*=F(d)**r
    square=isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator
    return {'exponents':rs,'weight':str(weight),
            'sum_dr':sum(d*r for d,r in rs.items()),
            'sum_Nd_r':sum((N//d)*r for d,r in rs.items()),
            'trivial_character':weight.denominator==1 and weight%2==0 and square,
            'cusp_orders':{str(c):str(sum(F(N*gcd(c,d)**2*r,24*gcd(c,N//c)*c*d) for d,r in rs.items())) for c in range(1,N+1) if N%c==0}}

H={1:-5,2:1,3:-1,6:5}
Z={1:6,2:-3,3:-2,6:1}
U={1:12,2:-6,3:-4,6:2}
B={1:7,2:-5,3:-5,6:7}
V={1:2,2:-4,3:-6,6:12}
J={1:2,2:2,3:2,6:2}
R={1:-12,2:12,3:12,6:-12}
H2={1:-24,2:24}
S={1:-6,2:6,3:-6,6:6}
PARAMS={'a':(10,3,9),'Q':(-17,-6,72),'f':(7,2,-8)}


def modular_data(n,engine='log'):
    uu=eta(U,n,engine);bb=eta(B,n,engine);vv=eta(V,n,engine);jj=eta(J,n,engine)
    forms={'a':add(uu,scale(bb,18),scale(vv,72)), 'Q':add(uu,scale(vv,-72))}
    ts={name:mul(jj,mpow(inv(th,n),2,n),n) for name,th in forms.items()}
    return {'h':eta(H,n,engine),'J':jj,'forms':forms,'parameters':ts}
