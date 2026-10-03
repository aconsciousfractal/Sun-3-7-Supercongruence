"""Independent exact finite checks for Sun 3.7.
Pure standard library. No import of the recurrence or rational-pair engines.
Full-polynomial identities, original binomial sums, and rational-map point checks
have different scopes and are kept in separate output records.
"""
from __future__ import annotations
import json
from pathlib import Path
from math import comb,isqrt,factorial
from fractions import Fraction as Q

def need(v,m):
    if not v:raise ArithmeticError(m)
def primes(n):return [p for p in range(2,n+1) if all(p%d for d in range(2,isqrt(p)+1))]
def leg(a,p):
    v=pow(a%p,(p-1)//2,p)
    return -1 if v==p-1 else v
def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a or [0]
def add(a,b,p=None):
    o=[0]*max(len(a),len(b))
    for i,x in enumerate(a):o[i]+=x
    for i,x in enumerate(b):o[i]+=x
    return trim([x%p for x in o] if p else o)
def mul(a,b,p=None):
    o=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):o[i+j]+=x*y
    return trim([x%p for x in o] if p else o)
def scale(a,c,p=None):return trim([x*c%p if p else x*c for x in a])
def power(a,n,p=None):
    o=[1]
    while n:
        if n&1:o=mul(o,a,p)
        n//=2
        if n:a=mul(a,a,p)
    return o
def der(a,p):return trim([i*a[i]%p for i in range(1,len(a))])
def peval(a,x,p):
    o=0
    for c in a[::-1]:o=(o*x+c)%p
    return o
def invfrac(x,p):return x.numerator*pow(x.denominator,-1,p)%p
def pmul(a,b,p):return ((a[0]*b[0]+5*a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
def pvalue(poly,x,p):
    o=(0,0)
    for c in poly[::-1]:
        o=pmul(o,x,p);o=((o[0]+c)%p,o[1])
    return o
def legendre(n,p=None):
    # Rodrigues' explicit polynomial, not values fitted over a field.
    v=[Q(0)]*(n+1)
    for j in range(n//2+1):
        v[n-2*j]=Q((-1)**j*comb(n,j)*comb(2*n-2*j,n),2**n)
    return trim([invfrac(c,p) for c in v] if p else v)
def square_formula(n):
    v=[Q(0)]
    for j in range(n+1):
        v=add(v,scale(power([-1,0,1],j),Q(comb(n,j)*comb(n+j,j)*comb(2*j,j),4**j)))
    return v
def hasse(p):
    N=(p-1)//2;v=[0];iv2=pow(2,-1,p)
    aa=[-15*iv2%p,-9*iv2%p];bb=[7,9]
    for i in range(N+1):
        j=2*N-3*i;k=2*i-N
        if j>=0 and k>=0 and i+j+k==N:
            c=comb(N,i)*comb(N-i,j)%p
            v=add(v,scale(mul(power(aa,j,p),power(bb,k,p),p),c,p),p)
    return v
def hasse_translated(p):
    N=(p-1)//2;d=N//2;v=[0]
    for j in range(d+1):
        c=comb(N,2*j)*comb(2*j,j)*pow(6,N-2*j,p)*pow(9*pow(2,-1,p),j,p)%p
        v=add(v,scale(power([1,-1],j,p),c,p),p)
    return v

class Field:
    """F_p2=F_p[w]/(w^2-nu), represented by integer pairs."""
    def __init__(self,p):
        self.p=p
        self.nu=next(a for a in range(2,p) if leg(a,p)==-1)
        self.zero=(0,0);self.one=(1,0)
        self.elts=[(a,b) for a in range(p) for b in range(p)]
        self.sqr={}
        for a in self.elts:self.sqr.setdefault(self.mul(a,a),a)
    def con(self,a):
        if isinstance(a,Q):return (invfrac(a,self.p),0)
        return (a%self.p,0)
    def add(self,a,b):return ((a[0]+b[0])%self.p,(a[1]+b[1])%self.p)
    def neg(self,a):return ((-a[0])%self.p,(-a[1])%self.p)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        return ((a[0]*b[0]+self.nu*a[1]*b[1])%self.p,(a[0]*b[1]+a[1]*b[0])%self.p)
    def inv(self,a):
        d=pow((a[0]*a[0]-self.nu*a[1]*a[1])%self.p,-1,self.p)
        return (a[0]*d%self.p,-a[1]*d%self.p)
    def div(self,a,b):return self.mul(a,self.inv(b))
    def scale(self,a,c):return self.mul(a,self.con(c))
    def eval(self,poly,x):
        v=self.zero
        for c in poly[::-1]:v=self.add(self.mul(v,x),c)
        return v
    def diff(self,a):return [self.scale(a[i],i) for i in range(1,len(a))]
    def combine(self,a,b,root):return self.add(self.con(a),self.scale(root,b))

def curve_checks(p):
    F=Field(p);s0=F.sqr[F.con(5)];delta=F.sqr[F.con(-5)];rows=[]
    for s in [s0,F.neg(s0)]:
        C=lambda a,b=0:F.combine(a,b,s)
        aa=C(Q(-15,2),Q(-9,4));bb=C(7,Q(9,2))
        k=[C(Q(13,4),Q(21,10)),C(Q(-5,2),Q(-3,2)),C(1)]
        n=[C(Q(-1395,4),-135),C(Q(7045,16),Q(897,4)),
           C(Q(-1055,4),Q(-405,4)),C(60,Q(45,2)),C(-5,-3),C(1)]
        f=[bb,aa,F.zero,F.one]
        def on(P):
            return P is None or F.mul(P[1],P[1])==F.eval(f,P[0])
        def plus(P,Qp):
            if P is None:return Qp
            if Qp is None:return P
            x,y=P;xx,yy=Qp
            if x==xx:
                if y==F.neg(yy):return None
                m=F.div(F.add(F.scale(F.mul(x,x),3),aa),F.scale(y,2))
            else:m=F.div(F.sub(yy,y),F.sub(xx,x))
            nx=F.sub(F.sub(F.mul(m,m),x),xx)
            return nx,F.sub(F.mul(m,F.sub(x,nx)),y)
        def times(P,n):
            if n<0:
                if P is not None:P=(P[0],F.neg(P[1]))
                n=-n
            R=None
            while n:
                if n&1:R=plus(R,P)
                n//=2
                if n:P=plus(P,P)
            return R
        def phi(P):
            if P is None:return None
            x,y=P;kv=F.eval(k,x)
            if kv==F.zero:return None
            nv=F.eval(n,x);np=F.eval(F.diff(n),x);kp=F.eval(F.diff(k),x)
            rx=F.div(F.neg(nv),F.scale(F.mul(kv,kv),5))
            # r'=-(n'k-2nk')/(5k^3).
            dx=F.div(F.neg(F.sub(F.mul(np,kv),F.scale(F.mul(nv,kp),2))),
                     F.scale(F.mul(F.mul(kv,kv),kv),5))
            return rx,F.div(F.mul(y,dx),delta)
        points=[None]
        for x in F.elts:
            rhs=F.eval(f,x)
            if rhs in F.sqr:
                y=F.sqr[rhs];points.append((x,y))
                if y!=F.zero:points.append((x,F.neg(y)))
        kernel=0
        for P in points:
            need(on(P),"input curve")
            pp=phi(P);need(on(pp),"image outside curve")
            need(phi(pp)==times(P,-5),"phi^2 != [-5]")
            kernel+= pp is None
        homchecks=0
        for i in range(min(12,len(points))):
            P=points[i];Qp=points[-1-i]
            need(phi(plus(P,Qp))==plus(phi(P),phi(Qp)),"homomorphism sample")
            homchecks+=1
        rows.append({"p":p,"field_nonsquare":F.nu,"s":s,"delta":delta,
          "curve_points":len(points),"trace_over_p2":p*p+1-len(points),
          "map_compositions":len(points),"rational_kernel_points":kernel,
          "homomorphism_samples":homchecks})
    return rows

def run(output):
    polys=[]
    for p in primes(199):
        if p<=5:continue
        d=p//4;P=legendre(d,p);H=hasse(p);T=hasse_translated(p)
        need(H==T,"Hasse translation "+str(p))
        need(H==scale(P,leg(6,p),p),"Hasse/Legendre "+str(p))
        # The polynomial-square formula is also checked over every F_p here.
        sq=mul(P,P,p);rhs=[0]
        for j in range(d+1):
            cc=comb(d,j)*comb(d+j,j)*comb(2*j,j)*pow(4,-j,p)%p
            rhs=add(rhs,scale(power([-1,0,1],j,p),cc,p),p)
        need(sq==rhs,"square polynomial "+str(p))
        point=(0,pow(2,-1,p))
        value=pvalue(P,point,p);hp=pvalue(H,point,p);dp=pvalue(der(P,p),point,p)
        inert=leg(-5,p)==-1
        need((value==(0,0))==inert,"factor 4u^2-5 criterion")
        normdp=(dp[0]*dp[0]-5*dp[1]*dp[1])%p
        if inert:need(normdp!=0,"Hasse root not simple")
        polys.append({"p":p,"d":d,"hasse_coefficients":H,"legendre_coefficients":P,
              "hasse_at_s_over2":hp,"legendre_at_s_over2":value,
              "legendre_derivative_at_s_over2":dp,
              "inert":inert,"sqrt5_in_Fp":leg(5,p)==1})
    sqexact=[]
    for n in range(21):
        L=mul(legendre(n),legendre(n));R=square_formula(n)
        need(L==R,"exact Legendre square "+str(n))
        sqexact.append({"n":n,"coefficients":len(L),"ok":True})
    curves=[]
    for p in primes(43):
        if p>5:curves+=curve_checks(p)
    # Negative control: the finite transform requires t to be p-adically integral.
    p=7;c=7;t=Q(1,7);gate=1-4*c*t
    lhs=comb(12,6)*t**6;rhs=gate**3*comb(12,6)*(t/gate)**6
    need(lhs-rhs==Q(176,21609),"domain counterexample")
    # Original definitions, kept separate from theorem-level certificates.
    f=[sum(comb(n,k)**3 for k in range(n+1)) for n in range(199)]
    a=[sum(comb(n,k)**2*comb(2*k,k) for k in range(n+1)) for n in range(199)]
    qs=[sum(comb(n,k)*(-8)**(n-k)*f[k] for k in range(n+1)) for n in range(199)]
    reg=[];lifts=0
    for p in primes(199):
        if p in (2,5):continue
        m=p*p;av=[];bv=[]
        av=[comb(2*n,n)*a[n] for n in range(p)]
        bv=[comb(2*n,n)*qs[n] for n in range(p)]
        aa=peval(av,pow(20,-1,m),m);bb=peval(bv,pow(-16,-1,m),m)
        need(aa==bb,"two original values")
        eps=(-1)**((p-1)//2)
        if leg(-5,p)==-1:
            need(aa==bb==0,"inert original sums")
            for poly,t0 in [(av,pow(20,-1,m)),(bv,pow(-16,-1,m))]:
                need(peval(der(poly,p),t0,p)==0,"ordinary double root")
                need(peval(der(der(poly,p),p),t0,p)!=0,"ordinary triple root")
                if p<=97:
                    for r in range(p):
                        need(peval(poly,t0+r*p,m)==0,"lift invariance")
                        lifts+=1
        else:
            mult=1 if p%20 in (1,9) else 2
            reps=[(isqrt(mult*p-5*y*y),y) for y in range(1,isqrt(mult*p//5)+1)
                  if isqrt(mult*p-5*y*y)**2+5*y*y==mult*p]
            need(bool(reps),"missing representation")
            x,y=reps[0]
            target=((4 if mult==1 else 2)*x*x-2*p)%m
            need(eps*aa%m==target and eps*bb%m==target,"split values")
        reg.append({"p":p,"a":aa,"Q":bb,"eps":eps})
    summary={"ok":True,"polynomial_Hasse_rows":len(polys),
       "polynomial_Hasse_coefficients":sum(len(t['hasse_coefficients']) for t in polys),
       "exact_Legendre_square_rows":len(sqexact),
       "exact_Legendre_square_coefficients":sum(t["coefficients"] for t in sqexact),
       "finite_curves":len(curves),"finite_curve_primes":len(curves)//2,
       "finite_point_compositions":sum(t["map_compositions"] for t in curves),
       "finite_homomorphism_samples":sum(t["homomorphism_samples"] for t in curves),
       "fresh_original_prime_rows":len(reg),"ordinary_parameter_lifts":lifts,
       "missing_t_integrality":{"p":7,"c":7,"t":"1/7","weights":"delta_{j,6}",
                               "gate":"-3","lhs":str(lhs),"rhs":str(rhs),
                               "difference":"176/21609","v7_difference":-4}}
    for name,obj in [('finite_summary',summary),('Hasse_polynomials',polys),
          ('finite_curves',curves),('exact_Legendre_squares',sqexact),('original_sums',reg)]:
        (output/(name+'.json')).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    return summary
