"""Exact algebraic certificate over Q[s]/(s^2-5).
This standard-library engine uses rational pairs and polynomial convolutions.
Coefficients are rational pairs; polynomials are ascending coefficient lists.
"""
from __future__ import annotations
import json
from fractions import Fraction as Q
from pathlib import Path

R = tuple[Q,Q]
Z: R = (Q(0),Q(0))
O: R = (Q(1),Q(0))
S: R = (Q(0),Q(1))

def r(a=0,b=0)->R: return (Q(a),Q(b))
def ra(x:R,y:R)->R:return x[0]+y[0],x[1]+y[1]
def rn(x:R)->R:return -x[0],-x[1]
def rm(x:R,y:R)->R:return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def ri(x:R)->R:
    n=x[0]*x[0]-5*x[1]*x[1]
    if not n: raise ZeroDivisionError("zero algebraic norm")
    return x[0]/n,-x[1]/n
def rs(x:R,c)->R:return rm(x,r(c))
def norm(x:R)->Q:return x[0]*x[0]-5*x[1]*x[1]

def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==Z:a.pop()
    return a or [Z]
def add(*aa):
    out=[Z]*max(map(len,aa))
    for a in aa:
        for i,c in enumerate(a):out[i]=ra(out[i],c)
    return trim(out)
def neg(a):return [rn(c) for c in a]
def scale(a,c):return trim([rm(x,c) for x in a])
def mul(a,b):
    out=[Z]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x!=Z:
            for j,y in enumerate(b):
                if y!=Z:out[i+j]=ra(out[i+j],rm(x,y))
    return trim(out)
def power(a,n):
    if n<0:raise ValueError("negative polynomial exponent")
    out=[O]
    while n:
        if n&1:out=mul(out,a)
        n//=2
        if n:a=mul(a,a)
    return out
def der(a):return trim([rs(a[i],i) for i in range(1,len(a))])
def divmodp(a,b):
    a=trim(a);b=trim(b)
    if b==[Z]:raise ZeroDivisionError("zero polynomial")
    q=[Z]*max(1,len(a)-len(b)+1)
    while a!=[Z] and len(a)>=len(b):
        n=len(a)-len(b);c=rm(a[-1],ri(b[-1]));q[n]=ra(q[n],c)
        a=add(a,[Z]*n+neg(scale(b,c)))
    return trim(q),trim(a)
def homogeneous(poly,A,B,n):
    out=[Z]
    for i,c in enumerate(poly):
        out=add(out,scale(mul(power(A,i),power(B,n-i)),c))
    return out
def enc(a):return [[str(x),str(y)] for x,y in a]
def need(test,message):
    if not test:raise ArithmeticError(message)

def data():
    aa=r(Q(-15,2),Q(-9,4));bb=r(7,Q(9,2))
    k=[r(Q(13,4),Q(21,10)),r(Q(-5,2),Q(-3,2)),O]
    n=[r(Q(-1395,4),-135),r(Q(7045,16),Q(897,4)),
       r(Q(-1055,4),Q(-405,4)),r(60,Q(45,2)),r(-5,-3),O]
    f=[bb,aa,Z,O]
    return aa,bb,k,n,f

def run(output):
    aa,bb,k,n,f=data();X=[Z,O];A=neg(n);B=scale(power(k,2),r(5))
    records=[]
    def zero(name,a,b,slots=None):
        out=add(a,neg(b));need(out==[Z],name+": "+str(enc(out)))
        records.append({"name":name,"coefficient_pairs":slots or max(len(a),len(b)),
                        "residual":"0"})
    # phi=(A/B, y (A/B)'/delta), delta^2=-5.
    H=add(mul(der(A),B),neg(mul(A,der(B))))
    rhs=add(mul(power(A,3),B),scale(mul(A,power(B,3)),aa),scale(power(B,4),bb))
    zero("curve map: f H^2+5(A^3 B+a A B^3+b B^4)",mul(f,power(H,2)),
         scale(rhs,r(-5)))
    # Standard division polynomial construction for multiplication by 5.
    P3=add(scale(power(X,4),r(3)),scale(power(X,2),rs(aa,6)),
           scale(X,rs(bb,12)),[rn(rm(aa,aa))])
    F4=add(power(X,6),scale(power(X,4),rs(aa,5)),
           scale(power(X,3),rs(bb,20)),scale(power(X,2),rs(rm(aa,aa),-5)),
           scale(X,rs(rm(aa,bb),-4)),
           [rn(ra(rm(rm(aa,aa),aa),rs(rm(bb,bb),8)))])
    P5=add(scale(mul(power(f,2),F4),r(32)),neg(power(P3,3)))
    den5=power(P5,2)
    num5=add(mul(X,den5),
             scale(mul(mul(mul(f,F4),P3),add(P5,scale(power(F4,2),r(-4)))),r(-8)))
    AC=homogeneous(A,A,B,5);BC=homogeneous(B,A,B,5)
    zero("self-composition numerator equals numerator x([5])",AC,num5)
    zero("self-composition denominator equals psi5^2",BC,den5)
    # Exact y-coordinate normalization for [5]; do not infer equality of maps
    # in characteristic p from their differentials alone.
    P7=add(mul(P5,power(P3,3)),scale(mul(power(f,2),power(F4,3)),r(-128)))
    Om5=add(scale(mul(P7,power(F4,2)),r(4)),
            neg(mul(power(P3,3),power(add(P5,scale(power(F4,2),r(-4))),2))))
    zero("derivative of x([5]) equals 5 omega5/y",
         add(mul(der(num5),P5),scale(mul(num5,der(P5)),r(-2))),scale(Om5,r(5)))
    # Exact nondegeneration. Since k is monic, res(k,n) is the norm of n mod k.
    rem=divmodp(n,k)[1];need(len(rem)==2,"linear remainder")
    m,l=rem
    resultant=ra(ra(rm(rm(l,l),k[0]),rn(rm(rm(l,m),k[1]))),rm(m,m))
    wanted=r(Q(3**10,10),Q(-11*3**10,250))
    need(resultant==wanted,"resultant")
    need(norm(resultant)==Q(3**20,5**5),"norm resultant")
    dk=ra(rm(k[1],k[1]),rs(k[0],-4))
    need(dk==r(Q(9,2),Q(-9,10)) and norm(dk)==Q(81,5),"kernel discriminant")
    delta=rs(ra(rs(rm(rm(aa,aa),aa),4),rs(rm(bb,bb),27)),-16)
    need(delta==r(-1458,729) and norm(delta)==-(729**2),"curve discriminant")
    # All coefficients that will be specialized are integral away from 2 and 5.
    for poly in [A,B,k,P3,F4,P5,num5,den5]:
        for pair in poly:
            for c in pair:
                d=c.denominator
                for q in [2,5]:
                    while d%q==0:d//=q
                need(d==1,"unaccounted prime in denominator")
    # Deliberate faults: wrong scaling, wrong curve sign, altered input.
    mutation=[]
    for label,am,bm,sign in [
        ("omit divisor 5",A,power(k,2),5),
        ("replace delta^2=-5 by +5",A,B,-5),
        ("change numerator constant by 1",add(A,[O]),B,5)]:
        hh=add(mul(der(am),bm),neg(mul(am,der(bm))))
        rr=add(mul(power(am,3),bm),scale(mul(am,power(bm,3)),aa),scale(power(bm,4),bb))
        out=add(mul(f,power(hh,2)),scale(rr,r(sign)))
        need(out!=[Z],"mutation not rejected: "+label)
        first=next(i for i,c in enumerate(out) if c!=Z)
        mutation.append({"name":label,"first_nonzero_degree":first,
                         "coefficient":enc([out[first]])[0]})
    certificate={"field":"Q[s]/(s^2-5)","order":"ascending",
      "a":enc([aa])[0],"b":enc([bb])[0],"k":enc(k),"n":enc(n),
      "A":enc(A),"B":enc(B),"D=A'B-AB'":enc(H),"psi3":enc(P3),
      "F4=psi4/(4y)":enc(F4),"psi5":enc(P5),"multiplication5_numerator":enc(num5),
      "multiplication5_denominator":enc(den5),"self_composition_numerator":enc(AC),
      "self_composition_denominator":enc(BC),"omega5_over_y":enc(Om5)}
    summary={"ok":True,"engine":"standard-library rational-pair polynomial arithmetic",
      "identities":records,"residual_coefficient_pairs":sum(x["coefficient_pairs"] for x in records),
      "resultant_n_k":enc([resultant])[0],"resultant_norm":str(norm(resultant)),
      "kernel_discriminant":enc([dk])[0],"kernel_discriminant_norm":str(norm(dk)),
      "curve_discriminant":enc([delta])[0],"curve_discriminant_norm":str(norm(delta)),
      "denominator_primes":[2,5],"nondegenerate_characteristics":"p>5",
      "mutation_controls":mutation}
    for name,obj in [('isogeny_certificate',certificate),('isogeny_summary',summary)]:
        (output/(name+'.json')).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    return summary

