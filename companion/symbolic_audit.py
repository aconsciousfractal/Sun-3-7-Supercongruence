"""Second exact engine for four isogeny identities and three auxiliary identities.
SymPy is used here; neither the rational-pair nor the finite-field engine is imported.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as s

def need(v,m):
    if not v:raise ArithmeticError(m)

def run(output):
    x=s.symbols('x');a=-(30+9*s.sqrt(5))/4;b=(14+9*s.sqrt(5))/2
    K=s.QQ.algebraic_field(s.sqrt(5));p=lambda e:s.Poly(e,x,domain=K)
    k=p(x*x-(5+3*s.sqrt(5))*x/2+(65+42*s.sqrt(5))/20)
    n=p(x**5-(5+3*s.sqrt(5))*x**4+(120+45*s.sqrt(5))*x**3/2
        -(1055+405*s.sqrt(5))*x*x/4+(7045+3588*s.sqrt(5))*x/16
        -(1395+540*s.sqrt(5))/4)
    A=-n;B=5*k*k;f=p(x**3+a*x+b);D=A.diff()*B-A*B.diff()
    rows=[]
    def zero(name,expr):
        ok=expr.is_zero if isinstance(expr,s.Poly) else s.cancel(s.expand(expr))==0
        need(ok,name);rows.append({"name":name,"residual":"0"})
    zero('curve',f*D**2+5*(A**3*B+p(a)*A*B**3+p(b)*B**4))
    F4=p(x**6+5*a*x**4+20*b*x**3-5*a*a*x*x-4*a*b*x-a**3-8*b*b)
    p3=p(3*x**4+6*a*x*x+12*b*x-a*a)
    p5=32*f*f*F4-p3**3
    mnum=p(x)*p5*p5-8*f*F4*p3*(p5-4*F4*F4)
    def hom(P):
        return sum((p(c)*A**e*B**(5-e) for (e,),c in P.terms()),p(0))
    zero('composition numerator',hom(A)-mnum)
    zero('composition denominator',hom(B)-p5*p5)
    p7=p5*p3**3-128*f*f*F4**3
    om5=4*p7*F4**2-p3**3*(p5-4*F4*F4)**2
    zero('y-coordinate normalization of multiplication by 5',
         mnum.diff()*p5-2*mnum*p5.diff()-5*om5)
    res=s.simplify(s.resultant(n.as_expr(),k.as_expr(),x))
    need(s.simplify(res-3**10*(25-11*s.sqrt(5))/250)==0,'resultant')
    # Polynomial square identity derived from the Legendre differential equation.
    u,l=s.symbols('u lambda');d=1-u*u;w=s.Function('w')(u)
    T=lambda y:-d*s.diff(y,u)/(2*u)
    lhs=u*u*T(T(T(w)))-s.Rational(3,2)*d*T(T(w))-d*(s.Rational(1,2)-l)*T(w)+d*l*w/2
    sym=d*d*s.diff(w,u,3)-6*u*d*s.diff(w,u,2)+(6*u*u-2+4*l*d)*s.diff(w,u)-4*l*u*w
    zero('Clausen operator change z=1-u^2',lhs+d*sym/(8*u))
    y0,y1,y2,y3=s.symbols('y0 y1 y2 y3')
    jets={w:y0*y0,s.diff(w,u):2*y0*y1,
          s.diff(w,u,2):2*y1*y1+2*y0*y2,
          s.diff(w,u,3):6*y1*y2+2*y0*y3}
    symj=sym.xreplace(jets)
    y2v=(2*u*y1-l*y0)/d
    y3v=s.diff(y2v,u)+s.diff(y2v,y0)*y1+s.diff(y2v,y1)*y2v
    zero('symmetric square of Legendre',symj.subs({y3:y3v,y2:y2v},simultaneous=True))
    # Exact cubic translation used for Hasse, not a formula inferred on F_p points.
    z,v=s.symbols('z v')
    cubic=x**3-s.Rational(3,2)*(3*v+5)*x+9*v+7
    zero('Hasse cubic translation',cubic.subs(x,z+2)-z*(z*z+6*z+s.Rational(9,2)*(1-v)))
    # Differential multipliers alone do not identify
    # maps in characteristic p. The explicit y-coordinate identity above fixes the sign.
    out={"ok":True,"engine":"SymPy exact algebraic polynomial domain",
         "identities":rows,"resultant":str(res)}
    (output/'symbolic_summary.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n', encoding='utf-8', newline='\n')
    return out
