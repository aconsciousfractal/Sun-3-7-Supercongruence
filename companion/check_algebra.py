"""Symbolic certificate replay. Requires SymPy; no floats or conjectural outputs."""
from __future__ import annotations
from pathlib import Path
import sympy as s
from .core import need,save


def run(output):
    records=[]
    def zero(name,expr):
        reduced=s.factor(s.cancel(s.expand(expr)))
        need(reduced==0,name+': '+str(reduced));records.append({'name':name,'residual':'0'})
    t,a,b,c=s.symbols('t a b c');y0,y1,y2,y3=s.symbols('Y Y1 Y2 Y3')
    D=1-4*a*t+16*c*t*t
    I=D*(2*y0*y2-y1*y1)+t*s.diff(D,t)*y0*y1+4*t*(-b+3*c*t)*y0*y0
    thI=t*s.diff(I,t)+y1*s.diff(I,y0)+y2*s.diff(I,y1)+y3*s.diff(I,y2)
    L=y3-2*t*(2*a*y3+3*a*y2+(a+2*b)*y1+b*y0)+4*c*t*t*(4*y3+12*y2+11*y1+3*y0)
    zero('theta I=2Y L',thI-2*y0*L)
    Yp,Ypp,Yppp=s.symbols('Yp Ypp Yppp')
    ordinary=s.expand(L.subs({y1:t*Yp,y2:t*t*Ypp+t*Yp,y3:t**3*Yppp+3*t*t*Ypp+t*Yp},simultaneous=True))
    zero('ordinary highest coefficient',ordinary.coeff(Yppp)-t**3*D)
    zero('ordinary second coefficient',ordinary.coeff(Ypp)-(3*t*t*D+s.Rational(3,2)*t**3*s.diff(D,t)))
    j,z,ell=s.symbols('j z ell')
    trans0=j*j*y0;trans1=j**4*y1+2*z*j**3*y0;trans2=j**6*y2+6*z*j**5*y1+6*z*z*j**4*y0
    zero('weight4 covariance',trans0*(j*j*ell+4*z*j)-2*trans1-j**4*(y0*ell-2*y1))
    zero('weight8 covariance',2*trans0*trans2-3*trans1**2-j**8*(2*y0*y2-3*y1*y1))
    bad=s.expand(2*trans0*trans2-4*trans1**2-j**8*(2*y0*y2-4*y1*y1))
    need(bad!=0,'wrong differential coefficient not detected')
    p=s.symbols('p',integer=True,positive=True);N=p-1
    top=16*c*(2*N*N-N*N)+32*c*N+12*c
    zero('top invariant coefficient',top-4*c*(2*p-1)*(2*p+1))
    h=s.symbols('h');P=h*(1+8*h)*(1+9*h);C=1+18*h+72*h*h;K=1-72*h*h
    R=(1+8*h)*(1+9*h)/h;S=h*(1+9*h)/(1+8*h)
    ta=P/C**2;tq=P/K**2;phi2=-(1+8*h)/(8*(1+9*h));phi6=1/(72*h)
    for label,phi in [('W2',phi2),('W6',phi6)]:
        zero('involution '+label,phi.subs(h,phi)-h)
        zero('t_a invariant '+label,ta.subs(h,phi)-ta)
        zero('t_Q invariant '+label,tq.subs(h,phi)-tq)
        zero('S product '+label,S.subs(h,phi)*S-s.Rational(1,64))
    zero('R W2 inverse',R.subs(h,phi2)*R-1)
    zero('R W6 fixed',R.subs(h,phi6)-R)
    zero('parameter Mobius relation',tq-ta/(1-36*ta))
    sq=s.sqrt(5);delta=s.I*sq;tau=(1+delta)/6
    zero('CM quadratic',6*tau*tau-2*tau+1)
    zero('inversion tau',-1/tau-(6*tau-2))
    zero('inversion 2tau',-1/(2*tau)-(3*tau-1))
    h0=(sq-2+s.I*(sq+2))/36
    zero('eta algebraic h',h0-(1-s.I)*tau*tau/2)
    def simplify_at(expr):return s.simplify(s.expand_complex(expr.subs(h,h0)))
    zero('t_a at tau',simplify_at(ta)-s.Rational(1,20))
    zero('t_Q at tau',simplify_at(tq)+s.Rational(1,16))
    v=h*(1+9*h)**3/(1+8*h);v0=-(9+4*sq)/64
    zero('v at tau',simplify_at(v)-v0)
    j0=632000+282880*sq
    zero('CM j from v',(1+256*v0)**3/v0-j0)
    u=s.symbols('u');EA=-s.Rational(3,2)*(3*u+5);EB=9*u+7
    disc=-16*(4*EA**3+27*EB**2);je=1728*4*EA**3/(4*EA**3+27*EB**2)
    zero('elliptic j equals CM j',je.subs(u,sq/2)-j0)
    zero('elliptic discriminant',disc.subs(u,sq/2)-729*(sq-2))
    zero('discriminant norm',729*(sq-2)*729*(-sq-2)+729**2)
    zero('a target D',D.subs({a:10,b:3,c:9,t:s.Rational(1,20)})+s.Rational(16,25))
    zero('Q target D',D.subs({a:-17,b:-6,c:72,t:-s.Rational(1,16)})-s.Rational(5,4))
    zero('f target D',D.subs({a:7,b:2,c:-8,t:s.Rational(1,16)})+s.Rational(5,4))
    xx=(2+sq)/256
    zero('quartic Clausen parameter',xx*(1-64*xx)+s.Rational(1,1024))
    # Exact generic Hecke conjugations. Integrality follows from these integer polynomials.
    A,B,Cc,E=s.symbols('A B C E');M=s.Matrix([[A,B],[6*Cc,E]])
    W2=s.Matrix([[2,-1],[6,-2]]);W6=s.Matrix([[0,-1],[6,0]])
    conjugations={}
    for name,W in [('W2',W2),('W6',W6)]:
        mm=W.inv()*M*W
        need(all(s.Poly(v,A,B,Cc,E).domain.is_ZZ for v in mm),'integer Hecke conjugation '+name)
        need(s.Poly(mm[1,0]/6,A,B,Cc,E).domain.is_ZZ,'level Hecke conjugation '+name)
        zero('conjugated determinant '+name,mm.det()-M.det())
        conjugations[name]=[[str(mm[i,k]) for k in range(2)] for i in range(2)]
    for den in [C,K]:
        need(s.gcd(P,den)==1 and s.discriminant(den,h)!=0,'parameter no cancellation / simple zeros')
    # j identity after the E4 identity: E4=eta1^16/eta2^8*(1+256v).
    u1,u2=s.symbols('eta1 eta2',nonzero=True)
    zero('j eta formula algebra',(u1**16/u2**8+256*u2**16/u1**8)**3/u1**24-(1+256*u2**24/u1**24)**3/(u2**24/u1**24))
    out={'ok':True,'identities':records,'identity_count':len(records),'wrong_covariance_residual':str(bad),'hecke_conjugations':conjugations,'CM_exact':{'tau':str(tau),'h':str(h0),'v':str(v0),'j':str(j0),'curve_discriminant':str(729*(sq-2))},'no_float':True}
    save(output/'algebra.json',out)
    return out

