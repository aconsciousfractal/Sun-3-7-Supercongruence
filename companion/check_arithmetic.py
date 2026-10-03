"""Recurrence-based prime regressions, original sequence checks, CM matrices,
Hasse bridges, and exact finite invariants.
All modular divisions are checked inverses; no floating point; no assertions.
"""
from __future__ import annotations
import csv
from fractions import Fraction as F
from math import comb,isqrt,factorial,gcd
from pathlib import Path
from .core import *


def bezout(a,b):
    aa,bb=a,b;u,v,s,t=1,0,0,1
    while bb:
        q=aa//bb;aa,bb=bb,aa-q*bb;u,s=s,u-q*s;v,t=t,v-q*t
    if aa<0:aa,u,v=-aa,-u,-v
    need(aa==1,'Bezout not primitive');need(a*u+b*v==1,'Bezout arithmetic')
    return u,v


def qadd(x,y):return x[0]+y[0],x[1]+y[1]
def qmul(x,y,d):return x[0]*y[0]+d*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def qinv(x,d):
    norm=x[0]**2-d*x[1]**2;need(norm!=0,'quadratic denominator zero');return F(x[0],norm),F(-x[1],norm)
def qdiv(x,y,d):return qmul(x,qinv(y,d),d)
def qscale(x,c):return x[0]*c,x[1]*c


def represent(p):
    if p%20 in (11,13,17,19):return 0,None
    m=1 if p%20 in (1,9) else 2
    for y in range(1,isqrt(m*p//5)+1):
        x=isqrt(m*p-5*y*y)
        if x*x+5*y*y==m*p:return m,(x,y)
    raise ArithmeticError('No expected representation for '+str(p))


def cm(p,x,y,m):
    tau=(F(1,6),F(1,6));d=-5
    if m==1:
        a,b=bezout(x+y,6*y);M=[[a,b],[-6*y,x+y]];k=b*(x-y)-a*y;character=1
    else:
        a,b=bezout(x+y,3*y);M=[[2*a,b],[-6*y,x+y]];k=b*(x-y)//2-a*y;character=-1
        need(b*(x-y)%2==0,'nonintegral CM k')
        coset=[[-2*a-3*b,a+b],[3*(y-x),x-2*y]]
        need(coset[0][0]*coset[1][1]-coset[0][1]*coset[1][0]==1 and coset[1][0]%6==0,'W2 coset')
    det=M[0][0]*M[1][1]-M[0][1]*M[1][0];need(det==m and M[1][0]%6==0,'CM determinant')
    denom=qadd(qscale(tau,M[1][0]),(M[1][1],0))
    image=qdiv(qadd(qscale(tau,M[0][0]),(M[0][1],0)),denom,d)
    need(image==qscale(qadd(tau,(k,0)),F(1,p)),'CM fractional identity')
    root=qscale(qinv(qmul(denom,denom,d),d),F(det*p*p,character))
    expected=qscale(qmul((x,y),(x,y),d),F(character,m))
    need(root==expected,'exact root from character and determinant')
    norm=root[0]**2+5*root[1]**2;trace=2*root[0]
    need(norm==p*p,'CM norm');need(trace.denominator==1,'rational integral trace')
    # Construct one actual unit embedding modulo p^2 by a one-step Hensel lift.
    s=next(z for z in range(p) if (z*z+5)%p==0)
    s=(s+p*(-((s*s+5)//p)*pow(2*s,-1,p)%p))%(p*p)
    need((s*s+5)%(p*p)==0,'Hensel lift')
    val=residue(root[0],p*p)+residue(root[1],p*p)*s
    if val%p==0:s=(-s)%(p*p);val=residue(root[0],p*p)+residue(root[1],p*p)*s
    need(val%p!=0,'unit root embedding');need((val-int(trace))%(p*p)==0,'unit root versus trace')
    return {'p':p,'multiplier':m,'x':x,'y':y,'M':M,'k':k,'character':character,
            'root_coefficients':list(map(str,root)),'trace':int(trace),'norm':int(norm),'sqrt_minus5_mod_p2':s,'unit_root_mod_p2':val%(p*p)}


def invpoly(U,a,b,c):
    n=2*len(U);DU=theta(U);DDU=theta(DU)
    D=[1,-4*a,16*c]
    return add(mul(D,add(scale(mul(U,DDU,n),2),scale(mul(DU,DU,n),-1)),n),
               mul(theta(D),mul(U,DU,n),n),mul([0,-4*b,12*c],mul(U,U,n),n))


def hasse(p):
    N=(p-1)//2;out=[0]*(p//4+1)
    A=[-15*pow(2,-1,p)%p,-9*pow(2,-1,p)%p];B=[7,9]
    for i in range((N+1)//2,2*N//3+1):
        j=2*N-3*i;k=2*i-N
        need(i+j+k==N and 3*i+j==2*N and min(i,j,k)>=0,'Hasse support')
        coefficient=comb(N,i)*comb(N-i,j)
        term=modmul(modpow(A,j,p),modpow(B,k,p),p)
        need(len(term)<=len(out),'Hasse degree bound')
        for l,v in enumerate(term):out[l]=(out[l]+coefficient*v)%p
    return out


def legendre_coefficients(n,p):
    if n==0:return [1]
    prev,current=[1],[0,1]
    for k in range(1,n):
        following=add(scale([0]+current,2*k+1),scale(prev,-k))
        following=[v*pow(k+1,-1,p)%p for v in following]
        prev,current=current,following
    return current


def qmodmul(x,y,m,d=5):return ((x[0]*y[0]+d*x[1]*y[1])%m,(x[0]*y[1]+x[1]*y[0])%m)
def qmodev(poly,x,m):
    out=(0,0)
    for a in reversed(poly):
        out=qmodmul(out,x,m);out=((out[0]+a)%m,out[1])
    return out


def finite_bridges(p,g):
    N=(p-1)//2;eps=(-1)**N
    f16=sum_mod(g['f'],p,F(1,16),1)
    aa=sum_mod(g['a'],p,F(1,20),1);qq=sum_mod(g['Q'],p,F(-1,16),1)
    need((aa-leg(5,p)*f16)%p==0 and (qq-eps*f16)%p==0,'finite transform modulo p')
    # Exact inner Vandermonde identity used in the mod-p Jacobi interchange.
    for j in range(N+1):
        lhs=sum(comb(N,k)**2*comb(k,j)*comb(N-k,j) for k in range(j,N-j+1))
        rhs=comb(N,j)**2*comb(2*N-2*j,N-2*j) if 2*j<=N else 0
        need(lhs==rhs,'exact Vandermonde inner identity')
    jacobi=0
    for k in range(N+1):
        term=total=1
        for j in range(k):
            term=term*(-k+j)*(2*k+2*j+1)*pow(-8*(j+1)**2,-1,p)%p
            total=(total+term)%p
        jacobi=(jacobi+comb(2*k,k)**2*pow(16,-k,p)*total)%p
    need(jacobi==f16,'Dixon bucket Jacobi identity modp')
    rv=sum(comb(2*k,k)**2*comb(4*k,2*k)*pow(-1024,-k,p) for k in range(p))%p
    need((jacobi-eps*rv)%p==0,'Jacobi RV modp bridge')
    H=hasse(p);PL=legendre_coefficients(p//4,p)
    need(PL==[leg(6,p)*x%p for x in H],'polynomial Hasse/Legendre identity')
    rhalf=(0,pow(2,-1,p));hval=qmodev(H,rhalf,p);pval=qmodev(PL,rhalf,p)
    square=qmodmul(pval,pval,p)
    need(square==(rv,0),'RV Legendre square over quadratic algebra')
    # Full polynomial Clausen specialization modulo p^2 in Z/p^2[r]/(r^2-5).
    M=p*p; xx=(2*pow(256,-1,M)%M,pow(256,-1,M));power=(1,0);cs=(0,0)
    for k in range(p):
        coeff=comb(2*k,k)*comb(4*k,2*k)%M
        cs=((cs[0]+coeff*power[0])%M,(cs[1]+coeff*power[1])%M)
        power=qmodmul(power,xx,M)
    rvp2=sum(comb(2*k,k)**2*comb(4*k,2*k)*pow(-1024,-k,M) for k in range(p))%M
    need(qmodmul(cs,cs,M)==(rvp2,0),'finite Clausen at algebraic parameter modp2')
    if leg(-5,p)==-1:
        need(hval==(0,0),'CM Hasse vanishing inert')
        need(aa==qq==f16==rv==0,'all modp inert values')
    return {'p':p,'Hasse_coefficients':H,'Legendre_coefficients':PL,'Hasse_at_sqrt5_over2':hval,
            'Legendre_at_sqrt5_over2':pval,'F16_modp':f16,'a20_modp':aa,'Qminus16_modp':qq,
            'RV_modp':rv,'RV_modp2':rvp2,'Clausen_sum_modp2':cs,'inert':leg(-5,p)==-1,
            'Vandermonde_rows':N+1}


def run(out, bound=1999):
    need(type(bound) is int and bound>=199,'bound must be an integer including controls through 199')
    u={name:sequence(*par,bound) for name,par in PARAMS.items()};g={name:central(v) for name,v in u.items()}
    direct=[]
    for n in range(199):
        f=sum(comb(n,k)**3 for k in range(n+1))
        a=sum(comb(n,k)**2*comb(2*k,k) for k in range(n+1))
        q=sum(comb(n,k)*(-8)**(n-k)*u['f'][k] for k in range(n+1))
        need((f,a,q)==(u['f'][n],u['a'][n],u['Q'][n]),'original sequence definitions')
        need(a==sum(comb(n,k)*u['f'][k] for k in range(n+1)),'exact binomial a transform')
        need(f==sum(comb(n,k)**2*comb(2*k,n) for k in range(n+1) if 2*k>=n),'exact Dixon identity')
        direct.append({'n':n,'f':f,'a':a,'Q':q})
    save(out/'original_sequences.json',direct)
    rows=[];matrices=[]
    for p in primes(bound):
        if p in (2,5):continue
        m,rep=represent(p);eps=(-1)**((p-1)//2)
        target=(4*rep[0]**2-2*p if m==1 else 2*rep[0]**2-2*p)%(p*p) if rep else 0
        A=sum_mod(g['a'],p,F(1,20));Q=sum_mod(g['Q'],p,F(-1,16));f16=sum_mod(g['f'],p,F(1,16))
        need(eps*A%(p*p)==target and eps*Q%(p*p)==target,'Sun 3.7 exact finite regression p='+str(p))
        need(g['a'][p]%(p*p)==6%(p*p) and g['Q'][p]%(p*p)==(-12)%(p*p),'endpoint coefficients')
        if p>5 and rep:matrices.append(cm(p,*rep,m))
        rows.append({'p':p,'branch':'inert' if m==0 else 'principal' if m==1 else 'nonprincipal','x':rep[0] if rep else '', 'y':rep[1] if rep else '', 'epsilon':eps,'a20':A,'Qminus16':Q,'F16':f16,'signed_target':target})
    with (out/'prime_regression.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\r\n');w.writeheader();w.writerows(rows)
    save(out/'CM_matrices.json',matrices)
    bridges=[];polys=[];numcoef=0
    for p in primes(199):
        if p<=5:continue
        bridges.append(finite_bridges(p,g))
        for name,t0 in [('a',F(1,20)),('Q',F(-1,16)),('f',F(1,16))]:
            a,b,c=PARAMS[name];M=p*p
            U0=[v%M for v in g[name][:p]]
            I=invpoly(U0,a,b,c)
            need(not any(v%M for v in I),'finite invariant '+name+str(p))
            # Compute exact integer resonant coefficients independently from the full unreduced I.
            exact=invpoly(g[name][:p],a,b,c)
            need(exact[p]==-2*p*p*g[name][p],'resonance degree p exact')
            need(exact[2*p]==4*c*(2*p-1)*(2*p+1)*g[name][p-1]**2,'resonance degree 2p exact')
            D=residue(1-4*a*t0+16*c*t0*t0,p)
            need(D!=0,'target non-singular for p>5')
            val=sum_mod(g[name],p,t0)
            if val%p==0:
                first=sum_mod(theta(g[name]),p,t0,1);second=sum_mod(theta(theta(g[name])),p,t0,1)
                need(first==0 and second!=0 and val==0,'ordinary zero lift conditions')
            numcoef+=len(I)
            polys.append({'p':p,'form':name,'number_of_coefficients':len(I),'resonance_p_verified':True,'resonance_2p_verified':True,'discriminant_modp':D})
    save(out/'finite_Hasse_bridges.json',bridges);save(out/'finite_invariants.json',polys)
    negative={}
    # Deliberate corruption of the source, never used as a theorem counterexample.
    control=next(r for r in matrices if r['p']==7)
    negative['p7_CM']={'correct_raw':control['unit_root_mod_p2'],'wrong_character':(-control['unit_root_mod_p2'])%49,'omit_half':(2*control['unit_root_mod_p2'])%49}
    need(len(set(negative['p7_CM'].values()))==3,'CM negative controls discriminate')
    deltas=[]
    for p in (7,11,13):
        raw=sum(comb(2*n,n)*pow(20,-n,p*p) for n in range(p))%(p*p)
        target=leg(5,p)%(p*p)
        if raw!=target:deltas.append({'p':p,'a20_delta_weight':raw,'Legendre_delta_comparison':target})
    need(deltas,'generic modp2 binomial lift control');negative['false_generic_lift']=deltas
    cubics=[]
    for p in (7,11,13,19,29):
        A=sum_mod(g['a'],p,F(1,20),3);Q=sum_mod(g['Q'],p,F(-1,16),3)
        if A!=Q:cubics.append({'p':p,'a20_modp3':A,'Qminus16_modp3':Q})
    need(cubics,'no blanket modp3 equality control');negative['modp3_nonequality']=cubics
    negative['p3']={k:sum_mod(g[k],3,F(1,20) if k=='a' else F(-1,16)) for k in ('a','Q')}
    negative['mutated_a20_plus1_rejected']=all((r['epsilon']*(r['a20']+1)-r['signed_target'])%(r['p']**2)!=0 for r in rows)
    need(negative['mutated_a20_plus1_rejected'],'corrupted sum rejected')
    save(out/'negative_controls.json',negative)
    summary={'ok':True,'max_prime':rows[-1]['p'],'primes':len(rows),'branch_counts':{b:sum(r['branch']==b for r in rows) for b in ('principal','nonprincipal','inert')},'direct_sequence_values':3*len(direct),'CM_matrices':len(matrices),'finite_invariants':len(polys),'invariant_coefficients':numcoef,'exact_resonant_coefficients':2*len(polys),'Hasse_polynomial_identities':len(bridges),'Vandermonde_integer_rows':sum(r['Vandermonde_rows'] for r in bridges)}
    save(out/'arithmetic_summary.json',summary)
    return summary

