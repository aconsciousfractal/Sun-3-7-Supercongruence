"""Exact finite modular certificates and bounded Hecke regressions.
The text proves modularity and valence bounds; finite coefficients alone are not a proof.
Polynomial reconstruction never uses CM target values to choose coefficients.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import comb
from pathlib import Path
from .core import *


def certificates(out):
    n=24;one=[1]+[0]*n
    d=modular_data(n,'product');e=modular_data(n,'log')
    need(d==e,'independent eta expansions mismatch')
    h=d['h'];hh=mul(h,h,n);j=d['J']
    forms=d['forms'];ts=d['parameters']
    mr={name:eta_metadata(r) for name,r in [('h',H),('U',U),('B',B),('V',V),('J',J),('R',R),('v',H2),('S',S)]}
    for name,row in mr.items():
        need(row['sum_dr']%24==row['sum_Nd_r']%24==0 and row['trivial_character'], 'Newman '+name)
    # qR is computed from its eta product, not from h.
    qR=eta(R,n,'product',remove_shift=True)
    need(qR==eta(R,n,'log',remove_shift=True),'two qR engines')
    Rh_rhs=add(one,scale(h,17),scale(hh,72))
    rres=add(mul(qR,h,n),scale([0]+Rh_rhs[:-1],-1))
    need(not any(rres),'R rational identity')
    v=eta(H2,n,'product'); ss=eta(S,n,'product')
    vres=add(mul(v,add(one,scale(h,8)),n),scale(mul(h,mpow(add(one,scale(h,9)),3,n),n),-1))
    sres=add(mul(ss,add(one,scale(h,8)),n),scale(mul(h,add(one,scale(h,9)),n),-1))
    need(not any(vres+sres),'v or S identity')
    logR=log_eta(R,n);logS=log_eta(S,n);logJ=log_eta(J,n)
    w2={'a':add(forms['a'],scale(logS,-1)),'Q':add(forms['Q'],logR)}
    need(not any(w2['a']+w2['Q']),'weight2 log identification')
    residues={};mutations=[]
    P=mul(h,Rh_rhs,n)
    for name,th in forms.items():
        a,b,c=PARAMS[name];W=add(mul(th,logJ,n),scale(theta(th),-2));th2=mpow(th,2,n)
        p1=add(mul(W,W,n),scale(mpow(th,4,n),-1),scale(mul(j,th2,n),4*a),scale(mul(j,j,n),-16*c))
        p2=add(scale(mul(th,theta(theta(th)),n),2),scale(mul(theta(th),theta(th),n),-3),scale(mul(j,th2,n),-4*b),scale(mul(j,j,n),12*c))
        need(not any(p1+p2),'weight8 '+name)
        poly=add(one,scale(h,18),scale(hh,72)) if name=='a' else add(one,scale(hh,-72))
        need(mul(ts[name],mpow(poly,2,n),n)==P,'rational parameter '+name)
        need(compose(central(sequence(a,b,c,n)),ts[name],n)==th,'period recurrence '+name)
        residues[name]={'first':p1,'second':p2,'required_through':8,'required_coefficients':18}
        for label,bad in [('a+1',add(p1,scale(mul(j,th2,n),4))),('b+1',add(p2,scale(mul(j,th2,n),-4))),('c+1',add(p1,scale(mul(j,j,n),-16)))]:
            need(any(bad[:9]),'parameter mutation not detected')
            mutations.append({'form':name,'mutation':label,'first_nonzero':next(i for i,x in enumerate(bad) if x),'value':next(x for x in bad if x)})
    # Classical j formula: both summands below are holomorphic weight 4 on Gamma_0(2).
    e4=[1]+[240*sum(d**3 for d in range(1,k+1) if k%d==0) for k in range(1,n+1)]
    jres=add(e4,scale(eta({1:16,2:-8},n),-1),scale(eta({1:-8,2:16},n),-256))
    need(not any(jres),'E4 identity')
    save(out/'modular_certificates.json',{'eta_metadata':mr,'series':d,'qR':qR,'v':v,'S':ss,
        'R_residual_after_multiplication_by_q':rres,'v_residual':vres,'S_residual':sres,
        'weight2':w2,'weight8':residues,'E4_residual':jres,'parameter_mutations':mutations,
        'required_coefficient_counts':{'weight8':36,'weight2':6,'E4':2,'v_identity':6},
        'note':'Exact finite arithmetic; the manuscript supplies the uniform closure, divisor, covariance and cusp-regularity arguments.'})
    return {'forms':2,'weight8_identities':4,'weight8_required_coefficients':36,'weight2_required_coefficients':6,'E4_required_coefficients':2,'series_through':n,'mutations_detected':len(mutations)}


def exact_div(a,k):
    need(all(x%k==0 for x in a),'Newton division '+str(k));return [x//k for x in a]


def quadratic_remainder(poly,trace,norm):
    a=list(map(F,poly))
    for k in range(len(a)-1,1,-1):
        v=a[k];a[k]=0;a[k-1]+=v*trace;a[k-2]-=v*norm
    return a[:2]


def hecke(p,name):
    n=p+3; degree_s=p*n
    ths=modular_data(degree_s,'log')['forms'][name]
    th=ths[:n+1];t=modular_data(n,'log')['parameters'][name]
    traces=[[0]*(n+1)];current=[1]+[0]*degree_s
    for k in range(1,p+1):
        current=mul(current,ths,degree_s)
        traces.append([p*current[p*i] for i in range(n+1)])
    el=[[1]+[0]*n]
    for k in range(1,p+1):
        val=[0]*(n+1)
        for j in range(1,k+1):val=add(val,scale(mul(el[k-j],traces[j],n),(-1)**(j-1)))
        el.append(exact_div(val,k))
    first=[th[k//p] if k%p==0 else 0 for k in range(n+1)]
    inverse=inv(th,n);inverse_powers=[[1]+[0]*n]
    for k in range(1,p+2):inverse_powers.append(mul(inverse_powers[-1],inverse,n))
    polys=[]
    for k in range(p+2):
        total=scale(el[k],p**(2*(p-k))) if k<=p else [0]*(n+1)
        if k:total=add(total,scale(mul(first,el[k-1],n),p**(2*(p-k+1))))
        series=scale(mul(total,inverse_powers[k],n),(-1)**k)
        rem=series[:];poly=[];tp=[1]+[0]*n
        for d in range(k//2+1):
            v=rem[d];poly.append(v);rem=add(rem,scale(tp,-v));tp=mul(tp,t,n)
        need(not any(rem),f'Hecke reconstruction p={p} form={name} coefficient={k}')
        polys.append(poly)
        if k<p:need(all(v%(p*p)==0 for v in poly),'Hecke low divisibility')
    need(polys[-1][0]%p==1 and all(v%p==0 for v in polys[-1][1:]),'leading modp')
    g=central(sequence(*PARAMS[name],p))
    finite=add(polys[p],mul(polys[p+1],g[:p],len(polys[p+1])+p-2))
    need(not any(x%(p*p) for x in finite),'finite polynomial identity')
    t0=F(1,20) if name=='a' else F(-1,16)
    specialized=[ev(poly,t0) for poly in polys]
    row={'p':p,'form':name,'coefficients':polys,'specialized':[str(x) for x in specialized],'t0':str(t0),'q_terms':n+1,'integer_coefficients':sum(map(len,polys))}
    if leg(-5,p)==1:
        m=1 if p%4==1 else 2
        xy=next(( (x,y) for y in range(1,isqrt(m*p//5)+1) for x in [isqrt(m*p-5*y*y)] if x>0 and x*x+5*y*y==m*p),None)
        need(xy is not None,'CM representation')
        x,y=xy;trace=4*x*x-2*p if m==1 else 2*p-2*x*x
        rem=quadratic_remainder(specialized,trace,p*p)
        bad=quadratic_remainder(specialized,-trace,p*p)
        need(rem==[0,0],'exact CM factor');need(any(bad),'wrong-sign CM factor must fail')
        row.update({'representation':xy,'representation_multiplier':m,'raw_trace':trace,'factor':[p*p,-trace,1],'remainder':[str(x) for x in rem],'wrong_sign_remainder':[str(x) for x in bad]})
    return row


def run(out):
    summary=certificates(out);hp=[]
    for p in (7,11,29):
        need(p in primes(p) and p>5,'Hecke prime')
        for name in ('a','Q'):
            hp.append(hecke(p,name))
    save(out/'hecke_polynomials.json',hp)
    summary.update({'hecke_polynomials':len(hp),'hecke_integer_coefficients':sum(r['integer_coefficients'] for r in hp),'exact_CM_factors':sum('factor' in r for r in hp),'ok':True})
    save(out/'modular_summary.json',summary)
    return summary

