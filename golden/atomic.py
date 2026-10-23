"""CHIANTI 11.0.2 C IV 1548+1551 A: coronal excitation and radiative cascades.

Low-density, Maxwellian, ionization-equilibrium approximation. This does not
claim non-equilibrium C populations or optically thick chromospheric transfer.
Uses all supplied ground-excitation transitions and radiative decay branches,
not a hand-chosen Gaussian temperature response.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.interpolate import splrep,splev
HERE=Path(__file__).resolve().parent;DATA=HERE/'data/atomic'

def coefficient():
    out=DATA/'CIV-coefficient.npz'
    if out.exists():
        d=np.load(out);return d['logT'],d['G_erg_cm3_s_sr']
    lines=(DATA/'c_4.elvlc').read_text().splitlines();levels={}
    for s in lines:
        if s.strip()=='-1':break
        # CHIANTI fixed-width energy-level format; first comma suffix excluded.
        a=s.split(',')[0].split();levels[int(a[0])]=(float(a[-2]) if float(a[-2])>=0 else float(a[-1]),2*float(a[-3])+1)
    N=max(levels);energy=np.zeros(N);weight=np.zeros(N)
    for k,(e,g) in levels.items():energy[k-1]=e;weight[k-1]=g
    A=np.zeros((N,N));wavelength={}
    for s in (DATA/'c_4.wgfa').read_text().splitlines():
        if s.strip()=='-1':break
        a=s.split();lo,up=int(a[0])-1,int(a[1])-1
        if up<N and lo<N:A[up,lo]+=float(a[4])
        wavelength[(lo,up)]=abs(float(a[2]))
    lt=np.linspace(4,7,1201);T=10**lt;exc=np.zeros((N,len(T)));s=(DATA/'c_4.scups').read_text().splitlines();k=0
    while k<len(s) and s[k].strip()!='-1':
        a=s[k].split();lo,up=int(a[0])-1,int(a[1])-1;de=float(a[2]);typ=int(a[6]);C=float(a[7]);xs=np.fromstring(s[k+1],sep=' ');ys=np.fromstring(s[k+2],sep=' ');k+=3
        if lo!=0:continue
        kte=8.617333262145e-5*T/(de*13.605693122994)
        st=1-np.log(C)/np.log(kte+C) if typ in [1,4] else kte/(kte+C)
        su=splev(st,splrep(xs,ys,s=0))
        if typ==1:ups=su*np.log(kte+np.e)
        elif typ==2:ups=su
        elif typ==3:ups=su/(kte+1)
        elif typ==4:ups=su*np.log(kte+C)
        elif typ==5:ups=su/kte
        elif typ==6:ups=10**su
        else:raise ValueError(typ)
        # Observed energy difference in Boltzmann factor, matching ChiantiPy.
        exc[up]=8.629132e-6*np.maximum(0,ups)/weight[0]/np.sqrt(T)*np.exp(-(energy[up]-energy[0])*1.4387768775039338/T)
    direct=exc.copy();photon=np.zeros_like(exc);arrival=exc.copy()
    for up in np.argsort(energy)[::-1]:
        if up==0:continue
        total=A[up].sum()
        if total<=0:continue
        pop=arrival[up]/total;photon[up]=pop*A[up,0]
        targets=np.flatnonzero(A[up]);arrival[targets]+=A[up,targets,None]*pop[None,:]
    io=(DATA/'chianti.ioneq').read_text().splitlines();temps=np.fromstring(io[1],sep=' ')
    for s in io[2:]:
        a=s.split()
        if len(a)>3 and a[:2]==['6','4']:
            fraction=np.asarray(a[2:],float);break
    fi=np.interp(lt,temps,fraction,left=0,right=0)
    G=np.zeros_like(T);Gdirect=np.zeros_like(T)
    for up in [1,2]:
        erg=1.9864458571489286e-8/wavelength[(0,up)]
        G+=photon[up]*erg;Gdirect+=direct[up]*erg
    # photospheric carbon abundance 12+log10(nC/nH)=8.43; H/ne=1/1.2
    G*=fi*10**(8.43-12)/1.2/(4*np.pi)
    np.savez(out,logT=lt,G_erg_cm3_s_sr=G)
    rec=dict(CHIANTI_version=(DATA/'VERSION').read_text().strip(),levels=N,formula='j = ne^2 G(T), C IV 1548.187 + 1550.772 A; low-density Maxwellian ground excitation plus radiative cascades',assumptions=['Equilibrium ion fractions from chianti.ioneq','nH/ne=1/1.2 in the C IV emitting phase','Carbon abundance 8.43 on logH=12 scale','No time-dependent C ionization, resonance scattering, or line absorption'],peak_logT=float(lt[G.argmax()]),peak_G=float(G.max()),cascade_multiplier_at_peak=float((Gdirect*0+G/G.max())[G.argmax()]) if False else float((sum(photon[u,G.argmax()]*1.9864458571489286e-8/wavelength[(0,u)] for u in [1,2]))/Gdirect[G.argmax()]),files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.suffix in ['.scups','.elvlc','.wgfa','.ioneq']})
    (DATA/'CIV-coefficient.json').write_text(json.dumps(rec,indent=2)+'\n');return lt,G

if __name__=='__main__':
    l,g=coefficient();print((DATA/'CIV-coefficient.json').read_text())
