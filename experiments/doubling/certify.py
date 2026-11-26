import json
import math
import subprocess
import tempfile
from pathlib import Path
from experiment import HERE,exact,fmpz
from eisenstein import local_form
from structural import pure_e4

def counterexample(p,q):
    N=p*q
    checks={}
    for s,t in ((p,q),(q,p)):
        assert fmpz(s).is_prime()
        h=(t-1)%(s-1)
        a,b,c,d=local_form(h)
        assert a%s==0 and b%s and c%s and d%s
        checks[s]=dict(weight=h,hasse_exponent=(N-1-h)//(s-1),
                       G_scaled=a*pow(b,-1,s)%s,constant=c*pow(d,-1,s)%s)
    B=exact(N)[0]
    assert math.gcd(B,N)==1
    with tempfile.TemporaryDirectory() as tmp:
        binary=Path(tmp)/'direct'
        subprocess.run(['c++','-O3','-std=c++17',str(HERE/'direct.cpp'),'-o',str(binary)],check=True)
        out=subprocess.check_output([str(binary),str(N)],text=True)
        assert list(map(int,out.split()))==[N,B]
    return dict(N=N,B=B,gcd=1,locals=checks)

if __name__=='__main__':
    pole=json.loads((HERE/'eisenstein.json').read_text())['pole']
    a,b=pure_e4(100)
    hasse=a*pow(b,-1,101)%101
    pole.update(H100=hasse,G_scaled=pole['G68_CM']*pow(hasse,271,101)%101)
    result=dict(pole=pole,false_positive=counterexample(397,2755601))
    (HERE/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)
