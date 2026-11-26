import json
import gmpy2
from experiment import HERE, factorial
from structural import block_product

def validate(row):
    p,q,N = row['p'],row['q'],row['N']
    d = (q-p)//4
    a = block_product(N,1,d)
    b = block_product(N,3,d)
    A = a*pow(b,-1,N)%N
    U = row['B']*pow(4,-(N//4),N)%N
    character = 1
    for s,trace_half in ((p,row['ap']),(q,row['aq'])):
        assert U%s == 4*trace_half**2*A % s
        D = 4*d
        odd = factorial(D,s)*pow(factorial(D//2,s),-1,s)*pow(2,-(D//2),s)%s
        character *= int(gmpy2.jacobi(odd,s))
    assert character == int(gmpy2.jacobi(row['B'],N))
    return dict(N=N,p=p,q=q,B=row['B'],U=U,A=A,d=d,J=character)

if __name__ == '__main__':
    out=[]
    for size in (32,48,64):
        with (HERE/f'rows{size}.jsonl').open() as stream:
            for line in stream:
                row=json.loads(line)
                if row['c']==1:
                    out.append(validate(row))
    (HERE/'gap.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Validated',len(out),'example',out[0])
