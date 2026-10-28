import re, sys
root='/Volumes/git/factoring-connection/notes/'
exec(open('scratchpad/audit2L/make_blind.py').read().split('# remove parentheticals')[0].split("root='/Volumes")[0])
parts=eval(open('scratchpad/audit2L/make_blind.py').read().split('parts=')[1].split('\n]')[0]+'\n]')
kw=re.compile(r'(codex|Codex|critic|reviewer|Reviewer|\bround\b|thmAB_\d|c12rat_\d|pc3_\d|uprime_\d|multidescent_\d|c12alg_\d|\bv[123](\.\d)?\b|both reviewers)')
def strip_groups(t):
    out=[];i=0;n=len(t)
    while i<n:
        if t[i]=='(':
            depth=0;j=i
            while j<n and j-i<600:
                if t[j]=='(':depth+=1
                elif t[j]==')':
                    depth-=1
                    if depth==0:break
                j+=1
            if j<n and depth==0 and kw.search(t[i:j+1]):
                # drop group and preceding space
                while out and out[-1]==' ': out.pop()
                i=j+1;continue
        out.append(t[i]);i+=1
    return ''.join(out)
def L(f,a,b): return open(root+f).read().split('\n')[a-1:b]
out=[open('scratchpad/audit2L/header.md').read()]
for f,a,b,title in parts:
    body='\n'.join(L(f,a,b))
    body=strip_groups(body)
    body='\n'.join(l for l in body.split('\n') if not re.search(r'\d+ ?%',l))
    out.append('\n\n---\n\n# '+title+'\n\n'+body)
open(sys.argv[1],'w').write('\n'.join(out))
