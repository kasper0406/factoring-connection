import re, sys
root='/Volumes/git/factoring-connection/notes/'
def lines(f,a,b):
    L=open(root+f).read().split('\n')
    return L[a-1:b]
parts=[
 ('thmAB_results.md',85,986,'A. Theorem AB/AB2 and the assembly (source: notes/thmAB_results.md §§1-14)'),
 ('c12rat_results.md',91,760,'B. Theorem RAT (source: notes/c12rat_results.md §§1-8)'),
 ('c12rat_results.md',805,812,'B′. RAT precision margins (source: notes/c12rat_results.md §9.2)'),
 ('pc3_results.md',105,204,'C. Prop DC / DC_e, Lemma 0, Lemma LAT (source: notes/pc3_results.md §§1-3)'),
 ('pc3_results.md',288,331,'C′. Design DW, Lemma DW(a), DW-L (source: notes/pc3_results.md §6.2)'),
 ('multidescent_results.md',234,346,'D. Size model, Prop UB, Theorem 2L (source: notes/multidescent_results.md §§3.1-3.3)'),
 ('uprime_results.md',96,375,'E. Hyp U′: Lemma 8.5H, level-1 moments, CS′, SQF, recursion (source: notes/uprime_results.md §§1-6)'),
 ('c12alg_results.md',133,531,'F. Window design, Lemma η, Prop PB, Lemma VD, DW lemmas, Lemma TP (source: notes/c12alg_results.md §§1-3.2)'),
 ('c12alg_results.md',872,1029,'G. Prop CLB and the conditional assembly 2L-DW (source: notes/c12alg_results.md §§4-5)'),
]
# remove parentheticals mentioning reviews/attributions
kw=r'(codex|Codex|critic|reviewer|Reviewer|round|thmAB_\d|c12rat_\d|pc3_\d|uprime_\d|multidescent_\d|c12alg_\d|\bv[123](\.\d)?\b|both reviewers|proof critic)'
def strip_paren(s):
    # iteratively remove innermost (...) containing keyword
    changed=True
    while changed:
        changed=False
        for m in re.finditer(r'\(([^()]*)\)',s):
            if re.search(kw,m.group(1)):
                s=s[:m.start()].rstrip()+s[m.end():]
                changed=True
                break
    return s
out=[]
out.append(open('/Volumes/git/factoring-connection/scratchpad/audit2L/header.md').read())
for f,a,b,title in parts:
    out.append('\n\n---\n\n# '+title+'\n')
    for l in lines(f,a,b):
        l2=strip_paren(l)
        if re.search(r'\d+ ?%',l2): continue
        out.append(l2)
txt='\n'.join(out)
open(sys.argv[1],'w').write(txt)
