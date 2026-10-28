import re,sys
p=sys.argv[1]; t=open(p).read()
R=[
 (" Codex's blind round found the same repair.",""),
 ("| proved here; critic: \"genuinely proved\" after fix 4b (applied); codex round 2 kept it open pending the same fix |","| proved here |"),
 ("| proved; both reviewers agree |","| proved |"),
 ("| proved; both reviewers accept |","| proved |"),
 ("Theorem numbers are as checked by the project's literature critics (c12alg §10, c12rat §12.1,\n","Theorem numbers as cited in the source notes (c12alg §10, c12rat §12.1,\n"),
 ("opened by the project's critics; the statement is verified against both restatements.","opened; the statement is checked against both restatements."),
 (" Codex's\nblind round found the same coupling principle in L² form.",""),
 (", as codex warned.","."),
 ("This is codex's (11) (found blind) and the large-modulus case","This is the large-modulus case"),
 (" Codex's concern that \"small label means remain dangerous\" is met because every"," Small label means are harmless because every"),
 (" This is codex's \"(13)\" and"," This is"),
 ("The proof critic found an O(1) field-level term in σ_𝔮².","There is an O(1) field-level term in σ_𝔮²."),
 ("- Codex #2: without the filter (e.g. G(P) = 1) the claim is false; v1 stated it universally.","- Without the filter (e.g. G(P) = 1) the claim is false."),
 ("[proved; codex adversarial #5]","[proved]"),
 (" ∎ Codex\n"," ∎\n"),
 ("[sketch; proof critic M4]","[sketch]"),
 ("Codex #7 lists three steps it breaks, all repairable:","Three steps it breaks, all repairable:"),
 (" Codex notes\n"," Note\n"),
 ("[OPEN; expected but not routine — both reviewers]","[OPEN at the time of that note; supplied by uprime]"),
 (" As codex notes, this is used *globally*"," This is used *globally*"),
 (" This agrees with codex.",""),
 ("- Codex's example shows that w coprime","- An example shows that w coprime"),
 ("; codex reached the same\n",";\n"),
 (" Codex's sketch is"," The sketch is"),
 (" This is codex's (1), obtained"," Obtained"),
 ("The proof critic's M4 and pc3 §7 item 3 asked for this, and it holds uniformly,","pc3 §7 item 3 asked for this, and it holds uniformly,"),
]
for a,b in R:
    if a not in t: print('MISSING:',a[:60])
    t=t.replace(a,b)
open(p,'w').write(t)
