import json
import subprocess
import tempfile
from pathlib import Path
from experiment import HERE

rows=[json.loads(l) for l in (HERE/'rows32.jsonl').read_text().splitlines()]
subset={r['N']:r['B'] for r in rows if r['index']<2}
with tempfile.TemporaryDirectory() as temp:
    binary=Path(temp)/'direct'
    subprocess.run(['c++','-O3','-std=c++17',str(HERE/'direct.cpp'),'-o',str(binary)],check=True)
    proc=subprocess.run([str(binary),*map(str,subset)],check=True,capture_output=True,text=True)
    got=dict(tuple(map(int,l.split())) for l in proc.stdout.splitlines())
    assert got==subset
    (HERE/'direct.json').write_text(json.dumps(got,indent=2)+'\n')
    print(got)
