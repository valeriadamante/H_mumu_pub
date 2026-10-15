# Compares the summary.json of the reproducibility test with the reference numbers (sync/reference/repro_TAMU-Pisa_ggH.json).
# usage: python3 sync/python/check_reproducibility.py <output dir of the test> [--write-reference]
import json,sys
from pathlib import Path
out=Path(sys.argv[1]); ref_path=Path(__file__).resolve().parents[1]/"reference"/"repro_TAMU-Pisa_ggH.json"
d=json.load(open(out/"TAMU-Pisa"/"ggH"/"summary.json"))
cell=d["cells"]["Z_sideband"]["inclusive"]; v=d["distributions"]["Z_sideband"]["inclusive"]["all_common_weights"]["variables"]
now={"entries":d["entries"],"unique_events":d["unique_events"],"common_events":d["common_events"],
     "selected_TAMU":cell["selected"]["TAMU"],"selected_Pisa":cell["selected"]["Pisa"],"both":cell["both"],
     "TAMU_only_absent_from_Pisa":cell["TAMU_only_absent_from_Pisa"],"TAMU_only_present_in_Pisa":cell["TAMU_only_present_in_Pisa"],
     "Pisa_only_absent_from_TAMU":cell["Pisa_only_absent_from_TAMU"],"Pisa_only_present_in_TAMU":cell["Pisa_only_present_in_TAMU"],
     **{f"shape_distance_{n}":round(v[n]["shape_distance"],4) for n in ("m_mumu","mu1_pt","jet1_pt","jet2_pt","njets")}}
if "--write-reference" in sys.argv: json.dump(now,open(ref_path,"w"),indent=1); print("reference written"); sys.exit(0)
ref=json.load(open(ref_path)); bad=0
for k,r in ref.items():
    ok=(now[k]==r) if isinstance(r,list) else abs(now[k]-r)<=(2e-3 if k.startswith("shape") else 0)
    bad+=not ok; print(f"{'OK  ' if ok else 'FAIL'} {k:34s} expected {r}  got {now[k]}")
print("REPRODUCED" if not bad else f"{bad} DIFFERENCES"); sys.exit(1 if bad else 0)
