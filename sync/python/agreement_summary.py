# Summarises the agreement between the groups from the summary.json files written by compare_run3_sync_groups.py
# (usage: python3 sync/python/agreement_summary.py > tables.md).  Samples: data (Data_Run2024H) and MC (ggH, VBFH), reported separately.
import json,glob,os,statistics as st,datetime,collections
O=os.environ.get("OUT") or os.path.expandvars("$HMUMU_WWW")+"/skim_v4/studies/sync/groups_comparison_13Oct"
MINT=float(os.environ.get("MINT","17.4"))   # ignore summaries written before this hour of the day (stale runs)
CH=["ggH","VBF","inclusive","VBF_onlyLeading","ggH_noVBF_onlyLeading"]
KIND={"Data_Run2024H":"Data","ggH":"MC","VBFH":"MC"}
MODES=["unweighted","all_common_weights","genWeights"]
EXCL={"mu1_charge","mu2_charge"}
rows=collections.defaultdict(list); vard=collections.defaultdict(lambda: collections.defaultdict(list)); pairs=collections.defaultdict(set); stale=[]
for pd in sorted(glob.glob(O+"/*-*/")):
    name=os.path.basename(pd.rstrip("/")); fs=glob.glob(pd+"*/summary.json")
    if not fs: continue
    if not all((lambda t:t.hour+t.minute/60>=MINT)(datetime.datetime.fromtimestamp(os.path.getmtime(f))) for f in fs): stale.append(name); continue
    for f in fs:
        d=json.load(open(f)); g1,g2=d["groups"]; s=d["sample"]
        if s not in KIND: continue
        k=KIND[s]; pairs[name].add(s)
        for reg,rc in d["cells"].items():
            if s.startswith("Data") and reg=="Signal_Fit": continue
            for ch,c in rc.items():
                a,b=c["selected"][g1],c["selected"][g2]
                if min(a,b)<50: continue
                rows[(name,k)].append(dict(s=s,reg=reg,ch=ch,A=c["both"]/min(a,b),R=a/b))
        for reg,rc in d["distributions"].items():
            if s.startswith("Data") and reg=="Signal_Fit": continue
            for ch,cc in rc.items():
                m=next((x for x in MODES if x in cc),None)
                if not m: continue
                for v,x in cc[m]["variables"].items():
                    if v not in EXCL and min(x["entries"].values())>=200: vard[(name,k)][v].append(x["shape_distance"])
mean=lambda l: sum(l)/len(l) if l else float("nan")
fmt=lambda l: f"{100*mean([x['A'] for x in l]):.1f}% / {100*mean([abs(x['R']-1) for x in l]):.0f}%" if l else "-"
out=["Pairs used: "+", ".join(sorted(pairs))+".  Not yet regenerated (excluded): "+(", ".join(stale) or "none")+"."]
for k,title in (("Data","DATA (`Data_Run2024H`; Z and H sidebands, Signal_Fit blinded)"),("MC","MC (`ggH` and `VBFH` samples; Z sideband, H sideband, Signal_Fit)")):
    out.append(f"\n# {title}\n")
    prs=sorted({p for p,kk in rows if kk==k},key=lambda p:-mean([x["A"] for x in rows[(p,k)]]))
    tot={p:(mean([x["A"] for x in rows[(p,k)]]),mean([abs(x["R"]-1) for x in rows[(p,k)]])) for p in prs}
    out.append("## (1) Pairs ranked from best to worst agreement, per channel\n\nCells: agreement (common / smaller selection) / mean |yield ratio - 1|.\n")
    out.append("| rank | pair | "+" | ".join(CH)+" | all channels |"); out.append("|---:|---|"+"---|"*(len(CH)+1))
    for i,p in enumerate(prs,1):
        out.append(f"| {i} | {p} | "+" | ".join(fmt([x for x in rows[(p,k)] if x["ch"]==ch]) for ch in CH)+f" | **{100*tot[p][0]:.1f}%** / {100*tot[p][1]:.0f}% |")
    if k=="MC":
        out.append("\nSplit by sample (all channels and regions):\n"); out.append("| pair | ggH | VBFH |"); out.append("|---|---|---|")
        for p in prs: out.append(f"| {p} | "+" | ".join(fmt([x for x in rows[(p,k)] if x["s"]==s]) for s in ("ggH","VBFH"))+" |")
    for ch in ("ggH","VBF"):
        r=sorted(((mean([x["A"] for x in rows[(p,k)] if x["ch"]==ch]),p) for p in prs),reverse=True)
        out.append(f"\nChannel `{ch}`, best to worst: "+", ".join(f"{p} {100*a:.1f}%" for a,p in r))
    out.append("\n## (2) Group by group: partners ranked from best to worst\n\nFor each partner: agreement, mean |yield ratio - 1|, mean shape distance of the variables (0 identical).\n")
    vs={p:mean([st.median(v) for v in vard[(p,k)].values()]) for p in prs}
    for g in sorted({g for p in prs for g in p.split("-")}):
        l=sorted(((p,tot[p],vs[p]) for p in prs if g in p.split("-")),key=lambda t:-t[1][0])
        out.append(f"- **{g}**: "+"; ".join(f"{[x for x in p.split('-') if x!=g][0]} {100*a[0]:.1f}% (|R-1| {100*a[1]:.0f}%, dist {d:.3f})" for p,a,d in l))
    out.append("\n## (3) Variables per pair, best to worst\n\nMedian shape distance over the cells with at least 200 entries in both groups (muon charges excluded). Pairs ordered as in (1).\n")
    out.append("| pair | 5 best variables | 5 worst variables |"); out.append("|---|---|---|")
    for p in prs:
        r=sorted((st.median(v),n) for n,v in vard[(p,k)].items())
        out.append(f"| {p} | "+", ".join(f"`{n}` {m:.3f}" for m,n in r[:5])+" | "+", ".join(f"`{n}` {m:.3f}" for m,n in r[::-1][:5])+" |")
print("\n".join(out))
