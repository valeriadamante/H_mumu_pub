#!/bin/bash
# Local, sequential: scan per-job tuples, hadd (1 thread), verify entries, copy to the group EOS area.
# usage: merge_sync_update_local.sh <workdir> <final_name> <src_dir> [<src_dir2> ...]
work="$1"; final="$2"; shift 2
SRC=${SRC:-/eos/user/v/vdamante/H_mumu/Sync/update_after_Hmm_15Oct_v3}
DST=/eos/cms/store/group/phys_higgs/cmshmm/run3_sync/TAMU/update_after_Hmm_15Oct
mkdir -p "$work"; cd "$work" || exit 1
cat > scan.C <<'EOC'
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <memory>
#include <string>
#include <cstdio>
void scan(const char* list){
  std::ifstream in(list); std::string f;
  while(in>>f){ std::unique_ptr<TFile> t(TFile::Open(f.c_str()));
    if(!t||t->IsZombie()||t->TestBit(TFile::kRecovered)){printf("BAD %s\n",f.c_str());continue;}
    auto tr=(TTree*)t->Get("Events"); if(!tr){printf("BAD %s\n",f.c_str());continue;}
    printf("OK %lld %s\n",tr->GetEntries(),f.c_str()); }
}
EOC
: > files.txt
for d in "$@"; do ls $SRC/$d/job_*[0-9].root >> files.txt; done
root -l -b -q "scan.C(\"files.txt\")" 2>&1 | grep -E '^(OK|BAD)' > scan.out
echo "files OK: $(grep -c '^OK' scan.out) / $(wc -l < files.txt)"
grep '^BAD' scan.out && { echo "BAD FILES, abort"; exit 1; }
sum=$(awk '/^OK/{s+=$2}END{print s}' scan.out)
hadd -f -j 1 "$final.root" $(cat files.txt) 2>&1 | grep -v 'do not have the same data type' | tail -3 || exit 1
got=$(root -l -b -q -e "TFile f(\"$final.root\");printf(\"%lld\n\",((TTree*)f.Get(\"Events\"))->GetEntries());" | tail -1)
echo "sum=$sum merged=$got"
[ -n "$sum" ] && [ "$sum" = "$got" ] || { echo "ENTRY MISMATCH"; exit 1; }
cp "$final.root" "$DST/$final.root.tmp" && mv "$DST/$final.root.tmp" "$DST/$final.root" && rm -f "$final.root" && echo "DONE $final"
