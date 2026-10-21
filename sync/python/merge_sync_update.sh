#!/bin/bash
# Condor job: scan per-job tuples for truncation, hadd, verify entries, publish to the group EOS area.
# usage: merge_sync_update.sh <final_name> <src_dir> [<src_dir2> ...]
final="$1"; shift
SRC=${SRC:-/eos/user/v/vdamante/H_mumu/Sync/update_after_Hmm_15Oct_v3}
DST=/eos/cms/store/group/phys_higgs/cmshmm/run3_sync/TAMU/update_after_Hmm_15Oct
export X509_USER_PROXY="${X509_USER_PROXY:-}"
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh 2>/dev/null || true
cat > scan.C <<'EOC'
void scan(const char* list){
  std::ifstream in(list); std::string f; long tot=0;
  while(in>>f){ std::unique_ptr<TFile> t(TFile::Open(f.c_str()));
    if(!t||t->IsZombie()||t->TestBit(TFile::kRecovered)){printf("BAD %s\n",f.c_str());continue;}
    auto tr=(TTree*)t->Get("Events"); if(!tr){printf("BAD %s\n",f.c_str());continue;}
    printf("OK %lld %s\n",tr->GetEntries(),f.c_str()); }
}
EOC
: > all.txt
for d in "$@"; do ls $SRC/$d/job_*[0-9].root | sed 's#^#root://eoshome-v.cern.ch/#' >> all.txt; done
root -l -b -q "scan.C+(\"all.txt\")" 2>&1 | grep -E '^(OK|BAD)' > scan.out
grep -c '^OK' scan.out; grep '^BAD' scan.out && { echo "BAD FILES, abort"; exit 1; }
sum=$(awk '/^OK/{s+=$2}END{print s}' scan.out)
sed 's#root://eoshome-v.cern.ch/##' all.txt > files.txt
hadd -f -j 4 "$final.root" $(cat files.txt) || exit 1
got=$(root -l -b -q -e "TFile f(\"$final.root\");printf(\"%lld\n\",((TTree*)f.Get(\"Events\"))->GetEntries());" | tail -1)
echo "sum=$sum merged=$got"
[ "$sum" = "$got" ] || { echo "ENTRY MISMATCH"; exit 1; }
xrdcp -f "$final.root" "root://eoscms.cern.ch/$DST/$final.root" 2>&1
echo DONE
