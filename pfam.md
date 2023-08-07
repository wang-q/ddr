# PFAM-A domains

## Download

```shell
mkdir -p $HOME/data/ddr/pfam35
cd $HOME/data/ddr/pfam35

curl -LO https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/relnotes.txt
curl -LO https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/userman.txt

# hmm
aria2c -x 4 -s 2 -c https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz

# hmm data
aria2c -x 4 -s 2 -c https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.dat.gz

# seed
aria2c -x 4 -s 2 -c https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.seed.gz

```

## Descriptive stats

### tags

* #=GF <feature> <Generic per-File annotation, free text>
* ID Identification:      One word name for family
* AC Accession number:    Accession number in form PFxxxxx (Pfam) or RFxxxxx (Rfam)
* DE Definition:          Short description of family
* GA Gathering threshold: Search threshold to build the full alignment
* TP Type:                Type of family -- presently Family, Domain, Motif or Repeat for Pfam
* ML:                     (Domain length)
* CL Clan:                Clan accession
* NE Pfam accession:      Indicates a nested domain.

```shell
gzip -dcf Pfam-A.hmm.dat.gz |
    grep "^#=" |
    tsv-summarize -d" " -g 1,2 --count
#=GF ID 19632
#=GF AC 19632
#=GF DE 19632
#=GF GA 19632
#=GF TP 19632
#=GF ML 19632
#=GF CL 7769
#=GF NE 159

```

### type

```shell
gzip -dcf Pfam-A.hmm.dat.gz |
    grep "#=GF TP" |
    perl -nlp -e 's/#=GF\s+TP\s+//g' |
    tsv-summarize -g 1 --count |
    tsv-sort -k2,2nr |
    (echo -e "item\tcount" && cat) |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'


```

| item        | count |
|-------------|------:|
| Family      | 11710 |
| Domain      |  6909 |
| Repeat      |   624 |
| Coiled-coil |   177 |
| Motif       |   110 |
| Disordered  |   102 |

### clan

```shell
gzip -dcf Pfam-A.hmm.dat.gz |
    grep "#=GF CL" |
    perl -nlp -e 's/#=GF\s+CL\s+//g' |
    tsv-summarize -g 1 --count |
    tsv-summarize --quantile 2:0,0.25,0.5,0.75,1 |
    (echo "#quantile" && cat)
##quantile
#2       3       5       10      381

```

## Mapping files

```shell
cd $HOME/data/ddr/pfam35

# All useful fields
gzip -dcf Pfam-A.hmm.dat.gz |
    perl -nl -e '
        BEGIN { our %info = (); }

        if ( m(#=GF\s+(\w+)\s+(.+)$) ) {
            my $k = $1;
            my $v = $2;
            $info{$k} = $v;
        }
        if ($_ eq q(//)) {
            my @fields = ();
            for $k (qw(ID AC TP CL ML DE)) {
                push @fields, $info{$k};
            }
            print join qq(\t), @fields;
            %info = ();
        }
    ' \
    > fields.tsv

tsv-select -f 1,2 fields.tsv > ID-AC.tsv
tsv-select -f 2,1 fields.tsv > AC-ID.tsv

tsv-select -f 1,4 fields.tsv |
    tsv-filter --not-empty 2 \
    > ID-CL.tsv

# idx  name                 accession        nseq eff_nseq      M relent   info p relE compKL
#hmmstat Pfam-A.hmm.gz |
#    grep -v "^#" |
#    tr -s ' ' |
#    cut -d ' ' -f2,4,6

cat fields.tsv |
    tsv-filter --lt 5:20

cat fields.tsv |
    tsv-filter --gt 5:1000

cat fields.tsv |
    tsv-filter --char-len-gt 1:14

```

## Similarity between two domains

Goodness of seed sequences--profile alignments

### Extract HMM profiles and seed sequences

```shell
mkdir -p $HOME/data/ddr/pfam/profile
mkdir -p $HOME/data/ddr/pfam/fas
mkdir -p $HOME/data/ddr/pfam/seed

cd $HOME/data/ddr/pfam

cp ../pfam35/fields.tsv .

if [[ ! -f Pfam-A.hmm ]]; then
    gzip -dcf ../pfam35/Pfam-A.hmm.gz > Pfam-A.hmm
fi
hmmfetch --index Pfam-A.hmm

# hmmfetch can't be used in this seed file
if [[ ! -f Pfam-A.seed ]]; then
    gzip -dcf ../pfam35/Pfam-A.seed.gz > Pfam-A.seed
fi
esl-afetch --index Pfam-A.seed

cat fields.tsv |
tsv-select -f 1 |
parallel --no-run-if-empty --linebuffer -k -j 8 '
    if [ $(({#} % 100)) -eq "0" ]; then
        >&2 printf "."
    fi

    # a2m files created by esl are not well aligned
    # afa uses . as gaps
    esl-afetch --outformat afa Pfam-A.seed {} |
        perl -nl -e '\''
            $_ !~ /^>/ and s/\./-/g;
            print;
        '\'' \
        > fas/{}.fas

    # remove dashes
    faops filter -d -l 0 fas/{}.fas seed/{}.fasta

    # also as a list of names
    faops size seed/{}.fasta > seed/{}.sizes

    # the hmm profile
    hmmfetch Pfam-A.hmm {} > profile/{}.hmm
'

```

### Example of `hmmalign`

* The `#=GR PP` annotation
    * The posterior probability is encoded as 11 possible characters 0-9*+: 0.0 <= p < 0.05 is coded
      as 0, 0.05 <= p < 0.15 is coded as 1, (... and so on ...), 0.85 <= p < 0.95 is coded as 9, and
      0.95 <= p <= 1.0 is coded as ’*’. Gap characters appear in the PP line where no residue has
      been assigned.

```shell
cd $HOME/data/ddr/pfam

hmmalign --amino --trim profile/GATA.hmm seed/ZZ.fasta

hmmalign --amino --trim profile/ZZ.hmm seed/Prion_octapep.fasta

cat <<EOF > sample.domain.lst
AAA
ABC_tran
ABC_trans_N
Arena_RNA_pol
Fer2
GATA
Prion_octapep
ZZ
EOF

```

```text
seed: ZZ    profile: GATA
# STOCKHOLM 1.0

RSC8_YEAST/254-298           CHTCGNES-INVRYHnlRARDTNLCSRC----------
#=GR RSC8_YEAST/254-298   PP 55555433.4444333433334488888..........
YOY6_CAEEL/8-52              CFSCFTTK-------..---------------------
#=GR YOY6_CAEEL/8-52      PP 88899886..............................
DMD_MOUSE/3300-3345          ---------------..KHFNYDICQSC----------
#=GR DMD_MOUSE/3300-3345  PP .................33445555555..........
REF2P_DROME/121-165          ---------------..---NYDLCQKCELAHKH----
#=GR REF2P_DROME/121-165  PP ....................678*****999997....
ADA2_YEAST/1-46              ---------------..----YDLCVPC----------
#=GR ADA2_YEAST/1-46      PP .....................3466666..........
CBP1_CAEEL/1493-1532         ---------------..-------CNKC----------
#=GR CBP1_CAEEL/1493-1532 PP ........................5555..........
CBP_MOUSE/1702-1743          CINCYNTK-------..---------------------
#=GR CBP_MOUSE/1702-1743  PP 55555555..............................
#=GC PP_cons                 66666655.444433..33355577777999997....
#=GC RF                      xxxxxxxxxxxxxxx..xxxxxxxxxxxxxxxxxxxxx
//
```

### Compute

* 为了减少随机因素的影响，在此使用中位数

```shell
cd $HOME/data/ddr/pfam

cp ~/Scripts/ddr/bin/ss-p.sh .

# bash ss-p.sh ZZ GATA

cat sample.domain.lst |
while read S; do `# seed`
    cat fields.tsv |
    tsv-select -f 1 |
    while read P; do `# profile`
        echo -e "${S}=${P}"
    done
done \
    > sample.job.lst

touch s2d.tsv

cat sample.job.lst |
tsv-join -f s2d.tsv -k 1 -e |
parallel --colsep "=" --no-run-if-empty --linebuffer -k -j 8 '
    if [ $(({#} % 100)) -eq "0" ]; then
        >&2 printf "."
    fi
    bash ss-p.sh {1} {2}
    ' \
    > s2d.tmp

#while IFS='=' read -r S P; do
#    echo >&2 -e "seed: ${S}\tprofile: ${P}"
#
#    bash ss-p.sh ${S} ${P}
#done

# Combine new results with the old ones
cat s2d.tsv s2d.tmp |
    tsv-uniq |
    sort \
    > tmp.tsv
mv tmp.tsv s2d.tsv
rm *.tmp

cat s2d.tsv |
    tsv-filter --or --gt 2:0 --gt 3:0

```

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/ddr/ \
    wangq@202.119.37.251:data/ddr

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/ddr/ \
    wangq@58.213.64.36:data/ddr

# back
rsync -avP \
    -e 'ssh -p 8804' \
    wangq@58.213.64.36:data/ddr/ \
    ~/data/ddr

```

## Useless

### alfpy

Very fast, but unable to distinguish different domains from the similarity results

```shell
cd $HOME/data/ddr/pfam35

pip install alfpy

calc_bbc.py --fasta domains/ZZ.fasta -m protein --outfmt pairwise
calc_wmetric.py --fasta domains/ZZ.fasta --outfmt pairwise
calc_word.py --fasta domains/ZZ.fasta --word_size 3 --outfmt pairwise

calc_bbc.py --fasta <(cat domains/*.fasta) -m protein --outfmt pairwise |
    grep RSC8_YEAST

calc_wmetric.py --fasta <(cat domains/*.fasta) --outfmt pairwise |
    grep RSC8_YEAST

calc_word.py --fasta <(cat domains/*.fasta) --vector counts --distance google --word_size 1 --outfmt pairwise |
    grep RSC8_YEAST

```

### pseqsid

Unable to specify output file name and output format

```shell
cd $HOME/data/ddr/pfam35

cargo install pseqsid

pseqsid -n domains/ZZ.fas

```
