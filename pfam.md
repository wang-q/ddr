# PFAM

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
* ML
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
            for $k (qw(ID AC TP CL DE)) {
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

```

## Domain similarity

### Extract HMM profiles and seed sequences

```shell
cd $HOME/data/ddr/pfam35

mkdir domains

if [[ ! -f Pfam-A.hmm ]]; then
    gzip -dcf Pfam-A.hmm.gz > Pfam-A.hmm
fi
hmmfetch --index Pfam-A.hmm

# hmmfetch can't be used in this seed file
if [[ ! -f Pfam-A.seed ]]; then
    gzip -dcf Pfam-A.seed.gz > Pfam-A.seed
fi
esl-afetch --index Pfam-A.seed

cat <<EOF > sample.domain.lst
AAA
Fer2
GATA
ZZ
EOF

cat sample.domain.lst |
while read D; do
    # a2m files created by esl are not well aligned
    # afa uses . as gaps
    esl-afetch --outformat afa Pfam-A.seed ${D} |
        perl -nl -e '
            $_ !~ /^>/ and s/\./-/g;
            print;
        ' \
        > domains/${D}.fas
    faops filter -d -l 0 domains/${D}.fas domains/${D}.fasta

    # as a list of names
    faops size domains/${D}.fasta > domains/${D}.sizes

    # the hmm profile
    hmmfetch Pfam-A.hmm ${D} > domains/${D}.hmm

done

```

### Compute

* The `#=GR PP` annotation
    * The posterior probability is encoded as 11 possible characters 0-9*+: 0.0 <= p < 0.05 is coded
      as 0, 0.05 <= p < 0.15 is coded as 1, (... and so on ...), 0.85 <= p < 0.95 is coded as 9, and
      0.95 <= p <= 1.0 is coded as ’*’. Gap characters appear in the PP line where no residue has
      been assigned.

```text
seed ZZ domain GATA
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

```shell
cd $HOME/data/ddr/pfam35

cat sample.domain.lst |
while read S; do `# seed`
    cat sample.domain.lst |
    while read D; do `# domain`
        echo >&2 -e "seed: ${S}\tdomain: ${D}"

        hmmalign --amino --trim domains/${D}.hmm domains/${S}.fasta |
            T=5 perl -nl -e '
                if ( m(^#=GR\s+([\w/_-]+)\s+PP\s+(.+)$) ) {
                    my $name = $1;
                    my $pp = $2;
                    $pp =~ s/[0-$ENV{T}.]//g; # Removal of unqualified PP characters

                    print join qq(\t), $name, length $pp;
                }
            ' |
            tsv-join -f domains/${S}.sizes -k 1 -a 2 |
            perl -nla -F"\t" -e '
                @F == 3 or next;
                print $F[1] / $F[2];
            ' |
            tsv-summarize --median 1 |
            perl -nl -e 'printf qq(%.4f\n), $_' |
            (printf "${S}\t${D}\t" && cat)
    done
done

```

## Useless

### alfpy

Very fast, but unable to distinguish different structural domains from the similarity results

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