#!/usr/bin/env bash

BASH_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "${BASH_DIR}"

#----------------------------#
# Colors in term
#----------------------------#
GREEN=
RED=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    NC='\033[0m' # No Color
fi

log_warn () {
    echo >&2 -e "${RED}==> $@ <==${NC}"
}

log_info () {
    echo >&2 -e "${GREEN}==> $@${NC}"
}

log_debug () {
    echo >&2 -e "==> $@"
}

export -f log_warn
export -f log_info
export -f log_debug

#----------------------------#
# helper functions
#----------------------------#
set +e

# set stacksize to unlimited
if [[ "$OSTYPE" != "darwin"* ]]; then
    ulimit -s unlimited
fi

signaled () {
    log_warn Interrupted
    exit 1
}
trap signaled TERM QUIT INT

readlinkf () {
    perl -MCwd -l -e 'print Cwd::abs_path shift' "$1";
}

export -f readlinkf

#----------------------------#
# Usage
#----------------------------#
USAGE="
Usage: $0 <SEED> <PROFILE>

$ bash ss-p.sh ZZ GATA

"

if [ "$#" -lt 2 ]; then
    echo >&2 "$USAGE"
    exit 1
fi

SEED=$1
PROFILE=$2

#----------------------------#
# Run
#----------------------------#
P_LEN=$(
    cat fields.tsv |
        tsv-filter --str-eq "1:${PROFILE}" |
        tsv-select -f 5
)

hmmalign --amino --trim domains/${PROFILE}.hmm domains/${SEED}.fasta |
    T=5 perl -nl -e '
        if ( m(^#=GR\s+([\w/_-]+)\s+PP\s+(.+)$) ) {
            my $name = $1;
            my $pp = $2;
            $pp =~ s/[0-$ENV{T}.]//g; # Removal of unqualified PP characters

            print join qq(\t), $name, length $pp;
        }
    ' |
    tsv-join -f domains/${SEED}.sizes -k 1 -a 2 |
    P_LEN=${P_LEN} perl -nla -F"\t" -e '
        @F == 3 or next;
        print $F[1] / $F[2], qq(\t), $F[1] / $ENV{P_LEN};
    ' |
    tsv-summarize --median 1 --median 2 |
    perl -nla -F"\t" -e 'printf qq(%.4f\t%.4f\n), $F[0], $F[1]' |
    (printf "${SEED}=${PROFILE}\t" && cat)
