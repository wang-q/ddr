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
Usage: $0 <INFILE> <COL> [PATTERN]

$ bash z-score.sh s2d.tsv 2 '^ZZ='

"

if [ "$#" -lt 2 ]; then
    echo >&2 "$USAGE"
    exit 1
fi

INFILE=$1
COL=$2
PATTERN=$3

#----------------------------#
# Run
#----------------------------#
shopt -s lastpipe # let `read` modifies variables

if [ -n ${PATTERN} ]; then
    cat s2d.tsv |
        grep "${PATTERN}"
else
    cat s2d.tsv
fi |
    tsv-summarize --min ${COL} --max ${COL} |
    IFS=$'\t' read -r MIN MAX

#echo $MIN $MAX

if [ -n ${PATTERN} ]; then
    cat s2d.tsv |
        grep "${PATTERN}"
else
    cat s2d.tsv
fi |
    COL=${COL} MIN=${MIN} MAX=${MAX} perl -nla -F"\t" -e '
        my $col = $ENV{COL} - 1;
        my $x = $F[$col];
        next if $x == 0;
        my $z = ($x - $ENV{MIN}) / $ENV{MAX};
        print join qq(\t), (@F, sprintf(qq(%.4f), $z));
    '
