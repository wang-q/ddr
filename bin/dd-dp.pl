#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use FindBin;
use YAML::Syck qw();

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

dd-dp.pl - Domain distance via dynamic programming

=head1 SYNOPSIS

    cat <file> | perl dd-dp.pl [options]
      Options:
        --help                  Help message
        --ma            FLT     Match, default is 1.0
        --mm            FLT     Mismatch, default is 0.0
        --go            FLT     Gap opening, default is -0.01
        --ge            FLT     Gap extension, default is -0.001
        --spe       -s  STR     Separator of fields, default is "\t"
        --global    -g          Global DP, i.e., Needleman--Wunsch

=head1 EXAMPLE

    $ echo -e 'GCATGCU\nGATTACA' |
        perl dd-dp.pl --spe ""

    $ echo -e 'AAA\tFer2\tZZ\nFer2\tZZ' |
        perl dd-dp.pl

=cut

Getopt::Long::GetOptions(
    'help|?'   => sub { Getopt::Long::HelpMessage(0) },
    'ma=f'     => \( my $match         = 1.0 ),
    'mm=f'     => \( my $mismatch      = 0.0 ),
    'go=f'     => \( my $gap_opening   = -0.01 ),
    'ge=f'     => \( my $gap_extension = -0.001 ),
    'sep|s=s' => \( my $separator   = "\t" ),
    'global|g' => \( my $is_global ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
