#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use FindBin;
use YAML::Syck qw();

use List::Util qw();

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

=head1 EXAMPLE

    $ echo -e 'GCATGCG\nGATTACA' |
        perl dd-dp.pl --sep ""

    $ echo -e 'AAA\tFer2\tZZ\nFer2\tZZ' |
        perl dd-dp.pl

    $ echo -e 'a b c d e f g h\nb c h' |
        perl dd-dp.pl --sep " "

=head1 DESCRIPTION

The Needleman--Wunsch Algorithm

Reference:

L<https://pubmed.ncbi.nlm.nih.gov/5420325/>

L<https://bioboot.github.io/bimm143_W20/class-material/nw/>

L<https://gist.github.com/bibymaths/b1d649fe0ed6e641bf3948cc4d36ebe9>

L<https://metacpan.org/pod/Algorithm::NeedlemanWunsch>

=cut

Getopt::Long::GetOptions(
    'help|?'  => sub { Getopt::Long::HelpMessage(0) },
    'ma=f'    => \( my $opt_match         = 1.0 ),
    'mm=f'    => \( my $opt_mismatch      = 0.0 ),
    'go=f'    => \( my $opt_gap_opening   = -0.01 ),
    'ge=f'    => \( my $opt_gap_extension = -0.001 ),
    'sep|s=s' => \( my $opt_separator     = "\t" ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#

# Global Variables
# The alignment matrix
my @mat;

{
    chomp( my $l = <> );
    my $s = [ split( $opt_separator, $l ) ];
    chomp( $l = <> );
    my $t = [ split( $opt_separator, $l ) ];

    # print YAML::Syck::Dump([$s, $t]);

    # fill
    similarity( $s, $t );

    for my $x ( align( $s, $t ) ) {
        print join( $opt_separator, @{$x} ), "\n";
    }
}

#----------------------------------------------------------#
# Subs
#----------------------------------------------------------#
sub compare {

    # domains/bases/residues
    my $c1 = shift;
    my $c2 = shift;

    return ( $c1 eq $c2 ) ? $opt_match : $opt_mismatch;
}

# filling in the alignment matrix
sub similarity {

    # array ref
    my $s = shift;
    my $t = shift;

    for my $i ( 0 .. scalar @{$s} ) { $mat[$i][0] = $opt_gap_opening * $i; }
    for my $j ( 0 .. scalar @{$t} ) { $mat[0][$j] = $opt_gap_opening * $j; }

    for my $i ( 1 .. scalar @{$s} ) {
        for my $j ( 1 .. scalar @{$t} ) {
            my $p = compare( $s->[ $i - 1 ], $t->[ $j - 1 ] );

            $mat[$i][$j] = List::Util::max(
                0,
                $mat[ $i - 1 ][$j] + $opt_gap_opening,
                $mat[$i][ $j - 1 ] + $opt_gap_opening,
                $mat[ $i - 1 ][ $j - 1 ] + $p
            );
        }
    }

    # Returns nothing but fills @mat.
    return ( $mat[ scalar @{$s} ][ scalar @{$t} ] );
}

# Recursively reconstructs best alignment
# Returns two array refs
sub align {

    # array ref
    my $s = shift;
    my $t = shift;

    my ( $i, $j ) = ( scalar @{$s}, scalar @{$t} );
    return ( [ ("-") x $j ], $t )             if ( $i == 0 );
    return ( $s,             [ ("-") x $i ] ) if ( $j == 0 );
    my ( $s_last, $t_last ) = ( $s->[-1], $t->[-1] );

    if (
        $mat[$i][$j] == $mat[ $i - 1 ][ $j - 1 ] + compare( $s_last, $t_last ) )
    {
        # Case 1: last elements are paired in the best alignment
        my ( $sa, $ta ) =
          align( [ @{$s}[ 0 .. $i - 2 ] ], [ @{$t}[ 0 .. $j - 2 ] ] );
        push @{$sa}, $s_last;
        push @{$ta}, $t_last;
        return ( $sa, $ta );
    }
    elsif ( $mat[$i][$j] == $mat[ $i - 1 ][$j] + $opt_gap_opening ) {

        # Case 2: last element of the 1st array is paired with a gap
        my ( $sa, $ta ) = align( [ @{$s}[ 0 .. $i - 2 ] ], $t );
        push @{$sa}, $s_last;
        push @{$ta}, "-";
        return ( $sa, $ta );
    }
    else {
        # Case 3: last element of the 2nd array is paired with a gap
        my ( $sa, $ta ) = align( $s, [ @{$t}[ 0 .. $j - 2 ] ] );
        push @{$sa}, "-";
        push @{$ta}, $t_last;
        return ( $sa, $ta );
    }
}
